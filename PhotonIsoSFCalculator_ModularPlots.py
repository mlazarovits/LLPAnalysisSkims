#!/usr/bin/env python3
"""
FullSim/FastSim photon-isolation efficiency-SF framework with BR reweighting.

For each target BR, the available nominal gogoGZ samples are reweighted using
`evtwt * ReweightTool::DecayReweight` logic and the event topology flags
Evt_isZZ, Evt_isGZ, and Evt_isGG. The SF is the ratio of weighted selection
efficiencies relative to the preselection, not a ratio of unit-area histogram
integrals:

    SF = (sum(evtwt*w_decay) in SR / sum(evtwt*w_decay) in preselection)_FullSim
         / (sum(evtwt*w_decay) in SR / sum(evtwt*w_decay) in preselection)_FastSim

All target-BR actions for a sample are booked on the same RDataFrame graph
before materialization.
"""
from __future__ import annotations
import array
import argparse
import csv
import math
import re
import sys
import numpy as np
from collections import defaultdict
from pathlib import Path

import yaml
import ROOT
from ROOT import RDataFrame, TFile, TGraphErrors, TCanvas, TLegend, TLatex
from FileProcessor import FileProcessor


# ------------------------------- User settings -------------------------------

DEFAULT_CONFIG = "BigGuy_NonCompressed_FullRegions_AnalysisConfig.yaml"
DEFAULT_OUTPUT = "photon_iso_sf"

# These are the actual available parameter points in the uploaded FileProcessor.
DEFAULT_CTAUS = None
DEFAULT_PT_EDGES = [0, 100, 150, 200, 300, 500, 1000]

SR_BIN_NAMES = [
    "Ch10CReq1PhoTightIsoPrompt00",
    "Ch10SReq1PhoTightIsoPrompt10",
    "Ch10SReq1PhoTightIsoPrompt01",
    "Ch10SReq1PhoTightIsoPrompt02",
    "Ch10SReq1PhoTightIsoPrompt12",
    "Ch10SReq1PhoTightIsoPrompt11",
    "Ch12CReq2PhoTightIsoPrompt00",
    "Ch12SReq2PhoTightIsoPrompt10",
    "Ch12SReq2PhoTightIsoPrompt01",
    "Ch12SReq2PhoTightIsoPrompt11",
]
DEFAULT_BR_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0]

# C++ port of the supplied ReweightTool::DecayReweight. The ROOT binomial_pdf
# calls used by the original function evaluate to p^2 for k=n=2.
ROOT.gInterpreter.Declare(r"""
inline double PhotonIsoDecayReweight(bool isZZ, bool isGZ, bool isGG,
                                     double Z_old, double Z_new,
                                     double G_old, double G_new) {
    const double pZZ0 = Z_old * Z_old;
    const double pGG0 = G_old * G_old;
    const double pGZ0 = 1.0 - (pZZ0 + pGG0);
    double pZZ1 = (Z_new == 0.0 || Z_new == 1.0) ? Z_new : Z_new * Z_new;
    double pGG1 = (G_new == 0.0 || G_new == 1.0) ? G_new : G_new * G_new;
    const double pGZ1 = 1.0 - (pZZ1 + pGG1);
    double decaywt = 0.0;
    if (isZZ) decaywt = pZZ1 / pZZ0;
    if (isGZ) decaywt = pGZ1 / pGZ0;
    if (isGG) decaywt = pGG1 / pGG0;
    return decaywt;
}
""")


from ROOT import TLatex

def draw_fit_equation(canvas, nom_fit, up_fit, dn_fit,
                     x=0.15, y=0.85, text_size=0.035):

    # Extract parameters: [0] = intercept, [1] = slope
    nom_int, nom_slope = nom_fit[0], nom_fit[1]
    up_int,  up_slope  = up_fit[0],  up_fit[1]
    dn_int,  dn_slope  = dn_fit[0],  dn_fit[1]

    print(nom_int, nom_slope)
    print(up_int, dn_int)
    # Build TLatex expression
    latex_str = (
        f"Overall fit: "
        f"{nom_slope:.3g}^{{+{up_slope - nom_slope:.3g}}}_{{{dn_slope - nom_slope:.3g}}}x"
        f" + "
        f"{nom_int:.3g}^{{+{up_int - nom_int:.3g}}}_{{{dn_int - nom_int:.3g}}}"
    )

    canvas.cd()

    latex = TLatex()
    latex.SetNDC(True)  # Coordinates relative to canvas
    latex.SetTextSize(text_size)
    latex.DrawLatex(x, y, latex_str)

    canvas.Modified()
    canvas.Update()

    return latex

#calc standard error on linear fit
#do error bars on overall fit, n = # pts, y_i = beta0 + beta1 * x + e
#standard error on residuals: SE = sqrt(sum_i (y_i - hat{y}_i)^2 / (n-2))
#standard error on slope: SE_beta1 = SE/sqrt(sum_i (x_i - bar{x})^2 )
#standard error on intercept: SE_beta0 = SE*sqrt((1/n) + bar{x}^2 / sum_i (x_i - bar{x})^2)
def calc_fit_err(graph, fit_func):
    npts = graph.GetN()
    x_vals = []
    y_vals = []
    res_sq = []
    for i in range(npts):
        x, y = 0, 0
        #graph.GetPoint(i,x,y)
        x = float(graph.GetPointX(i))
        y = float(graph.GetPointY(i))
        x_vals.append(x)
        y_vals.append(y)
        yhat = float(fit_func.Eval(x))
        res_sq.append((y - yhat)**2)
    bar_x = np.mean(x_vals)
    x_diff_sq = [(i - bar_x)**2 for i in x_vals]
    se = sum(res_sq) / (npts - 2)
    se_slope = se/math.sqrt( sum(x_diff_sq) )
    se_int = se*math.sqrt( (1/npts) + (bar_x**2)/(sum(x_diff_sq)) )

    nom_slope = fit_func.GetParameter(1)
    nom_int = fit_func.GetParameter(0)
    slope_up = nom_slope + se_slope
    slope_dn = nom_slope - se_slope
    int_up = nom_int + se_int
    int_dn = nom_int - se_int

    print("nom_slope",nom_slope,"nom_int",nom_int)
    print("slope_up",slope_up,"int_up",int_up)
    print("slope_dn",slope_dn,"int_dn",int_dn)

    #finely sample fit function to draw error band
    n_finepts = 100
    errband_yup = []
    errband_ydn = []
    errband_x = np.linspace(min(x_vals), max(x_vals), n_finepts)
    errband_y = np.linspace(min(y_vals), max(y_vals), n_finepts)
    errband_xerr = [0 for i in range(n_finepts)]
    for i in range(n_finepts):
        x = errband_x[i]
        y = errband_y[i]
        
        y_up = x*slope_up + int_up
        y_dn = x*slope_dn + int_dn
        errband_yup.append(y_up - y)
        errband_ydn.append(y - y_dn)

    errband_gr = ROOT.TGraphAsymmErrors(n_finepts, array.array("d",errband_x), array.array("d",errband_y), array.array("d",errband_xerr), array.array("d",errband_xerr), array.array("d",errband_ydn), array.array("d",errband_yup))
    errband_gr.SetFillColorAlpha(ROOT.kBlue, 0.3)  # Blue with 30% opacity
    errband_gr.SetFillStyle(3001) 
    return errband_gr, (int_up, slope_up), (int_dn, slope_dn) 
        
 
def flatten_cuts(obj):
    """Flatten arbitrarily nested YAML cut lists into a list of expressions."""
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, (list, tuple)):
        out = []
        for item in obj:
            out.extend(flatten_cuts(item))
        return out
    if obj is None:
        return []
    raise TypeError(f"Unexpected cut structure in YAML: {type(obj).__name__}")


def combine_and(cuts):
    cuts = [f"({c})" for c in cuts if c and str(c).strip()]
    return " && ".join(cuts) if cuts else "true"


def extract_grid(fileprocessor):
    """Discover mass/lifetime tuples that exist in both FullSim and FastSim."""
    def parse(files):
        points = set()
        for f in files:
            name = f.rsplit("/", 1)[-1]
            m = re.search(
                #r"mGl-(\d+)_mN2-(\d+)_mN1-(\d+)_ct([^_]+)", name
        r"mGl-(\d+)_mN2-(\d+)_mN1-(\d+)_GZ_N2ctau-([^_]+)", name
            )
            if m:
                points.add((int(m.group(1)), int(m.group(2)),
                            int(m.group(3)), m.group(4)))
        return points

    full = parse(fileprocessor._procFiles.get("SMS_gogoGZ", []))
    fast = parse(fileprocessor._procFiles.get("SMS_gogoGZFast", []))
    return sorted(full & fast, key=lambda x: (x[0], x[1], x[2], x[3]))


def get_files_checked(fp, process, masses, ctau):
    files = fp.GetFiles(process, masses[0], masses[1], masses[2], ctau)
    if not files:
        return []
    return list(files)


def sample_available(fp, key):
    return key in getattr(fp, "_procFiles", {})


def get_common_cut_strings(config, sr_bin_names):
    """Return nominal preselection and SR cuts, reusing YAML-defined SR bins."""
    baseline = flatten_cuts(config.get("baseline_cuts", []))
    kin = flatten_cuts(config.get("kin", []))

    # Keep the same trigger-free cleaning used in the original calculator:
    # the provided script explicitly notes that FastSim has no trigger branches.
    cleaning = [
        "((Flag_MetFilters == 1) && (rjrPTS[0] < 200) && (!Flag_hemVeto))"
    ]
    preselection = combine_and(baseline + cleaning + kin)

    bins = config.get("bins", {})
    missing = [name for name in sr_bin_names if name not in bins]
    if missing:
        raise KeyError(
            "Requested SR bins are absent from the YAML: " + ", ".join(missing)
        )

    # Use the region-specific YAML cuts, but replace the trigger-bearing
    # Cleaning alias with the trigger-free FastSim-compatible expression above.
    # Baseline, kin, noSV, Ch8 SR selection and delayed-bin cuts remain YAML-driven.
    per_bin = []
    for name in sr_bin_names:
        cuts = flatten_cuts(bins[name].get("cuts", []))
        clean_yaml = set(flatten_cuts(config.get("Cleaning", [])))
        cuts = [c for c in cuts if c not in clean_yaml]
        cuts += cleaning
        per_bin.append(combine_and(cuts))
    sr = " || ".join(f"({c})" for c in per_bin)
    return preselection, sr


def decay_weight_expression(br_gamma):
    """Build the exact nominal-0.5/0.5 to target-BR reweight expression."""
    br_gamma = float(br_gamma)
    z_new = 1.0 - br_gamma
    return (
        f"evtFillWgt * PhotonIsoDecayReweight(Evt_isZZ, Evt_isGZ, Evt_isGG, "
        f"0.5, {z_new:.17g}, 0.5, {br_gamma:.17g})"
    )


def book_sample_actions(files, preselection, sr, pt_edges, score_column, br_values, point_tag):
    """Book weighted sums and score histograms for every BR on one RDF graph."""
    if not files:
        return None
    df = RDataFrame("kuSkimTree", files)
    pre = df.Filter(preselection, "nominal preselection")
    pre_photon = pre.Filter("nBaseLinePhotons > 0", "at least one baseline photon")
    selected = pre_photon.Filter(sr, "Ch8 late tight-isolation SR")
    result = {}
    for ib, br in enumerate(br_values):
        weight_col = f"decayWeight_{point_tag}_{ib}"
        # Define on the original dataframe so both denominator and numerator
        # use precisely the same event-level reweight.
        weighted_df = df.Define(weight_col, decay_weight_expression(br)).Define(f"{weight_col}_sq",f"{weight_col}*{weight_col}")
        weighted_pre = weighted_df.Filter(preselection, "nominal preselection")
        weighted_pre_photon = weighted_pre.Filter("nBaseLinePhotons > 0", "at least one baseline photon")
        weighted_sr = weighted_pre_photon.Filter(sr, "Ch8 late tight-isolation SR")
        actions = {
            "pre_sum": weighted_pre.Sum(weight_col),
            "pre_sum_sq": weighted_pre.Sum(f"{weight_col}_sq"),
            "sr_sum": weighted_sr.Sum(weight_col),
            "sr_sum_sq": weighted_sr.Sum(f"{weight_col}_sq"),
            "pre_raw": weighted_pre.Count(),
            "sr_raw": weighted_sr.Count(),
            "pt": [],
        }
        if abs(br - 0.5) < 1e-10:
            actions["score_pre"] = weighted_pre_photon.Histo1D(
                (f"h_score_pre_{point_tag}", ";Baseline photon iso score;Weighted events", 50, 0.2, 1.05),
                score_column, weight_col
            )
            actions["score_sr"] = weighted_sr.Histo1D(
                (f"h_score_sr_{point_tag}", ";Baseline photon iso score;Weighted events", 50, 0.2, 1.05),
                score_column, weight_col
            )
        for ipt, (lo, hi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
            pt_filter = f"(baseLinePhoton_Pt[0] >= {lo}) && (baseLinePhoton_Pt[0] < {hi})"
            pre_pt = weighted_pre_photon.Filter(pt_filter, f"pre pT {lo}-{hi}")
            sr_pt = weighted_sr.Filter(pt_filter, f"SR pT {lo}-{hi}")
            actions["pt"].append((pre_pt.Sum(weight_col), sr_pt.Sum(weight_col),
                                  pre_pt.Sum(f"{weight_col}_sq"), sr_pt.Sum(f"{weight_col}_sq"),
                                  pre_pt.Count(), sr_pt.Count()))
        result[br] = actions
    return result


def materialize_br_actions(actions, br_values):
    """Evaluate booked RDF actions and return plain Python values/histograms."""
    if actions is None:
        return None
    out = {}
    for br in br_values:
        a = actions[br]
        item = {
            "pre": float(a["pre_sum"].GetValue()),
            "sr": float(a["sr_sum"].GetValue()),
            "pre_sq": float(a["pre_sum_sq"].GetValue()),
            "sr_sq": float(a["sr_sum_sq"].GetValue()),
            "pre_raw": int(a["pre_raw"].GetValue()),
            "sr_raw": int(a["sr_raw"].GetValue()),
            "pt": [],
        }
        for pre_sum, sr_sum, pre_sq, sr_sq, pre_count, sr_count in a["pt"]:
            item["pt"].append((float(pre_sum.GetValue()), float(sr_sum.GetValue()),
                                float(pre_sq.GetValue()), float(sr_sq.GetValue()),
                                int(pre_count.GetValue()), int(sr_count.GetValue())))
        for hname in ("score_pre", "score_sr"):
            if hname in a:
                hist = a[hname].GetValue().Clone()
                hist.SetDirectory(0)
                item[hname] = hist
        out[br] = item
    return out


def safe_eff(numerator, denominator):
    return float(numerator) / denominator if denominator > 0 else math.nan


def efficiency_sf(full, fast):
    ef = safe_eff(full["sr"], full["pre"])
    ea = safe_eff(fast["sr"], fast["pre"])
    sf = ef / ea if ea > 0 and math.isfinite(ef) and math.isfinite(ea) else math.nan
    return sf, ef, ea


def _format_series_value(key, value):
    """Readable label for a grouping value."""
    if key == "ctau":
        try:
            return f"{float(str(value).replace('p', '.')):g}"
        except (TypeError, ValueError):
            return str(value)
    try:
        numeric = float(value)
        if numeric.is_integer():
            return str(int(numeric))
        return f"{numeric:g}"
    except (TypeError, ValueError):
        return str(value)


def _row_matches_filters(row, filters):
    for key, expected in filters.items():
        actual = row.get(key)
        try:
            if abs(float(actual) - float(expected)) > 1e-8:
                return False
        except (TypeError, ValueError):
            if str(actual) != str(expected):
                return False
    return True


def make_graph_plot(rows, x_key, x_label, output_prefix, root_file,
                    series_keys=None, filters=None, plot_name=None):
    """Plot SF vs a chosen variable, with curves grouped by arbitrary columns."""
    series_keys = series_keys or []
    filters = filters or {}
    plot_name = plot_name or x_key
    grouped = defaultdict(list)
    for row in rows:
        if not _row_matches_filters(row, filters):
            continue
        value = row.get(x_key)
        sf = row.get("sf")
        eff_full = float(row.get("eff_full"))
        eff_fast = float(row.get("eff_fast"))
        pre_full = float(row.get("pre_full"))
        pre_fast = float(row.get("pre_fast"))
        pre_sq_full = float(row.get("pre_sq_full"))
        pre_sq_fast = float(row.get("pre_sq_fast"))
        
        if value is None or sf is None:
            continue
        try:
            x, y = float(value), float(sf)
            #binomial error on efficiency is err = sqrt( e(1-e) / n) where n = n_eff = (sum_i w_i)**2 / sum_i (w_i**2) for weighted events
            neff_full = (pre_full**2) / pre_sq_full
            neff_fast = (pre_fast**2) / pre_sq_fast
            err_full = math.sqrt( eff_full*(1 - eff_full)/neff_full )
            err_fast = math.sqrt( eff_fast*(1 - eff_fast)/neff_fast )
            err_sf = float(sf)*math.sqrt( (err_full/eff_full)**2 + (err_fast/eff_fast)**2 )
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        if series_keys:
            label = ", ".join(
                f"{key}={_format_series_value(key, row.get(key, 'NA'))}"
                for key in series_keys
            )
        else:
            label = "inclusive"
        grouped[label].append((x, y, err_sf))

    if not grouped:
        print(f"NOTE: no valid points for plot '{plot_name}' after filters {filters}")
        return

    safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", plot_name).strip("_")
    c = TCanvas(f"c_{safe_name}", f"SF vs {x_label}", 900, 700)
    c.SetGrid()
    graphs = []
    colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kOrange + 7, ROOT.kCyan + 2, ROOT.kBlack, ROOT.kViolet + 1,
              ROOT.kAzure + 2, ROOT.kPink + 7]
    overall_pts = []
    for i, (series, points) in enumerate(sorted(grouped.items())):
        # Multiple rows can share x when an unfiltered variable is not included
        # in the series definition. Keep them visible rather than silently averaging.
        points.sort()
        graph = TGraphErrors(len(points))
        graph.SetName(f"g_sf_{safe_name}_{i}")
        graph.SetTitle(f";{x_label};FullSim/FastSim efficiency SF")
        for j, (x, y, err) in enumerate(points):
            graph.SetPoint(j, x, y)
            graph.SetPointError(j, 0.0, err)
            overall_pts.append((x,y,err))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetMarkerColor(colors[i % len(colors)])
        graph.SetMarkerStyle(20 + i % 5)
        graph.SetLineWidth(2)
        graphs.append((series, graph))

    if len([i[1] for i in graphs]) > 0:
        #do fit to all points if multiple series
        graph = TGraphErrors(len(overall_pts))
        graph.SetName(f"g_sf_{safe_name}_total")
        graph.SetTitle(f";{x_label};FullSim/FastSim efficiency SF")
        for j, (x, y, err) in enumerate(overall_pts):
            graph.SetPoint(j, x, y)
            graph.SetPointError(j, 0.0, err)
        graph.SetLineColor(ROOT.kBlack)
        graph.SetMarkerColor(ROOT.kBlack)
        graph.SetMarkerStyle(20 + (len(graphs)+1) % 5)
        graph.SetLineWidth(3)
        graphs.append(("All points", graph))
   
    


    first = True
    allfit = None
    for series, graph in graphs:
        graph.Fit("pol1")
        fit_func = graph.GetFunction("pol1")
        fit_func.SetLineColor(graph.GetLineColor())
        fit_func.SetLineWidth(graph.GetLineWidth())
        if series == "All points" and graph.GetN() > 2:
            #do error on overall fit
            fiterr_gr, upfit, dnfit = calc_fit_err(graph, fit_func)
            #graph.Draw("A" if first else "A SAME")
            fit_func.Draw("SAME")
            fiterr_gr.Draw("3SAME")
            allfit = (fit_func, upfit, dnfit)
        else:
            graph.Draw("AP" if first else "P SAME")
        if first:
            graph.GetYaxis().SetTitle("FullSim / FastSim efficiency SF")
            graph.GetXaxis().SetTitle(x_label)
            graph.SetMinimum(0.0)
            first = False
    leg = TLegend(0.56, 0.68, 0.89, 0.89)
    for series, graph in graphs:
        if series == "All points":
            leg.AddEntry(graph, series, "l")
            continue
        leg.AddEntry(graph, series, "lp")
    leg.Draw()
    #draw text of allfit
    draw_fit_equation(c, allfit[0], allfit[1], allfit[2],
                     x=0.15, y=0.85, text_size=0.035)

    c.Modified()
    c.Update()
    c.Write()
    c.SaveAs(f"{output_prefix}_vs_{safe_name}.png")


def print_plot_spec_options(columns):
    """Print the available columns and syntax for custom plot specifications."""
    preferred_order = [
        "mGl", "mN2", "mN1", "deltaM_gl_N2", "deltaM_N2_N1",
        "ctau", "ctau_numeric", "BRgamma", "BRZ", "topology",
        "pre_full", "sr_full", "pre_fast", "sr_fast",
        "eff_full", "eff_fast", "sf", "pt_low", "pt_high", "pt_center",
    ]
    ordered = [key for key in preferred_order if key in columns]
    ordered.extend(sorted(columns - set(ordered)))
    print("\nAvailable --plot-spec options:")
    print("  x=<column>       Choose any output CSV column for the x-axis.")
    print("  series=<columns> Split into separate curves by one or more comma-separated columns.")
    print("  filter=<pairs>   Select rows using comma-separated column:value pairs.")
    print("  kind=inclusive|pt|both   Choose inclusive SF rows, pT-binned rows, or both.")
    print("  label=<text>     Optional output plot name (defaults to the x column).")
    print("  Available columns (valid for x, series, and filter where applicable):")
    print("    " + ", ".join(ordered))
    print("  Note: pt_low, pt_high, and pt_center are only populated for pT-binned rows;")
    print("        use kind=pt for those. ctau_numeric is the numeric form of ctau.")
    print("  Example: --plot-spec 'x=deltaM_N2_N1;series=mGl;filter=BRgamma:0.5,ctau_numeric:0.1;kind=inclusive'\n")


def parse_plot_spec(spec):
    """Parse e.g. 'x=deltaM_N2_N1;series=mGl,ctau;filter=BRgamma:0.5;kind=inclusive'."""
    parsed = {}
    for part in spec.split(";"):
        if not part.strip():
            continue
        if "=" not in part:
            raise ValueError(f"Invalid --plot-spec component '{part}' in '{spec}'")
        key, value = part.split("=", 1)
        parsed[key.strip().lower()] = value.strip()
    if "x" not in parsed:
        raise ValueError(f"--plot-spec needs x=...: {spec}")
    unknown = set(parsed) - {"x", "series", "filter", "kind", "label"}
    if unknown:
        raise ValueError(f"Unknown --plot-spec field(s): {', '.join(sorted(unknown))}")
    series_keys = [x.strip() for x in parsed.get("series", "").split(",") if x.strip()]
    filters = {}
    filter_text = parsed.get("filter", "")
    if filter_text:
        for item in filter_text.split(","):
            if ":" not in item:
                raise ValueError(f"Filter '{item}' must use column:value syntax")
            key, value = item.split(":", 1)
            filters[key.strip()] = value.strip()
    kind = parsed.get("kind", "inclusive").lower()
    if kind not in {"inclusive", "pt", "both"}:
        raise ValueError("plot kind must be inclusive, pt, or both")
    x_key = parsed["x"].strip()
    label = parsed.get("label", x_key)
    return {"x": x_key, "series": series_keys, "filters": filters,
            "kind": kind, "label": label}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--pt-edges", type=float, nargs="+", default=DEFAULT_PT_EDGES)
    parser.add_argument("--sr-bins", nargs="+", default=SR_BIN_NAMES,
                        help="YAML bin names whose union defines the signal region.")
    parser.add_argument("--score", default="baseLinePhoton_isoANNScore")
    parser.add_argument("--br-values", type=float, nargs="+", default=DEFAULT_BR_VALUES,
                        help="Target BR(N2->gamma N1); the Z BR is 1-BRgamma.")
    parser.add_argument("--max-points", type=int, default=0,
                        help="Optional debugging limit; 0 processes all matched points.")
    parser.add_argument(
        "--plot-spec", action="append", default=[],
        help=("Custom plot recipe; repeatable. Syntax: "
              "'x=deltaM_N2_N1;series=mGl;filter=BRgamma:0.5;kind=inclusive'. "
              "Columns can be any output CSV field. series accepts comma-separated columns; "
              "filter accepts comma-separated column:value pairs; kind is inclusive, pt, or both. "
              "If omitted, the framework's standard plots are made.")
    )
    args = parser.parse_args()

    if len(args.pt_edges) < 2 or any(b <= a for a, b in zip(args.pt_edges[:-1], args.pt_edges[1:])):
        parser.error("--pt-edges must be strictly increasing and contain at least two values")
    if not args.br_values or any(b < 0 or b > 1 for b in args.br_values):
        parser.error("--br-values must be non-empty and lie in [0, 1]")
    br_values = list(dict.fromkeys(float(x) for x in args.br_values))
    if 0.5 not in br_values:
        br_values.append(0.5)  # nominal weights also provide score-shape diagnostics

    with open(args.config, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    preselection, sr = get_common_cut_strings(config, args.sr_bins)
    fp = FileProcessor()
    points = extract_grid(fp)
    if args.max_points:
        points = points[:args.max_points]
    if not points:
        raise RuntimeError("No matching FullSim/FastSim gogoGZ mass/lifetime points found.")

    print(f"Matched {len(points)} FullSim/FastSim gogoGZ parameter points.")
    print(f"Preselection: {preselection}")
    print(f"Signal region is the union of: {', '.join(args.sr_bins)}")
    print("BR reweight uses Evt_isZZ/Evt_isGZ/Evt_isGG and the supplied C++ function.")
    print("The leading baseline photon baseLinePhoton_Pt[0] defines photon-pT bins.")

    booked = {}
    for ipoint, (mgl, mn2, mn1, ctau) in enumerate(points):
        point = (mgl, mn2, mn1, ctau)
        masses = (mgl, mn2, mn1)
        for sim, process in (("full", "SMS_gogoGZ"), ("fast", "SMS_gogoGZFast")):
            files = get_files_checked(fp, process, masses, ctau)
            if not files:
                continue
            tag = f"{sim}_{ipoint}_{mgl}_{mn2}_{mn1}_{ctau}".replace(".", "p")
            booked[(point, sim)] = book_sample_actions(
                files, preselection, sr, args.pt_edges, args.score, br_values, tag
            )

    # Trigger all the already-booked actions. Multiple BR weights share each
    # underlying dataframe graph, so ROOT can reuse the event loop.
    values = {key: materialize_br_actions(action, br_values)
              for key, action in booked.items()}

    rows = []
    nominal_score_hists = []
    for mgl, mn2, mn1, ctau in points:
        point = (mgl, mn2, mn1, ctau)
        full_by_br = values.get((point, "full"))
        fast_by_br = values.get((point, "fast"))
        if full_by_br is None or fast_by_br is None:
            print(f"WARNING: missing FullSim/FastSim pair for {point}; skipping.", file=sys.stderr)
            continue
        for br in br_values:
            full, fast = full_by_br[br], fast_by_br[br]
            sf, full_eff, fast_eff = efficiency_sf(full, fast)
            row_base = {
                "mGl": mgl, "mN2": mn2, "mN1": mn1,
                "deltaM_gl_N2": mgl - mn2, "deltaM_N2_N1": mn2 - mn1,
                "ctau": ctau, "BRgamma": br, "BRZ": 1.0 - br,
                "topology": "GZ sample reweighted", "series": f"BR(gamma)={br:g}",
                "pre_full": full["pre"], "sr_full": full["sr"],
                "pre_fast": fast["pre"], "sr_fast": fast["sr"],
                "pre_sq_full": full["pre_sq"], "sr_sq_full": full["sr_sq"],
                "pre_sq_fast": fast["pre_sq"], "sr_sq_fast": fast["sr_sq"],
                "pre_raw_full": full["pre_raw"], "sr_raw_full": full["sr_raw"],
                "pre_raw_fast": fast["pre_raw"], "sr_raw_fast": fast["sr_raw"],
                "eff_full": full_eff, "eff_fast": fast_eff, "sf": sf,
            }
            rows.append(row_base)
            if abs(br - 0.5) < 1e-10 and "score_pre" in full and "score_pre" in fast:
                nominal_score_hists.append((row_base, full["score_pre"], fast["score_pre"],
                                            full["score_sr"], fast["score_sr"]))
            for i, (lo, hi) in enumerate(zip(args.pt_edges[:-1], args.pt_edges[1:])):
                npre_f, nsr_f, npre_sq_f, nsr_sq_f, npre_f_raw, nsr_f_raw = full["pt"][i]
                npre_a, nsr_a, npre_sq_a, nsr_sq_a, npre_a_raw, nsr_a_raw = fast["pt"][i]
                ef = safe_eff(nsr_f, npre_f)
                ea = safe_eff(nsr_a, npre_a)
                pt_sf = ef / ea if ea > 0 and math.isfinite(ef) else math.nan
                rows.append({
                    **row_base, "series": f"BR={br:g}, mGl={mgl},mN2={mn2},mN1={mn1},ct={ctau}",
                    "pt_low": lo, "pt_high": hi, "pt_center": 0.5 * (lo + hi),
                    "pre_full": npre_f, "sr_full": nsr_f,
                    "pre_fast": npre_a, "sr_fast": nsr_a,
                    "pre_sq_full": npre_sq_f, "sr_sq_full": nsr_sq_f,
                    "pre_sq_fast": npre_sq_a, "sr_sq_fast": nsr_sq_a,
                    "pre_raw_full": npre_f_raw, "sr_raw_full": nsr_f_raw,
                    "pre_raw_fast": npre_a_raw, "sr_raw_fast": nsr_a_raw,
                    "eff_full": ef, "eff_fast": ea, "sf": pt_sf,
                })

    if not rows:
        raise RuntimeError("No SF rows were produced; check files, branches, and YAML selections.")
    fields = sorted({key for row in rows for key in row})
    csv_path = f"{args.output}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} yield/SF rows to {csv_path}")

    root_path = f"{args.output}.root"
    root_file = TFile.Open(root_path, "RECREATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not create ROOT output {root_path}")
    root_file.cd()
    iso_dir = root_file.mkdir("IsoScoreHists")
    iso_dir.cd()
    for row, hfp, hap, hfs, has in nominal_score_hists:
        label = f"mGl{row['mGl']}_mN2{row['mN2']}_mN1{row['mN1']}_ct{row['ctau']}"
        hfp.Write(f"isoScore_full_pre_{label}")
        hap.Write(f"isoScore_fast_pre_{label}")
        hfs.Write(f"isoScore_full_sr_{label}")
        has.Write(f"isoScore_fast_sr_{label}")
    root_file.cd()
    inclusive = [r for r in rows if "pt_center" not in r]
    pt_rows = [r for r in rows if "pt_center" in r]
    for row in inclusive:
        try:
            row["ctau_numeric"] = float(str(row["ctau"]).replace("p", "."))
        except (ValueError, TypeError):
            row["ctau_numeric"] = math.nan
    for row in pt_rows:
        try:
            row["ctau_numeric"] = float(str(row["ctau"]).replace("p", "."))
        except (ValueError, TypeError):
            row["ctau_numeric"] = math.nan

    x_labels = {
        "mGl": "m_{#tilde{g}} [GeV]",
        "mN2": "m_{#tilde{#chi}^{0}_{2}} [GeV]",
        "mN1": "m_{#tilde{#chi}^{0}_{1}} [GeV]",
        "deltaM_gl_N2": "m_{#tilde{g}} - m_{#tilde{#chi}^{0}_{2}} [GeV]",
        "deltaM_N2_N1": "m_{#tilde{#chi}^{0}_{2}} - m_{#tilde{#chi}^{0}_{1}} [GeV]",
        "BRgamma": "BR(N_{2} #rightarrow #gamma N_{1})",
        "BRZ": "BR(N_{2} #rightarrow Z N_{1})",
        "ctau_numeric": "c#tau [cm]",
        "pt_center": "Leading baseline photon p_{T} [GeV]",
        "pt_low": "Photon p_{T} bin low edge [GeV]",
        "pt_high": "Photon p_{T} bin high edge [GeV]",
    }
    known_columns = {key for row in rows for key in row}
    known_columns.add("ctau_numeric")
    # Always print the current options so users can discover plot-spec values
    # from the actual output schema without having to inspect the CSV manually.
    print_plot_spec_options(known_columns)
    if args.plot_spec:
        try:
            specs = [parse_plot_spec(s) for s in args.plot_spec]
        except ValueError as exc:
            parser.error(str(exc))
        for i, spec in enumerate(specs):
            if spec["x"] not in known_columns:
                parser.error(f"Unknown x column '{spec['x']}'. Check the output CSV column names.")
            for key in spec["series"]:
                if key not in known_columns:
                    parser.error(f"Unknown series column '{key}'. Check the output CSV column names.")
            data_sets = []
            if spec["kind"] in {"inclusive", "both"}:
                data_sets.append((inclusive, "inclusive"))
            if spec["kind"] in {"pt", "both"}:
                data_sets.append((pt_rows, "pt"))
            for data, suffix in data_sets:
                # Build a descriptive name from the plot specification
                name_parts = [spec["label"]]
                # Include the curve grouping variables
                if spec["series"]:
                    name_parts.append("series_" + "_".join(spec["series"]))
                else:
                    name_parts.append("series_inclusive")
                # Include the applied filters and their values
                if spec["filters"]:
                    filter_tag = "_".join(
                        f"{key}_{value.replace('.','p')}" for key, value in sorted(spec["filters"].items())
                    )
                    name_parts.append("filter_" + filter_tag)
                else:
                    name_parts.append("filter_none")
                # Include the plot kind (inclusive or pT-binned)
                name_parts.append("kind_" + suffix)
                plot_name = "_".join(name_parts)
                make_graph_plot(
                    data,
                    spec["x"],
                    x_labels.get(spec["x"], spec["x"]),
                    args.output,
                    root_file,
                    spec["series"],
                    spec["filters"],
                    plot_name=plot_name,
                )
    else:
        # Standard plots retain the earlier overview, while --plot-spec provides
        # arbitrary x-axis, curve grouping, and selection control.
        for key in ("mGl", "mN2", "mN1", "deltaM_gl_N2", "deltaM_N2_N1"):
            make_graph_plot(inclusive, key, x_labels[key], args.output, root_file,
                            ["BRgamma"])
        make_graph_plot(inclusive, "BRgamma", x_labels["BRgamma"], args.output,
                        root_file, ["mGl"])
        make_graph_plot(inclusive, "ctau_numeric", x_labels["ctau_numeric"], args.output,
                        root_file, ["BRgamma"])
        make_graph_plot(pt_rows, "pt_center", x_labels["pt_center"], args.output + "_pt",
                        root_file, ["BRgamma", "mGl"])
    root_file.Write()
    root_file.Close()
    print(f"Wrote ROOT output and plots with prefix {args.output}")


if __name__ == "__main__":
    main()
