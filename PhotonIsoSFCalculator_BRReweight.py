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

import argparse
import csv
import math
import re
import sys
import os
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
    "Ch8SRgeq1PhoNotBHLateTightIso00",
    "Ch8SRgeq1PhoNotBHLateTightIso10",
    "Ch8SRgeq1PhoNotBHLateTightIso01",
    "Ch8SRgeq1PhoNotBHLateTightIso11",
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


def extract_grid(fileprocessor, mgl = None, mn2 = None, mn1 = None, ct = None):
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
                if mgl is not None and str(m.group(1)) not in mgl:
                    continue
                if mn2 is not None and str(m.group(2)) not in mn2:
                    continue
                if mn1 is not None and str(m.group(3)) not in mn1:
                    continue
                if ct is not None and m.group(4) not in ct:
                    continue
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


def make_graph_plot(rows, x_key, x_label, output_prefix, root_file):
    """Save one TGraphErrors plot for a chosen parameter using ROOT."""
    grouped = defaultdict(list)
    for row in rows:
        value = row.get(x_key)
        sf = row.get("sf")
        eff_full = float(row.get("eff_full"))
        eff_fast = float(row.get("eff_fast"))
        pre_full = float(row.get("pre_full"))
        pre_fast = float(row.get("pre_fast"))
        pre_sq_full = float(row.get("pre_sq_full"))
        pre_sq_fast = float(row.get("pre_sq_fast"))
        if value is None or sf is None or not math.isfinite(float(sf)) or sf == 0.:
            continue
        try:
            x = float(value)
            y = float(sf)
            #binomial error on efficiency is err = sqrt( e(1-e) / n) where n = n_eff = (sum_i w_i)**2 / sum_i (w_i**2) for weighted events
            neff_full = (pre_full**2) / pre_sq_full
            neff_fast = (pre_fast**2) / pre_sq_fast
            err_full = math.sqrt( eff_full*(1 - eff_full)/neff_full )
            err_fast = math.sqrt( eff_fast*(1 - eff_fast)/neff_fast )
            err_sf = float(sf)*math.sqrt( (err_full/eff_full)**2 + (err_fast/eff_fast)**2 )
            print("eff_full",eff_full,"eff_fast",eff_fast,"sf",sf,"err_sf",err_sf)
        except (TypeError, ValueError):
            continue
        grouped[row.get("series", "inclusive")].append((x, y, err_sf))

    if not grouped:
        return

    c = TCanvas(f"c_{x_key}", f"SF vs {x_label}", 900, 700)
    c.SetGrid()
    graphs = []
    colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kOrange + 7, ROOT.kCyan + 2, ROOT.kBlack]
    for i, (series, points) in enumerate(sorted(grouped.items())):
        points.sort()
        graph = TGraphErrors(len(points))
        graph.SetName(f"g_sf_vs_{x_key}_{i}")
        graph.SetTitle(f";{x_label};FullSim/FastSim efficiency SF")
        for j, (x, y, err) in enumerate(points):
            graph.SetPoint(j, x, y)
            graph.SetPointError(j, 0.0, err)
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetMarkerColor(colors[i % len(colors)])
        graph.SetMarkerStyle(20 + i % 5)
        graph.SetLineWidth(2)
        graphs.append((series, graph))

    first = True
    for _, graph in graphs:
        #fit graph with linear function
        graph.Fit("pol1")
        graph.Draw("AP" if first else "P SAME")
        if first:
            graph.GetYaxis().SetTitle("FullSim / FastSim efficiency SF")
            graph.GetXaxis().SetTitle(x_label)
            graph.SetMinimum(0.0)
            first = False
    leg = TLegend(0.60, 0.70, 0.89, 0.89)
    for series, graph in graphs:
        leg.AddEntry(graph, series, "lp")
    leg.Draw()
    c.Modified()
    c.Update()
    c.Write()
    c.SaveAs(f"{output_prefix}_vs_{x_key}.png")


def main():
    #import kerberos credentials to conda env if not already there
    kerb = os.getenv("KRB5CCNAME")
    if(kerb is None):
        print("Setting kerebos credentials")
        os.environ["KRB5CCNAME"] = "API:"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG,help="yaml analysis config")
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--pt-edges", type=float, nargs="+", default=DEFAULT_PT_EDGES)
    parser.add_argument("--sr-bins", nargs="+", default=SR_BIN_NAMES,
                        help="YAML bin names whose union defines the signal region.")
    parser.add_argument("--score", default="baseLinePhoton_isoANNScore")
    parser.add_argument("--br-values", type=float, nargs="+", default=DEFAULT_BR_VALUES,
                        help="Target BR(N2->gamma N1); the Z BR is 1-BR_gamma.")
    parser.add_argument("--max-points", type=int, default=0,
                        help="Optional debugging limit; 0 processes all matched points.")
    parser.add_argument("--mGl",nargs="+",default=None,help="list of gluino masses to run over")
    parser.add_argument("--mN1",nargs="+",default=None,help="list of n2 masses to run over")
    parser.add_argument("--mN2",nargs="+",default=None,help="list of n1 masses to run over")
    parser.add_argument("--ctau",nargs="+",default=None,help="list of ctaus to run over")
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
    points = extract_grid(fp, args.mGl, args.mN2, args.mN1, args.ctau)
    print("points",points)
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
                "ctau": ctau, "BR_gamma": br, "BR_Z": 1.0 - br,
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
    for row, hfp, hap, hfs, has in nominal_score_hists:
        label = f"mGl{row['mGl']}_mN2{row['mN2']}_mN1{row['mN1']}_ct{row['ctau']}"
        hfp.Write(f"isoScore_full_pre_{label}")
        hap.Write(f"isoScore_fast_pre_{label}")
        hfs.Write(f"isoScore_full_sr_{label}")
        has.Write(f"isoScore_fast_sr_{label}")

    inclusive = [r for r in rows if "pt_center" not in r]
    pt_rows = [r for r in rows if "pt_center" in r]
    for key, label in (
        (f"mGl-{mgl}", "m_{#tilde{g}} [GeV]"),
        ("mN2", "m_{#tilde{#chi}^{0}_{2}} [GeV]"),
        ("mN1", "m_{#tilde{#chi}^{0}_{1}} [GeV]"),
        ("deltaM_gl_N2", "m_{#tilde{g}} - m_{#tilde{#chi}^{0}_{2}} [GeV]"),
        ("deltaM_N2_N1", "m_{#tilde{#chi}^{0}_{2}} - m_{#tilde{#chi}^{0}_{1}} [GeV]"),
        ("BR_gamma", "BR(N_{2} #rightarrow #gamma N_{1})"),
    ):
        make_graph_plot(inclusive, key, label, args.output, root_file)
    ctau_rows = []
    for row in inclusive:
        try:
            val = float(str(row["ctau"]).replace("p", "."))
        except ValueError:
            continue
        ctau_rows.append({**row, "ctau_numeric": val})
    make_graph_plot(ctau_rows, "ctau_numeric", "c#tau [cm]", args.output, root_file)
    make_graph_plot(pt_rows, "pt_center", "Leading baseline photon p_{T} [GeV]", args.output + "_pt", root_file)
    pt_mass_rows = [{**r, "series": f"BR={r['BR_gamma']:g}, pT {r['pt_low']:g}-{r['pt_high']:g} GeV"}
                    for r in pt_rows]
    make_graph_plot(pt_mass_rows, "mGl", "m_{#tilde{g}} [GeV]", args.output + "_ptbins", root_file)
    root_file.Write()
    root_file.Close()
    print(f"Wrote ROOT output and plots with prefix {args.output}.root")


if __name__ == "__main__":
    main()
