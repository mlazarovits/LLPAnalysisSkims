import os
import re
from ROOT import gInterpreter, kAzure, kRed, TLegend, TLatex
import cmsstyle as CMS

FULLSIM_COLOR = kAzure + 1
FASTSIM_COLOR = kRed - 6

FULLSIM_MARKER = 20
FASTSIM_MARKER = 21


dist = {
	'nSelPhotons': {'xlabel': 'Number of selected photons',	'nbins': 5,  'xmin': 0,	'xmax': 5},
	'nIsoPho': {'xlabel': 'Number of signal photons', 'nbins': 5,  'xmin': 0,	'xmax': 5},
	'selPhoEta': {'xlabel': 'Selected Photon Pseudorapidity #eta',	 'nbins': 50,  'xmin': -3,	'xmax': 3},
	'isoPhoEta': {'xlabel': 'Signal Photon Pseudorapidity #eta',  'nbins': 50,  'xmin': -3,	'xmax': 3},
	'rawSeedTime': {'xlabel': 'selected seed time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'rawSeedTimeIso': {'xlabel': 'signal seed time [ns]',	  'nbins': 50,  'xmin': -3,   'xmax': 3},
	'rawSeedTimeLead': {'xlabel': 'lead selected seed time [ns]',	   'nbins': 50,  'xmin': -3,   'xmax': 3},
	'rawSeedTimeLeadIso': {'xlabel': 'lead signal seed time [ns]',	  'nbins': 50,  'xmin': -3,   'xmax': 3},
	'rawSeedTimeSublead': {'xlabel': 'sublead selected seed time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'rawSeedTimeSubleadIso': {'xlabel': 'sublead signal seed time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTime': {'xlabel': 'selected weighted time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTimeIso': {'xlabel': 'signal weighted time [ns]',	  'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTimeLead': {'xlabel': 'lead selected weighted time [ns]',	   'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTimeLeadIso': {'xlabel': 'lead signal weighted time [ns]',	  'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTimeSublead': {'xlabel': 'sublead selected weighted time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'weightedTimeSubleadIso': {'xlabel': 'sublead signal weighted time [ns]',	 'nbins': 50,  'xmin': -3,   'xmax': 3},
	'timeSig': {'xlabel': 'selected S_t',	 'nbins': 50,  'xmin': -5,	'xmax': 5},
	'timeSigIso': {'xlabel': 'signal S_t',	 'nbins': 50,  'xmin': -5,	'xmax': 5},
	'timeSigLead': { 'xlabel': 'lead selected S_t',	   'nbins': 50,  'xmin': -5,	'xmax': 5},
	'timeSigLeadIso': { 'xlabel': 'lead signal S_t',	   'nbins': 50,  'xmin': -5,	'xmax': 5},
	'timeSigSublead': {'xlabel': 'sublead selected S_t', 	   'nbins': 50,  'xmin': -5,	'xmax': 5},
	'timeSigSubleadIso': {'xlabel': 'sublead signal S_t',	   'nbins': 50,  'xmin': -5,	'xmax': 5},
	
}
def draw_cms_labels():
    """Draw standard CMS labels on plots"""
    cms_label = TLatex()
    cms_label.SetNDC()
    cms_label.SetTextSize(0.065)
    cms_label.SetTextFont(61)
    cms_label.DrawLatex(0.122, 0.945, "CMS")
    cms_label.SetTextFont(52)
    cms_label.SetTextSize(0.05)
    cms_label.DrawLatex(0.23, 0.945, "Simulation Preliminary")
    return cms_label


def parse_signal_name(stem):
    # Extract parameters using regex
    mgl_match = re.search(r'mGl_(\d+)', stem)
    mn2_match = re.search(r'mN2_(\d+)', stem)
    mn1_match = re.search(r'mN1_(\d+)', stem)
    ct_match = re.search(r'ct-?(\d+(?:p\d+)?)', stem)

    if all([mgl_match, mn2_match, mn1_match, ct_match]):
        mgl = mgl_match.group(1)
        mn2 = mn2_match.group(1)
        mn1 = mn1_match.group(1)
        ct = ct_match.group(1).replace('p', '.')

        return f"m_{{#tilde{{g}}}}({mgl})-m_{{#tilde{{#chi}}_{{2}}^{{0}}}}({mn2})-m_{{#tilde{{#chi}}_{{1}}^{{0}}}}({mn1}), c#tau={ct} m"
    else:
        return stem

def create_comparison_canvas(dist_key, h_full, h_fast, signal_label, logy=False, labels = ["FullSim", "FastSim"]):
	"""Create a single comparison canvas for one distribution.

	Both histograms are normalized to unit area.
	Returns (canvas, list_of_objects_to_keep_alive).
	"""
	key = dist_key
	xlabel = dist[dist_key]['xlabel']
	nbins =  dist[dist_key]['nbins']
	xmin =   dist[dist_key]['xmin']
	xmax =   dist[dist_key]['xmax']
	ylabel = 'Fraction of entries'


	canvas_name = f'c_{key}{signal_label[:-1]}'
	canvas = CMS.cmsCanvas(canvas_name, xmin, xmax, 0, 1., xlabel, ylabel,
						   square=False, extraSpace=0.01, iPos=0)
	canvas.SetCanvasSize(800, 600)
	canvas.SetLeftMargin(0.12)
	canvas.SetRightMargin(0.05)
	canvas.SetGrid()
	if logy:
		canvas.SetLogy()

	# Normalize to unit area
	if h_full.Integral() > 0:
		h_full.Scale(1.0 / h_full.Integral())
	if h_fast.Integral() > 0:
		h_fast.Scale(1.0 / h_fast.Integral())

	# Styling
	h_full.SetTitle("")
	h_full.GetXaxis().SetTitle(xlabel)
	h_full.GetYaxis().SetTitle(ylabel)
	h_full.SetLineColor(FULLSIM_COLOR)
	h_full.SetLineWidth(2)
	h_full.SetMarkerColor(FULLSIM_COLOR)
	h_full.SetMarkerStyle(FULLSIM_MARKER)
	h_full.SetMarkerSize(0.8)
	h_full.SetFillColorAlpha(FULLSIM_COLOR, 0.15)

	h_fast.SetTitle("")
	h_fast.GetXaxis().SetTitle(xlabel)
	h_fast.GetYaxis().SetTitle(ylabel)
	h_fast.SetLineColor(FASTSIM_COLOR)
	h_fast.SetLineWidth(2)
	h_fast.SetLineStyle(2)
	h_fast.SetMarkerColor(FASTSIM_COLOR)
	h_fast.SetMarkerStyle(FASTSIM_MARKER)
	h_fast.SetMarkerSize(0.8)

	# Y-axis range
	y_max = 1.4 * max(h_full.GetMaximum(), h_fast.GetMaximum())
	if y_max == 0:
		y_max = 1.0
	if logy:
		y_max = 1.2
		if dist_key == "nSelPhotons" or dist_key == "nIsoPho":
			y_max = 1.5

	h_full.SetStats(0)
	h_full.SetMaximum(y_max)
	h_full.SetMinimum(1e-4)
	h_full.GetXaxis().SetTitleSize(0.05)
	h_full.GetXaxis().SetLabelSize(0.05)
	h_full.GetXaxis().SetTitleOffset(1.25)
	h_full.GetXaxis().CenterTitle(True)
	h_full.GetYaxis().SetTitleSize(0.05)
	h_full.GetYaxis().SetLabelSize(0.05)
	h_full.GetYaxis().SetTitleOffset(1.25)
	h_full.GetYaxis().CenterTitle(True)

	hh_full = h_full.DrawCopy("HIST")
	hh_fast = h_fast.DrawCopy("HIST SAME")

	# Legend
	legend = TLegend(0.60, 0.75, 0.93, 0.88)
	legend.SetBorderSize(0)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)
	legend.AddEntry(hh_full, labels[0], "lf")
	legend.AddEntry(hh_fast, labels[1], "l")
	legend.Draw()

	# CMS labels
	cms_label = draw_cms_labels()

	# Signal sample label
	signal_label = parse_signal_name(signal_label)
	signal_latex = TLatex()
	signal_latex.SetNDC()
	signal_latex.SetTextSize(0.032)
	signal_latex.SetTextFont(42)
	signal_latex.DrawLatex(0.14, 0.88, signal_label)

	# Entries label
	entries_latex = TLatex()
	entries_latex.SetNDC()
	entries_latex.SetTextSize(0.028)
	entries_latex.SetTextFont(42)
	n_full = int(h_full.GetEntries())
	n_fast = int(h_fast.GetEntries())
	entries_latex.DrawLatex(0.14, 0.84, f"{labels[0]} entries: {n_full}, {labels[1]} entries: {n_fast}")

	objects = [h_full, h_fast, legend, cms_label, signal_latex, entries_latex]
	return canvas, objects

