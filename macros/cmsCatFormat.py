import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import ROOT
import uproot
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal.*") #suppress warning about 0 being the smallest subnormal number (not relevant for plotting)

def get_plots_from_file(file_name, hists):
    infile = ROOT.TFile.Open(file_name)
    histos = [infile.Get(hist) for hist in hists]
    return histos


def make_plot(inhist, inhist_name):
    fig, ax = plt.subplots()
    inhist.plot2d(ax = ax, cbarextend = True, norm=LogNorm(), rasterized=True) #rasterized option reduces visual gaps between bins
    fig.get_axes()[0].set_ylabel('eta', fontsize=20) #yaxis
    fig.get_axes()[0].set_xlabel('time [ns]', fontsize=20) #xaxis
    fig.get_axes()[-1].set_ylabel('a.u.', fontsize=20, labelpad = 0.) #zaxis
    
    hep.cms.label(llabel="Preliminary",rlabel="(13 TeV)")
    
    plottitle = inhist_name+".pdf"
    print("Saving histogram",plottitle)
    fig.savefig(plottitle)

def main():
    #hep.style.use("CMS")
    plt.style.use(hep.style.CMS)
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    
    hists_to_plot = ["SC_TimeCenter_EtaCenter","SC_TimeCenter_EtaCenter_BHFilter"]

    file_name = "plots/condor_superclusters_defaultv8_beta0-1e-5_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NperGeV-0p0333333_emAlpha-1e-5_MET_R22_AL1NpSC_v31_MET_detbkg.root"

    hists = get_plots_from_file(file_name, hists_to_plot) 

    for h, hist in enumerate(hists):
        infile = uproot.open(file_name)
        hist_name = hists_to_plot[h]
        hist = infile[hist_name]
        hist = hist.to_hist()
        make_plot(hist, hist_name)


if __name__ == "__main__":
    main()
