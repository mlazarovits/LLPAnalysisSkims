import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import ROOT
import uproot
import argparse
import os
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal.*") #suppress warning about 0 being the smallest subnormal number (not relevant for plotting)

def get_plots_from_file_TFile(file_name, hists):
    infile = ROOT.TFile.Open(file_name)
    histos = [infile.Get(hist) for hist in hists]
    return histos


def get_plots_from_file(file_name, hists):
    infile = uproot.open(file_name)
    histos = [infile[hist].to_hist() for hist in hists]
    return histos

def make_plot_2d(inhist, inhist_name, pathname):
    fig, ax = plt.subplots(figsize=(12,10))
    inhist.plot2d(ax = ax, cbarextend = True, norm=LogNorm(), rasterized=True) #rasterized option reduces visual gaps between bins
    fig.get_axes()[0].set_ylabel('eta', fontsize=20) #yaxis
    fig.get_axes()[0].set_xlabel('time [ns]', fontsize=20) #xaxis
    fig.get_axes()[-1].set_ylabel('a.u.', fontsize=20, labelpad = 0.) #zaxis
    
    hep.cms.label(llabel="Preliminary",rlabel="(13 TeV)")
    
    plottitle = pathname+"/"+inhist_name+".pdf"
    print("Saving histogram",plottitle)
    fig.savefig(plottitle)

def main():
    hep.style.use("CMS")
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputFile",'-i',help='input file with histograms to be formatted',required=True)
    parser.add_argument("--hists1D",help='1D hists to format',nargs="+",default=[])
    parser.add_argument("--hists2D",help='2D hists to format',nargs="+",default=[])
    args = parser.parse_args()

    if(len(args.hists1D) < 1 and len(args.hists2D) < 0):
        print("No histograms specified.")
        return

    pathname = os.path.basename(args.inputFile)
    pathname = pathname[:pathname.rfind('.')]
    pathname = "formatted_plots/"+pathname
    if not os.path.exists(pathname):
        os.mkdirs(pathname)


    #add formatting for 1d hists

    hists2d = get_plots_from_file(args.inputFile, args.hists2D) 
    for h, hist in enumerate(hists2d):
        hist_name = args.hists2D[h]
        make_plot_2d(hist, hist_name, pathname)


if __name__ == "__main__":
    main()
