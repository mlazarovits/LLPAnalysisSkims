import numpy as np
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


def get_2d_plots_from_file(file_name, hists):
    infile = uproot.open(file_name)
    histos = []
    for hist in hists:
        h = infile[hist]
        counts, edges_x, edges_y = h.to_numpy()
        tot = np.sum(counts)
        vals = h.values()
        norm_vals = vals * (1/float(tot))
        h.values()[:] = norm_vals
        histos.append(h.to_hist())
    #histos = [infile[hist].to_hist() for hist in hists]
    return histos

def get_1d_plots_from_file(file_name, hists):
    infile = uproot.open(file_name)
    histos = []
    for hist in hists:
        h = infile[hist]
        counts, edges_x = h.to_numpy()
        tot = np.sum(counts)
        vals = h.values()
        norm_vals = vals * (1/float(tot))
        #h.values()[:] = norm_vals
        histos.append(h.to_hist())
    #histos = [infile[hist].to_hist() for hist in hists]
    return histos

def make_plot_2d(inhist, inhist_name, pathname, xlab, ylab, log):
    fig, ax = plt.subplots(figsize=(12,10))
    normlog = None
    if log:
        normlog = LogNorm()
    inhist.plot2d(ax = ax, cbarextend = True, norm=normlog, rasterized=True) #rasterized option reduces visual gaps between bins
    #inhist.plot2d(ax = ax, cbarextend = True, rasterized=True) #rasterized option reduces visual gaps between bins
    fig.get_axes()[0].set_ylabel(ylab, fontsize=20) #yaxis
    fig.get_axes()[0].set_xlabel(xlab, fontsize=20) #xaxis
    fig.get_axes()[-1].set_ylabel('a.u.', fontsize=20, labelpad = 0.) #zaxis
    
    hep.cms.label(llabel="Preliminary",rlabel="(13 TeV)")
    
    plottitle = pathname+"/"+inhist_name+".pdf"
    print("Saving histogram",plottitle)
    fig.savefig(plottitle)

def make_plot_1d(inhist, inhist_name, pathname, xlab, ylab, log, wide):
    fig = None
    ax = None
    if wide:
        fig, ax = plt.subplots(figsize=(20,6.5))
    else:
        fig, ax = plt.subplots(figsize=(12,10))
    if log:
        ax.set_yscale('log')
    inhist.plot(ax = ax)
    fig.get_axes()[0].set_ylabel(ylab, fontsize=20) #yaxis
    fig.get_axes()[0].set_xlabel(xlab, fontsize=20) #xaxis
    
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
    parser.add_argument("--log",help='log zaxis (2D) or yaxis (1D)',default=False,action='store_true')
    parser.add_argument("--wide",help='create wide figure size',default=False,action='store_true')
    args = parser.parse_args()

    if(len(args.hists1D) < 1 and len(args.hists2D) < 0):
        print("No histograms specified.")
        return

    pathname = os.path.basename(args.inputFile)
    pathname = pathname[:pathname.rfind('.')]
    pathname = "formatted_plots/"+pathname
    if not os.path.exists(pathname):
        os.makedirs(pathname)


    #add formatting for 1d hists
    hists1d = get_1d_plots_from_file(args.inputFile, args.hists1D)
    hists2d = get_2d_plots_from_file(args.inputFile, args.hists2D)
    

    #ylab = r"$\frac{E_{SC}}{p_{track}}$"
    for h, hist in enumerate(hists2d):
        hist_name = args.hists2D[h]
        xlab = ""
        ylab = ""
        if "map" in hist_name or "Map" in hist_name:
            xlab = "local ieta"
            ylab = "local iphi"
        log = False
        if args.log:
            log = True
        make_plot_2d(hist, hist_name, pathname, xlab, ylab, log)

    for h, hist in enumerate(hists1d):
        hist_name = args.hists1D[h]
        xlab = ""
        ylab = "a.u."
        if "timeSig" in hist_name:
            xlab = "Photon time significance"
        make_plot_1d(hist, hist_name, pathname, xlab, ylab, args.log, args.wide)


if __name__ == "__main__":
    main()
