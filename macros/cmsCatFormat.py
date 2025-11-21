import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import uproot
import argparse
import os
import re
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal.*") #suppress warning about 0 being the smallest subnormal number (not relevant for plotting)

def get_sms_label_from_file(file):
	#strip path
	filename = file.split("/")[-1]
	#strip extension
	filename = file.split(".")[0]
	#get mass info and ctau info only
	match = "_AODSIM_"
	filename = filename[filename.find(match)+len(match):filename.find("_rjrskim")]
	filename = filename.replace("_m","-m")
	filename = filename.replace("_c","-c")
	#gluino mass
	sms_name = filename.replace("mGl",r"$m_{\tilde{g}}$")
	##n2 mass
	sms_name = sms_name.replace("mN2",r"$m_{\tilde{\chi}^0_2}$")
	#n1 mass
	sms_name = sms_name.replace("mN1",r"$m_{\tilde{\chi}^0_1}$")
	sms_name = sms_name.replace("-ct",r", $c\tau$ = ")
	sms_name = sms_name.replace("p",".")
	return sms_name

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
    
    hep.cms.label(llabel="Preliminary",com="13")
    
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
    ax.set_ylabel(ylab, fontsize=20) #yaxis
    ax.set_xlabel(xlab, fontsize=20) #xaxis
    hep.cms.label(llabel="Preliminary",com="13")
   
    
 
    plottitle = pathname+"/"+inhist_name+".pdf"
    print("Saving histogram",plottitle)
    fig.savefig(plottitle)


def make_multi_plot_1d(inhists, inlabels, colors, inhist_name, pathname, xlab, ylab, log, wide):
    fig = None
    ax = None
    if wide:
        fig, ax = plt.subplots(figsize=(20,6.5))
    else:
        fig, ax = plt.subplots(figsize=(12,10))
    if log:
        ax.set_yscale('log')
    for h, hist in enumerate(inhists):
        hist.plot(ax = ax,label=inlabels[h], color = colors[h])
    ax.set_ylabel(ylab, fontsize=20) #yaxis
    ax.set_xlabel(xlab, fontsize=20) #xaxis
    plt.legend(loc='upper right')
    plt.ylim(1,5e6) 
    plt.xlim(-8,20) 
 
    hep.cms.label(llabel="Preliminary",com="13")
    
    plottitle = pathname+"/"+inhist_name+".pdf"
    print("Saving histogram",plottitle)
    fig.savefig(plottitle)


def main():
    hep.style.use("CMS")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputFile",'-i',help='input file with histograms to be formatted')
    parser.add_argument("--hists1D",help='1D hists to format',nargs="+",default=[])
    parser.add_argument("--hists2D",help='2D hists to format',nargs="+",default=[])
    parser.add_argument("--log",help='log zaxis (2D) or yaxis (1D)',default=False,action='store_true')
    parser.add_argument("--wide",help='create wide figure size',default=False,action='store_true')
    parser.add_argument("--multifile",help='plot multiple curves at once from given files',nargs='+')
    args = parser.parse_args()

    if(len(args.hists1D) < 1 and len(args.hists2D) < 0):
        print("No histograms specified.")
        return

    if args.inputFile is None:
        if args.multifile is None:
             print("Please provide input file.")
             exit()
        else:
             pathname = args.multifile[0].split("/")[-1]
             pathname = pathname.split(".")[0]
             pathname = pathname[:pathname.find("_AODSIM")]
             pathname = "formatted_plots/"+pathname
             if not os.path.exists(pathname):
                 os.makedirs(pathname)
    else:
         pathname = os.path.basename(args.inputFile)
         pathname = pathname[:pathname.rfind('.')]
         pathname = "formatted_plots/"+pathname
         if not os.path.exists(pathname):
             os.makedirs(pathname)




    if len(args.hists1D) == 1 and len(args.multifile) > 0:
        hist_name = args.hists1D[0]
        xlab = ""
        if "rjr_ms" in hist_name:
        	xlab = r"$M_{S}$ [GeV]"
        if "rjr_rs" in hist_name:
        	xlab = r"$R_{S}$ [GeV]"
        if "timesig" in hist_name:
                xlab = "Photon Time Significance"
        ylab = "a.u."
        hists = []
        labels = []
       	colors = [] 
        for file in args.multifile:
            hist = get_1d_plots_from_file(file, args.hists1D)[0]
            hists.append(hist)
            label = ""
            color = []
            label = get_sms_label_from_file(file)
            print("label",label)
            if not "SMS" in file:
                label = "Background"
                color = "Blue"
            if "mGl-2000_mN2-500_mN1-250_ct0p1" in file: 
                #label = r"$m_{\tilde{g}}$-2000-$m_{\tilde{\chi}^0_2}$-500-$m_{\tilde{\chi}^0_1}$-250, $c\tau$=0.1 m"
                color = "Magenta"
            if "mGl-2000_mN2-1000_mN1-500_ct0p1" in file: 
                #label = r"$m_{\tilde{g}}$-2000-$m_{\tilde{\chi}^0_2}$-1000-$m_{\tilde{\chi}^0_1}$-500, $c\tau$=0.1 m"
                color = "Green"
            if "mGl-2000_mN2-1950_mN1-1900_ct0p1" in file: 
                #label = r"$m_{\tilde{g}}$-2000-$m_{\tilde{\chi}^0_2}$-1950-$m_{\tilde{\chi}^0_1}$-1900, $c\tau$=0.1 m"
                color = "Orange"
            labels.append(label)
            colors.append(color)
        make_multi_plot_1d(hists, labels, colors, hist_name, pathname, xlab, ylab, args.log, args.wide)
        return

    #ylab = r"$\frac{E_{SC}}{p_{track}}$"
    hists2d = get_2d_plots_from_file(args.inputFile, args.hists2D)
    for h, hist in enumerate(hists2d):
        hist_name = args.hists2D[h]
        xlab = ""
        ylab = ""
        if "map" in hist_name or "Map" in hist_name:
            xlab = "local ieta"
            ylab = "local iphi"
        if "EtaCenter_TimeCenter" in hist_name:
            xlab = "time [ns]"
            ylab = "eta"
        if "subclRelGeoAvgVar_subclRelEnergy" in hist_name:
        	xlab = r"$\sqrt[3]{ \tilde{\sigma}^2_{time} \times \tilde{\sigma}^2_{\eta} \times \tilde{\sigma}^2_{\phi}}$"
        	ylab = r"$\tilde{E}$"
        log = False
        if args.log:
            log = True
        make_plot_2d(hist, hist_name, pathname, xlab, ylab, log)

    hists1d = get_1d_plots_from_file(args.inputFile, args.hists1D)
    for h, hist in enumerate(hists1d):
        hist_name = args.hists1D[h]
        xlab = ""
        ylab = "a.u."
        if "timeSig" in hist_name:
            xlab = "Photon time significance"
        make_plot_1d(hist, hist_name, pathname, xlab, ylab, args.log, args.wide)


if __name__ == "__main__":
    main()
