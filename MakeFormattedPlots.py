import argparse
from PlotFormatter import PlotFormatter, PlotFormatHelper

#for file processing, etc
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="format RJR plots")
    parser.add_argument("--inputFile","-i",help="file with histograms for formatting",required=True)
    parser.add_argument("--ofileName","-o",help="name for output file")
    parser.add_argument("--format","-f",help="what to save plots as",default="png")
    args = parser.parse_args()   
    plot_format = "."+args.format

    if args.ofileName is not None:
        ofilename = args.ofileName
    else: 
        ofilename = args.inputFile[:args.inputFile.find(".root")]
        ofilename += "_formatted.root"
  
     
    helper = PlotFormatHelper(args.inputFile)
    formatter = PlotFormatter()

    #make yield Ms-Rs plots
    #nonIsoEECR
    yields_hists = helper.GetHists(obs="yields", channel = "ge2pho", region="nonIsoEECR")
    yields_hists += helper.GetHists(obs="yields",channel = "ge2pho", region="verynonIsoEECR")
    yields_hists += helper.GetHists(obs="yields",channel = "1pho", region="nonIsoEECR")
    yields_hists += helper.GetHists(obs="yields",channel = "1pho", region="verynonIsoEECR")
    unrolled_hists = []
    for hist in yields_hists:
        print(hist.GetName())
        #format here
        #unrolls st x-ax are x-bins and y-ax are groups of x-bins
        #gevtotev tells which axis to transform from gev to tev for bin labels (0 = x, 1 = y)
        unrolled_msrs_hist, totalbins = helper.UnrollHist(hist, normalize = True, gevtotev = 0)
        unrolled_hists.append(unrolled_msrs_hist)
   
    #make 1D unrolled normalized shape Ms-Rs plots
    canvas_name = "can_unrolled_msrs_nonIsoEECRs"
    hist_labels = [hist.GetName() for hist in unrolled_hists]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label)        
    canvas.SaveAs(canvas_name+plot_format)
    
    #BH CRs
    yields_hists = helper.GetHists(obs="yields", channel = "ge2pho", region="earlyBHCR")
    yields_hists += helper.GetHists(obs="yields",channel = "ge2pho", region="lateBHCR")
    yields_hists += helper.GetHists(obs="yields",channel = "1pho", region="earlyBHCR")
    yields_hists += helper.GetHists(obs="yields",channel = "1pho", region="lateBHCR")
    unrolled_hists = []
    for hist in yields_hists:
        print(hist.GetName())
        #format here
        #unrolls st x-ax are x-bins and y-ax are groups of x-bins
        #gevtotev tells which axis to transform from gev to tev for bin labels (0 = x, 1 = y)
        unrolled_msrs_hist, totalbins = helper.UnrollHist(hist, normalize = True, gevtotev = 0)
        unrolled_hists.append(unrolled_msrs_hist)
   
    #make 1D unrolled normalized shape Ms-Rs plots
    canvas_name = "can_unrolled_msrs_BHCRs"
    hist_labels = [hist.GetName() for hist in unrolled_hists]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label)        
    canvas.SaveAs(canvas_name+plot_format)


    

    #make 2D Ms-Rs plots
    msrs_hists = helper.GetHists(obs="MsRs")
    #iterate through and format each
    #make 1D Ms plots, with processes overlaid
    procs = None #can specify which processes to overlay
    ch = "1pho"
    reg = "nonIsoEECR"
    ms_hists = helper.GetHists("Ms",procs, ch, reg) 
    #overlay each process and format
    
    #make 1D Rs plots, with processes overlaid
    rs_hists = helper.GetHists("Rs",procs, ch, reg) 
    #overlay each process and format
    




    
    #make 1D time sig plots, with processes overlaid
