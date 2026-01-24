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
    
    #TODO
    #make 2D Ms-Rs plots - nonIsoEECR (data and MC bkg), BHCR (data)
    channel = "1pho"
    regions = ["earlyBHCR","lateBHCR"]
    #iterate through and format each
    for region in regions:
        msrs_hist = helper.GetHists(process="METPD",obs="MsRs",channel=channel,region=region,histtype="2d")[0]
        hist_labels = [msrs_hist.GetName()]
        hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
        #remove sample, goes in sample_label
        group_labels = group_label.split(", ")
        sample_label = group_labels[0]
        group_label = ", ".join(group_labels[1:])
        print("group_label",group_label,"sample_label",sample_label)
        hist_tev = helper.GeVtoTeV(msrs_hist, axis = 0)
        canvas_name = "can_msrs_2d_METPD_"+channel+"_"+region
        canvas = formatter.format_2d_hist(canvas_name, hist_tev, sample_label, "M_{S} [TeV]", 0., 5, "R_{S}", 0., 1., normalize=True,globallabel = group_label,sample_label_x_pos=0.6)
        canvas.SaveAs(canvas_name+plot_format)

    exit() 
    #make 1D Ms plots, with processes overlaid
    channel = "1pho"
    ms_hist_data = helper.GetHists(process="METPD", obs="Ms", channel=channel, region="nonIsoEECR",histtype="stack")[0]
    ms_hists_mc = helper.GetHists(process="GJets", obs="Ms", channel=channel, region="nonIsoEECR",histtype="stack") 
    ms_hists_mc += helper.GetHists(process="QCD", obs="Ms", channel=channel, region="nonIsoEECR",histtype="stack")
    hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    ms_hists_mc_tev = []
    for hist in ms_hists_mc:
        ms_hists_mc_tev.append(helper.GeVtoTeV(hist))
    ms_hist_data_tev = helper.GeVtoTeV(ms_hist_data)
    canvas_name = "can_ms_datamccomp_nonIsoEECRs"
    x_min = 0.
    x_max = 5
    canvas = formatter.format_stack_hists_datamc(canvas_name, ms_hist_data_tev, ms_hists_mc_tev, "M_{S} [TeV]",x_min,x_max, normalize = True, globallabel=group_label) 
    canvas.SaveAs(canvas_name+plot_format)

    #make 1D Rs plots, with processes overlaid
    channel = "ge2pho"
    ms_hist_data = helper.GetHists(process="METPD", obs="Rs", channel=channel, region="nonIsoEECR",histtype="stack")[0]
    ms_hists_mc = helper.GetHists(process="GJets", obs="Rs", channel=channel, region="nonIsoEECR",histtype="stack") 
    ms_hists_mc += helper.GetHists(process="QCD", obs="Rs", channel=channel, region="nonIsoEECR",histtype="stack")
    canvas_name = "can_rs_datamccomp_nonIsoEECRs"
    hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    x_min = 0.
    x_max = 1.
    canvas = formatter.format_stack_hists_datamc(canvas_name, ms_hist_data, ms_hists_mc, "R_{S}",x_min,x_max, normalize = True, globallabel=group_label) 
    canvas.SaveAs("test"+plot_format)
    


    
    #make 1D time sig plots, with processes overlaid

    exit()

    #make yield Ms-Rs plots for data only
    process = "METPD"
    #nonIsoEECR
    yields_hists = helper.GetHists( process=process,obs="yields", channel = "ge2pho", region="nonIsoEECR")
    yields_hists += helper.GetHists(process=process,obs="yields",channel = "ge2pho", region="verynonIsoEECR")
    yields_hists += helper.GetHists(process=process,obs="yields",channel = "1pho", region="nonIsoEECR")
    yields_hists += helper.GetHists(process=process,obs="yields",channel = "1pho", region="verynonIsoEECR")
    unrolled_hists = []
    if len(yields_hists) == 0:
        print("no nonIsoEECR yields hists found for process",process)
        exit()
    for hist in yields_hists:
        print(hist.GetName())
        #format here
        #unrolls st x-ax are x-bins and y-ax are groups of x-bins
        #gevtotev tells which axis to transform from gev to tev for bin labels (0 = x, 1 = y)
        unrolled_msrs_hist, totalbins = helper.UnrollHist(hist, normalize = True, gevtotev = 0)
        unrolled_hists.append(unrolled_msrs_hist)
   
    #make 1D unrolled normalized shape Ms-Rs plots for data only
    canvas_name = "can_unrolled_msrs_nonIsoEECRs"
    hist_labels = [hist.GetName() for hist in unrolled_hists]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label)        
    canvas.SaveAs(canvas_name+plot_format)
    
    #BH CRs for data only
    yields_hists = helper.GetHists( process="METPD",obs="yields", channel = "ge2pho", region="earlyBHCR")
    yields_hists += helper.GetHists(process="METPD",obs="yields",channel = "ge2pho", region="lateBHCR")
    yields_hists += helper.GetHists(process="METPD",obs="yields",channel = "1pho", region="earlyBHCR")
    yields_hists += helper.GetHists(process="METPD",obs="yields",channel = "1pho", region="lateBHCR")
    if len(yields_hists) == 0:
        print("no BH CR yields hists found for process",process)
        exit()
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


    
