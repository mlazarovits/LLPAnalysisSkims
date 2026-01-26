import argparse
from PlotFormatter import PlotFormatter, PlotFormatHelper
from ROOT import TFile


#for file processing, etc
if __name__ == "__main__":
    obs = ["MsRs","Ms","Rs","yields"]
    argobs = obs + ["all"]
    parser = argparse.ArgumentParser(description="format RJR plots")
    parser.add_argument("--inputFile","-i",help="file with histograms for formatting",required=True)
    parser.add_argument("--ofileName","-o",help="name for output file")
    parser.add_argument("--format","-f",help="what to save plots as",default="png")
    parser.add_argument("--obs",nargs='+',help='which obs/plots to make',choices=argobs)
    args = parser.parse_args()   
    plot_format = "."+args.format

    if args.ofileName is not None:
        ofilename = args.ofileName
    else: 
        ofilename = args.inputFile[:args.inputFile.find(".root")]
        ofilename += "_formatted.root"
    ofile = TFile(ofilename,"RECREATE");
    helper = PlotFormatHelper(args.inputFile)
    formatter = PlotFormatter()

    if 'all' in args.obs:
        args.obs = obs

    data_procs = ["METPD"]
    mc_procs = ["GJets","QCD","WJets","ZJets","TTXJets","DiPJBox"]

   
    #do data-only plots
    for proc in data_procs:
        #make 2D Ms-Rs plots - nonIsoEECR (data and MC bkg), BHCR (data)
        channel = "1pho"
        regions = ["earlyBHCR","lateBHCR"]
        if 'MsRs' in args.obs:
            #iterate through and format each
            for region in regions:
                msrs_hist = helper.GetHists(process=proc,obs="MsRs",channel=channel,region=region,histtype="2d")[0]
                hist_labels = [msrs_hist.GetName()]
                hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
                #remove sample, goes in sample_label
                group_labels = group_label.split(", ")
                sample_label = group_labels[0]
                group_label = ", ".join(group_labels[1:])
                print("group_label",group_label,"sample_label",sample_label)
                hist_tev = helper.GeVtoTeV(msrs_hist, axis = 0)
                canvas_name = "can_msrs_2d_"+proc+"_"+channel+"_"+region
                canvas = formatter.format_2d_hist(canvas_name, hist_tev, sample_label, "M_{S} [TeV]", 0., 5, "R_{S}", 0., 1., normalize=True,globallabel = group_label,sample_label_x_pos=0.6)
                #canvas.SaveAs(canvas_name+plot_format)
                ofile.cd()
                canvas.Write()
       
        if 'yields' in args.obs: 
            #nonIsoEECR
            yields_hists = helper.GetHists( process=proc,obs="yields", channel = "ge2pho", region="nonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "ge2pho", region="verynonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="nonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="verynonIsoEECR")
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
            canvas_name = "can_unrolled_msrs_nonIsoEECRs_"+proc
            hist_labels = [hist.GetName() for hist in unrolled_hists]
            hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
            canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label)       
            canvas.SaveAs(canvas_name+plot_format)
            ofile.cd()
            canvas.Write()
    
            #BH CRs for data only
            yields_hists = helper.GetHists( process=proc,obs="yields", channel = "ge2pho", region="earlyBHCR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "ge2pho", region="lateBHCR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="earlyBHCR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="lateBHCR")
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
            canvas_name = "can_unrolled_msrs_BHCRs_"+proc
            hist_labels = [hist.GetName() for hist in unrolled_hists]
            hist_labels, group_label = helper.MakeLegendLabels(hist_labels,option='long')
            canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label)       
            canvas.SaveAs(canvas_name+plot_format)
            ofile.cd()
            canvas.Write()

   
    regions = ["nonIsoEECR","verynonIsoEECR"] 
    channel = "1pho"
    if "Ms" in args.obs:
        for region in regions:
            #make 1D Ms plots, with processes overlaid
            ms_hist_data = helper.GetHists(process="METPD", obs="Ms", channel=channel, region=region,histtype="stack")[0]
            ms_hists_mc = []
            for proc in mc_procs:
                ms_hists_mc += helper.GetHists(process=proc, obs="Ms", channel=channel, region=region,histtype="stack")
            hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
            hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
            ms_hists_mc_tev = []
            for hist in ms_hists_mc:
                ms_hists_mc_tev.append(helper.GeVtoTeV(hist))
            ms_hist_data_tev = helper.GeVtoTeV(ms_hist_data)
            canvas_name = "can_ms_datamccomp_"+region
            x_min = 0.
            x_max = 3
            canvas = formatter.format_stack_hists_datamc(canvas_name, ms_hist_data_tev, ms_hists_mc_tev, "M_{S} [TeV]",x_min,x_max, normalize = True, globallabel=group_label) 
            #canvas.SaveAs(canvas_name+plot_format)
            ofile.cd()
            canvas.Write()

    channel = "ge2pho"
    if "Rs" in args.obs:
        for region in regions:
            #make 1D Rs plots, with processes overlaid
            ms_hist_data = helper.GetHists(process="METPD", obs="Rs", channel=channel, region=region,histtype="stack")[0]
            ms_hists_mc = []
            for proc in mc_procs:
                ms_hists_mc += helper.GetHists(process=proc, obs="Rs", channel=channel, region=region,histtype="stack")
            canvas_name = "can_rs_datamccomp_"+region
            hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
            hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
            x_min = 0.
            x_max = 1.
            canvas = formatter.format_stack_hists_datamc(canvas_name, ms_hist_data, ms_hists_mc, "R_{S}",x_min,x_max, normalize = True, globallabel=group_label) 
            #canvas.SaveAs(canvas_name+plot_format)
            ofile.cd()
            canvas.Write()

print("Wrote canvases to",ofilename)


