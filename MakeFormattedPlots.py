import argparse
from PlotFormatter import PlotFormatter, PlotFormatHelper
from ROOT import TFile

histtype_dict = {
    "Ms" : "stack", #data/mc comparisons
    "Rs" : "stack",
    "timesig" : "stack",
    "MsRs" : "2d",
    "RISRPtISR" : "2d",
    "yields" : "unrolled"
}

labels_dict = {
    "Ms" : "M_{S} [TeV]",
    "Rs" : "R_{S}",
}



def make_datamc_plot(ofile, helper, obs, data, mc, region, channel):
    if(obs not in histtype_dict):
        print("Observable",obs,"does not have a histtype defined. Not making plot.")
        return
    if(histtype[obs] != "stack"):
        print("Observable",obs,"does not histtype stack, has histtype", histtype_dict[obs]," Not making plot.")
        return
    #make 1D Ms plots, with processes overlaid
    hist_data = helper.GetHists(process=data, obs=obs, channel=channel, region=region,histtype=histtype_dict[obs])[0]
    hists_mc = []
    for proc in mc:
        hists_mc += helper.GetHists(process=proc, obs=obs, channel=channel, region=region,histtype=histtype_dict[obs])
    if len(ms_hists_mc) < 1:
        print("No MC hists for observable",obs,". Not writing stack histogram.")
        return
    hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    hists_mc_tev = []
    if obs == "Ms":
        for hist in ms_hists_mc:
            hists_mc_tev.append(helper.GeVtoTeV(hist))
        hist_data_tev = helper.GeVtoTeV(hist_data)
    else:
        hists_mc_tev = hists_mc
        hist_data_tev = hist_data
    canvas_name = "can_"+obs+"_datamccomp_"+region
    x_min = 0.
    x_max = 3
    canvas = formatter.format_stack_hists_datamc(canvas_name, hist_data_tev, hists_mc_tev, labels_dict[obs],x_min,x_max, normalize = True, globallabel=group_label)
    #canvas.SaveAs(canvas_name+plot_format)
    ofile.cd()
    canvas.Write()


#for file processing, etc
if __name__ == "__main__":
    obs = ["MsRs","Ms","Rs","yields","RISRPtISR","timesig"]
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
    formatter = PlotFormatter(helper)

    if 'all' in args.obs:
        args.obs = obs

    data_procs = ["METFullRunII"]
    mc_procs = ["GJets","QCD","WJets","ZJets","TTXJets","DiPJBox"]
    procs = args.inputFile[:args.inputFile.find("_rjrObs")] 
    procs = procs.split("_") 
 
    for proc in procs:
        #make 2D Ms-Rs plots - nonIsoEECR (data and MC bkg), BHCR (data)
        channels = ["1pho", "ge2pho"]
        regions = ["earlyBHCR","lateBHCR","loosenonIsoEECR"]
        if 'MsRs' in args.obs:
            #iterate through and format each
            for channel in channels:
                for region in regions:
                    if "SMS" in proc and "CR" in region:
                        continue
                    if region in channels and region != channel: #channels are mutually exclusive
                        continue
                    msrs_hist = helper.GetHists(process=proc,obs="MsRs",channel=channel,region=region,histtype="2d")
                    print("# hists for proc",proc,"channel",channel,"region",region,":",len(msrs_hist))
                    if len(msrs_hist) < 1:
                        continue
                    msrs_hist = msrs_hist[0]
                    hist_labels = [msrs_hist.GetName()]
                    hist_labels, group_labels = helper.MakeLegendLabels(hist_labels)
                    #remove sample, goes in sample_label
                    sample_label = group_labels[1]
                    group_labels.remove(sample_label)
                    group_label = ", ".join(group_labels)
                    print("group_label",group_label,"sample_label",sample_label)
                    #convert x-axis (Ms) to TeV
                    hist_tev = helper.GeVtoTeV(msrs_hist, axis = 0)
                    canvas_name = "can_msrs_2d_"+proc+"_"+channel+"_"+region
                    ms_max = 3.
                    rs_max = 1.
                    if "SMS" in proc:
                        normalize = True
                        helper.NormalizeHist(hist_tev)
                        formatter.SetLumi(0)
                    else:
                        normalize = False
                        formatter.SetLumi(138)
                    canvas = formatter.format_2d_hist(canvas_name, hist_tev, sample_label, "M_{S} [TeV]", 0., ms_max, "R_{S}", 0., rs_max, normalize=normalize,globallabel = group_label,sample_label_x_pos=0.6)
                    #canvas.SaveAs(canvas_name+plot_format)
                    ofile.cd()
                    canvas.Write()
        if 'RISRPtISR' in args.obs:
            #iterate through and format each
            for channel in channels:
                for region in regions:
                    if "SMS" in proc and "CR" in region:
                        continue
                    if region in channels and region != channel: #channels are mutually exclusive
                        continue
                    rsptisr_hist = helper.GetHists(process=proc,obs="RISRPtISR",channel=channel,region=region,histtype="2d")
                    #print("# hists for proc",proc,"channel",channel,"region",region,":",len(msrs_hist))
                    if len(rsptisr_hist) < 1:
                        continue
                    rsptisr_hist = rsptisr_hist[0]
                    hist_labels = [rsptisr_hist.GetName()]
                    hist_labels, group_labels = helper.MakeLegendLabels(hist_labels)
                    #remove sample, goes in sample_label
                    sample_label = group_labels[1]
                    group_labels.remove(sample_label)
                    group_label = ", ".join(group_labels)
                    print("group_label",group_label,"sample_label",sample_label)
                    canvas_name = "can_rsptisr_2d_"+proc+"_"+channel+"_"+region
                    if "SMS" in proc:
                        normalize = True
                        helper.NormalizeHist(rsptisr_hist)
                    else:
                        normalize = False
                    canvas = formatter.format_2d_hist(canvas_name, rsptisr_hist, sample_label, "R_{ISR}", 50., 1500, "p^{ISR}_{T} [GeV]", 0., 1., normalize=normalize,globallabel = group_label,sample_label_x_pos=0.6)
                    #canvas.SaveAs(canvas_name+plot_format)
                    ofile.cd()
                    canvas.Write()
       
        if 'yields' in args.obs: 
            #nonIsoEECR
            yields_hists = helper.GetHists( process=proc,obs="yields", channel = "ge2pho", region="looseNotTightnonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "ge2pho", region="tightnonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="looseNotTightnonIsoEECR")
            yields_hists += helper.GetHists(process=proc,obs="yields",channel = "1pho", region="tightnonIsoEECR")
            unrolled_hists = []
            if len(yields_hists) == 0:
                print("no nonIsoEECR yields hists found for process",proc)
                continue
            for hist in yields_hists:
                print("using yields hist",hist.GetName())
                #format here
                #unrolls st x-ax are x-bins and y-ax are groups of x-bins
                #gevtotev tells which axis to transform from gev to tev for bin labels (0 = x, 1 = y)
                unrolled_msrs_hist, totalbins = helper.UnrollHist(hist, normalize = True, gevtotev = 0)
                unrolled_hists.append(unrolled_msrs_hist)
   
            #make 1D unrolled normalized shape Ms-Rs plots for data only
            canvas_name = "can_unrolled_msrs_nonIsoEECRs_"+proc
            hist_labels = [hist.GetName() for hist in unrolled_hists]
            hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
            canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, globallabel=group_label[0]) 
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
            canvas = formatter.format_unrolled_hists(canvas_name, unrolled_hists, totalbins, "M_{S} [TeV]", "R_{S}", hist_labels, group_label[0])
            canvas.SaveAs(canvas_name+plot_format)
            ofile.cd()
            canvas.Write()

   
    regions = ["loosenonIsoEECR","looseNotTightnonIsoEECR", "tightnonIsoEECR"] 
    channels = ["1pho","ge2pho"]
    for channel in channels:
        for region in regions:
            if "Ms" in args.obs:
                #make 1D Ms plots, with processes overlaid
                ms_hist_data = helper.GetHists(process="METFullRunII", obs="Ms", channel=channel, region=region,histtype="stack")[0]
                ms_hists_mc = []
                for proc in mc_procs:
                    ms_hists_mc += helper.GetHists(process=proc, obs="Ms", channel=channel, region=region,histtype="stack")
                if len(ms_hists_mc) < 1:
                    print("No MC hists for observable Ms. Not writing stack histogram.")
                    break
                hist_labels = [hist.GetName() for hist in ms_hists_mc + [ms_hist_data]]
                hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
                ms_hists_mc_tev = []
                for hist in ms_hists_mc:
                    ms_hists_mc_tev.append(helper.GeVtoTeV(hist))
                ms_hist_data_tev = helper.GeVtoTeV(ms_hist_data)
                canvas_name = "can_ms_datamccomp_"+channel+"_"+region
                x_min = 0.
                x_max = 3
                canvas = formatter.format_stack_hists_datamc(canvas_name, ms_hist_data_tev, ms_hists_mc_tev, "M_{S} [TeV]",x_min,x_max, normalize = True, globallabel=group_label) 
                #canvas.SaveAs(canvas_name+plot_format)
                ofile.cd()
                canvas.Write()

            if "Rs" in args.obs:
                    #make 1D Rs plots, with processes overlaid
                    rs_hist_data = helper.GetHists(process="METFullRunII", obs="Rs", channel=channel, region=region,histtype="stack")[0]
                    rs_hists_mc = []
                    for proc in mc_procs:
                        rs_hists_mc += helper.GetHists(process=proc, obs="Rs", channel=channel, region=region,histtype="stack")
                    if len(rs_hists_mc) < 1:
                        print("No MC hists for observable Rs. Not writing stack histogram.")
                        break
                    canvas_name = "can_rs_datamccomp_"+channel+"_"+region
                    hist_labels = [hist.GetName() for hist in rs_hists_mc + [rs_hist_data]]
                    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
                    x_min = 0.
                    x_max = 1.
                    canvas = formatter.format_stack_hists_datamc(canvas_name, rs_hist_data, rs_hists_mc, "R_{S}",x_min,x_max, normalize = True, globallabel=group_label) 
                    #canvas.SaveAs(canvas_name+plot_format)
                    ofile.cd()
                    canvas.Write()
    
    #needs to be in kinematic sideband
    channel = "1pho"
    region = "MsCR"
    if "timesig" in args.obs:
        timesig_hist = helper.GetHists(process=proc, obs="timeSig", channel=channel, region=region,histtype="1d")[0]
        canvas_name = "can_timesig_"+proc+"_"+region
        hist_labels = [hist.GetName() for hist in [timesig_hist]]
        hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
        sample_label = group_label[1]
        xmin = 0.
        xmax = 1.
        color = 1179
        canvas = formatter.format_1d_hist(canvas_name, sample_label, timesig_hist, "S_{t}", color=color, style=1, xmin=xmin, xmax=xmax, normalize = True, logy=True) 
        #canvas.SaveAs(canvas_name+plot_format)
        ofile.cd()
        canvas.Write()

print("Wrote canvases to",ofilename)


