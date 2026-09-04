import argparse
from PlotFormatter import PlotFormatter, PlotFormatHelper
from ROOT import TFile, kBlack, kGreen



def make_datamc_plot(ofile, helper, obs, data, mc, region):
    if(obs not in histtype_dict):
        print("Observable",obs,"does not have a histtype defined. Not making plot.")
        return
    if(histtype_dict[obs] != "stack"):
        print("Observable",obs,"does not histtype stack, has histtype", histtype_dict[obs]," Not making plot.")
        return
    channel = "presel"
    #make 1D Ms plots, with processes overlaid
    hist_data = helper.GetHists(process=data, obs=obs, channel=channel, region=region,histtype=histtype_dict[obs])
    if len(hist_data) > 1:
        hist_data = hist_data[0]
    else:
        print("No hists found for",data,obs,channel,region)
        exit()
    hists_mc = []
    for proc in mc:
        hists_mc += helper.GetHists(process=proc, obs=obs, channel=channel, region=region,histtype=histtype_dict[obs])
    if len(hists_mc) < 1:
        print("No MC hists for observable",obs,". Not writing stack histogram.")
        return
    hist_labels = [hist.GetName() for hist in hists_mc + [hist_data]]
    hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
    hists_mc_tev = []
    if obs == "Ms":
        for hist in hists_mc:
            hists_mc_tev.append(helper.GeVtoTeV(hist))
        hist_data_tev = helper.GeVtoTeV(hist_data)
    else:
        hists_mc_tev = hists_mc
        hist_data_tev = hist_data
    canvas_name = "can_"+obs+"_datamccomp_"+region
    x_min = 0.
    x_max = 3
    canvas = formatter.format_stack_hists_datamc(canvas_name, hist_data_tev, hists_mc_tev, labels_dict[obs],x_min,x_max, normalize = True, globallabel=group_label[0])
    #canvas.SaveAs(canvas_name+plot_format)
    ofile.cd()
    canvas.Write()


#for file processing, etc
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="format RJR plots")
    parser.add_argument("--inputFile","-i",help="file with histograms for formatting",required=True)
    parser.add_argument("--ofileName","-o",help="name for output file")
    #parser.add_argument("--format","-f",help="what to save plots as",default="png")
    parser.add_argument("--obs",nargs='+',help='which obs/plots to make')
    parser.add_argument("--regions",nargs="+",help='which regions to look at')
    parser.add_argument("--procs",nargs="+",help='which procs to look at')
    args = parser.parse_args()   
    #plot_format = "."+args.format

    lumi = 30
    com = 13.6

    if args.ofileName is not None:
        ofilename = args.ofileName
    else: 
        ofilename = args.inputFile[:args.inputFile.find(".root")]
    ofilename += "_formatted.root"
    ofile = TFile(ofilename,"RECREATE");
    helper = PlotFormatHelper(args.inputFile)
    formatter = PlotFormatter(helper)
    formatter.SetLumi(lumi)
    formatter.SetEnergy(com)

    if 'all' in args.obs:
        args.obs = obs

	
    for obs in args.obs:
        hists, colors, styles = [], [], []
        canvas_name = f"can_{obs}"
        for proc, region in zip(args.procs, args.regions):
            #make 1D Ms plots, with processes overlaid
            #ms_hist_data = helper.GetHists(process=data_procs[0], obs="Ms", channel=channel, region=region,histtype="stack")[0]
            print("proc",proc,"region",region)
            processes = helper.GetAllProcNames(proc)
            for fullproc in processes:
                hists += helper.GetHists(process=fullproc, obs=obs, region=region,histtype="stack")
            canvas_name += f"_{proc}{region}"
        print("# hists",len(hists))
        for hist in hists:
            if "METPD" in hist.GetName():
                colors.append(kBlack)
            else:
                colors.append(kGreen+len(colors)-5)
            styles.append(1)
        hist_labels = [hist.GetName() for hist in hists]
        hist_labels, group_label = helper.MakeLegendLabels(hist_labels)
        if len(group_label) > 0:
            globallabel = group_label[0]
        else:
            globallabel = ""
        canvas = formatter.format_hists_1d(canvas_name, hists, colors, styles,obs,obs, normalize = True, globallabel=globallabel, logy=True, fillstyle=0) 
        ofile.cd()
        canvas.Write()

print("Wrote canvases to",ofilename)


