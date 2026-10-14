from ROOT import RDataFrame, TFile, gPad, TLatex, kFullCircle, TLine, kGray, kBlack
from FileProcessor import FileProcessor
import yaml
import cmsstyle
from PhotonIsoSFCalculator_BRReweight import extract_grid
#fastsim signal SFs for photon isolation 

def draw_stack_with_ratio(fastsim_hist, fullsim_hist, plot_info):
    """Make the plot with the histograms on the fly with a ratio panel.

    """
    cmsstyle.setCMSStyle()
    ret_decor = []
    fastsim_hist.SetLineStyle(2)
    if(plot_info['normalize']):
        fastsim_hist.Scale(1./fastsim_hist.Integral())
        fullsim_hist.Scale(1./fullsim_hist.Integral())

    hist_max = max(fastsim_hist.GetMaximum(),fullsim_hist.GetMaximum())
    ymax = hist_max*1.1 
    # set the global ROOT style to CMS style

    # tweak the number of divisions on the X axis
    # since the number of jets is an integer, disable secondary tick marks
    #cmsstyle.getCMSStyle().SetNdivisions(5, "X")

    # set the luminosity, the COM energy, the Run period to show in the canvases
    cmsstyle.SetLumi(None, run="Run 3")
    cmsstyle.SetEnergy(13.6)
    # default extra text is "Preliminary", set it to an empty string to remove it
    cmsstyle.SetExtraText('Preliminary Simulation')

    cmsstyle.getCMSStyle().SetTitleSize(0.05, "Y")
    cmsstyle.getCMSStyle().SetTitleSize(0.05, "X")
    cmsstyle.getCMSStyle().SetLabelSize(0.04, "Y")
    cmsstyle.getCMSStyle().SetLabelSize(0.04, "X")
    cmsstyle.getCMSStyle().SetTitleOffset(1.1, "Y")
    cmsstyle.getCMSStyle().SetTitleOffset(1.1, "X")
    # prepare the canvas, specifying a name, min and maxes on X and Y axes and axis titles.
    c = cmsstyle.cmsDiCanvas(plot_info['name'],plot_info['xmin'], plot_info['xmax'], plot_info['ymin'], 
                           ymax,
                           0.5, 1.5,
                           plot_info['label'],"a.u.", "Fullsim / Fastsim", iPos=0, square = False)


   
    c.cd(1)
    fullsim_hist.SetLineColor(602)
    fastsim_hist.SetLineColor(602)
    cmsstyle.cmsObjectDraw(fullsim_hist,"HIST")
    cmsstyle.cmsObjectDraw(fastsim_hist,"HIST")
    # prepare a legend and fill it
    # note that this plot has a slightly different size, so the legend bounds and the font size are different
    plotlegend = cmsstyle.cmsLeg(0.19,0.57,0.69,0.77, textSize=0.04)  # The legend!
    cmsstyle.addToLegend(plotlegend, (fullsim_hist, "fullsim","l"))
    cmsstyle.addToLegend(plotlegend, (fastsim_hist, "fastsim","l"))
    ret_decor.append(plotlegend)

    # draw the uncertainty band 
    #cmsstyle.cmsObjectDraw(htotal_prediction, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    # draw the data histogram 
    #cmsstyle.cmsObjectDraw(histos[0], "E", MarkerStyle=ROOT.kFullCircle)

    # add an extra label on the plot a position specified in normalized coordinates
    # make sure to set the proper font
    extraLabel = None
    if "extra_label" in plot_info.keys():
        extraLabel = TLatex(0.197, 0.82, plot_info["extra_label"])
        extraLabel.SetNDC()
        cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont)
        ret_decor.append(extraLabel)

    # force redraw axis
    cmsstyle.UpdatePad()

    c.cd(2)
    ratio = fullsim_hist.Clone()
    ratio.Divide(fastsim_hist)
    ratio.GetYaxis().SetTitleSize(0.01)
    ratio.GetYaxis().SetTitleOffset(0.48)
    ratio.SetLineColor(kBlack)
    #cmsstyle.cmsObjectDraw(prediction_ratio, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    cmsstyle.cmsObjectDraw(ratio, "E", MarkerStyle=kFullCircle)
    
    line = TLine(plot_info['xmin'], 1.0, plot_info['xmax'], 1.0)
    line.SetLineColor(kGray + 1) # Standard neutral color for CMS ratio references
    line.SetLineStyle(2)             # Dashed line style
    line.SetLineWidth(2)
    line.Draw("same")
    ret_decor.append(line)
    # Saving the result!
    cmsstyle.UpdatePad(c)
    return c, ret_decor 

#run on output of PhotonIsoSFCalculator_BRReweight.py
def get_hists(f, mgl, n2, mn1, ctau, sr = False):
    if sr:
        sel = "sr"
    else:
        sel = "pre"
    hfast_name = f"isoScore_fast_{sel}_mGl{mgl}_mN2{mn2}_mN1{mn1}_ct{ctau}"
    hfull_name = f"isoScore_full_{sel}_mGl{mgl}_mN2{mn2}_mN1{mn1}_ct{ctau}"
    if not f.GetListOfKeys().Contains(hfast_name) or not f.GetListOfKeys().Contains(hfull_name):
        print(f"Either fast or full sim hist for point mGl{mgl}_mN2{mn2}_mN1{mn1}_ct{ctau} with sel {sel} does not exist")
        return None
    hfast = f.Get(hfast_name)
    hfull = f.Get(hfull_name)

    return {"fast":hfast, "full":hfull}


def format_ctau(ctau):
    ret = ctau.replace("p",".")
    ret = float(ret)*100
    if ret > 1:
        ret = int(ret)
    return str(ret)


fileprocessor = FileProcessor()
grid = extract_grid(fileprocessor)

plot_info = {
        "normalize" : True,
        "xmin" : 0.2,
        "xmax" : 1.01,
        "ymax" : 0.2,
        "ymin" : 1e-3,
        "label": "Baseline Photon Isolation Score",
       } 
#grab plots from root file
root_file = "testPhotonIsoSF.root"
f = TFile.Open(root_file)
caninfos = []
ofilename = "fastfullcomp_photoniso.root"
ofile = TFile.Open(ofilename,"RECREATE")
for pt in grid:
    mgl = pt[0]
    mn2 = pt[1]
    mn1 = pt[2]
    ctau = pt[3]
    f.cd() 
    hists = get_hists(f,mgl,mn2,mn1,ctau)
    if hists is None:
        continue

    plot_info["name"] = f"baseLinePhoton_isoANNScore_FullFastRatio_mGl{mgl}_mN2{mn2}_mN1{mn1}_ct{ctau}"
    plot_info["extra_label"] = f"m_{{#tilde{{g}}}}({mgl})-m_{{#tilde{{#chi}}_{{2}}^{{0}}}}({mn2})-m_{{#tilde{{#chi}}_{{1}}^{{0}}}}({mn1}), c#tau={format_ctau(ctau)} cm"
    #create fastsim/fullsim plots with ratio bottom panel - do inclusively, for pt bins, for separate dMs, gluino masses, lifetimes, BRs
    can, decor = draw_stack_with_ratio(hists['fast'], hists['full'], plot_info)
    ofile.cd()
    can.Write(f"{plot_info['name']}_FullFastRatio".replace(" ", "_"))
print(f"Wrote plots to {ofilename}")
f.Close()
ofile.Close()

