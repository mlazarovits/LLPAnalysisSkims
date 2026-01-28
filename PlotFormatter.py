from LLPStandardPlots.src.plotter import Plotter1D, Plotter2D, PlotterDataMC
from LLPStandardPlots.src.style import StyleManager
from BasePlotMaker import UnrollMaker, ComparisonMaker 
import ROOT
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import uproot
import cmsstyle as CMS
import array
import numpy as np

ROOT.gStyle.SetOptStat(0000)
#make dictionary for label names for regions, channels, etc
def get_signal_label(inlabel):
    sig_model = inlabel[:inlabel.find("mGl")]
    mgl = inlabel[inlabel.find("mGl")+3:inlabel.find("mN2")]
    mn2 = inlabel[inlabel.find("mN2")+3:inlabel.find("mN1")]
    mn1 = inlabel[inlabel.find("mN1")+3:inlabel.find("ct")]
    ct = inlabel[inlabel.find("ct")+2:]
    ct = ct.replace("p",".")
    return f"m_{{#tilde{{g}}}}({mgl})-m_{{#tilde{{#chi}}_{{2}}^{{0}}}}({mn2})-m_{{#tilde{{#chi}}_{{1}}^{{0}}}}({mn1}), c#tau={ct} m"

def transform_to_final_state(inlabel):
    retlabel = ""
    labels = inlabel.split(", ")
    objlabel = [l for l in labels if "#gamma" in l]
    if len(objlabel) == 1: #1 channel
        retlabel = objlabel[0]
        reglabel = [l for l in labels if retlabel not in l]
        if len(reglabel) == 1: #1 region
            retlabel += f"self._labels_dict[reglabel[0]]"
            return retlabel
        else: #multi-region, not defined
            print("labels",labels,"parsing not defined")
            return inlabel
    else: #multi-channel, not defined
        print("labels",labels,"parsing not defined")
        return inlabel

def interval_to_label(interval):
    str_label = "["+str(interval[0])+", "+str(interval[1])
    if interval[1] == "inf":
        str_label += ")"
    else:
        str_label += "]"
    return str_label

class PlotFormatHelper:
    def __init__(self, filename):
        self._filename = filename
        self._file = ROOT.TFile.Open(self._filename)
        self._hist_index = None
        self._all_hists = None
        self.IndexHists()
        self._labels_dict = {
            "ge2pho" : "#geq 2 #gamma",
            "1pho" : "1 #gamma",
            "1pho1HadSV" : "1 #gamma + 1 SV_{hh}",
            "nonIsoEECR" : "Prompt EE nonIso CR",
            "verynonIsoEECR" : "Prompt EE veryNonIso CR",
            "BHCR_long" : "Beam Halo CR",
            "earlyBHCR_long" : "Early Beam Halo CR",
            "lateBHCR_long" : "Late Beam Halo CR",
            "isoEESR" : "Prompt Iso SR",
            "earlyPBCR" : "Early !Beam Halo CR",
            "latePBCR" : "Late !Beam Halo SR",
            "METPD" : "MET PD Run II",
            "METFullRunII" : "MET PD Run II",
            "MET18" : "MET PD 2018",
            "MET17" : "MET PD 2017",
            "MET16" : "MET PD 2016",
            "Prompt EE nonIso CR" : "_{t0}^{CR,!Iso}",
            "Prompt EE veryNonIso CR" : "_{t0}^{CR,V!Iso}",
            "earlyBHCR" : "_{t-}^{CR, BH}",
            "lateBHCR" : "_{t+}^{CR, BH}"
        }


    def IndexHists(self): #call once
        self._hist_index = defaultdict(set)
        self._all_hists = set()
        with uproot.open(self._filename) as f:
            for name, obj in f.items():
                if not hasattr(obj, "to_numpy"):  # TH1 / TH2 check
                    continue

                parts = name.split("_")
                if len(parts) != 4:
                    continue

                observable, process, channel, region = parts
                #clean up string
                region = region[:region.find(";")]
                self._all_hists.add(name)
                self._hist_index[("obs", observable)].add(name)
                self._hist_index[("process", process)].add(name)
                self._hist_index[("channel", channel)].add(name)
                self._hist_index[("region", region)].add(name)

    def GetHists(self, obs = None, process = None, channel = None, region = None, histtype = None):
        if obs is None and process is None and channel is None and region is None:
            print("Need to specify at least one of the following: obs, process, channel, region")
    
        selections = []
        if obs is not None:
            selections.append(self._hist_index[("obs", obs)])
        if process is not None:
            selections.append(self._hist_index[("process", process)])
        if channel is not None:
            selections.append(self._hist_index[("channel", channel)])
        if region is not None:
            selections.append(self._hist_index[("region", region)])
    
        names = (
            set.intersection(*selections)
            if selections else self._all_hists
        )
        ret_hists = [self._file[name] for name in names if self._file[name].Integral() > 0]
        if histtype is not None:
            for hist in ret_hists:
                name = hist.GetName()+"_"+histtype
                hist.SetName(name)
        return ret_hists 

    
    #normalize histogram
    def NormalizeHist(self, hist):
        if hist.Integral() > 0:
            hist.Scale(1/hist.Integral())

    def GeVtoTeV(self, hist,axis = None):
        if axis is None:
            hTeV = ROOT.TH1D(
                hist.GetName()+"tev",
                hist.GetTitle(),
                hist.GetNbinsX(),
                hist.GetXaxis().GetXmin() / 1000.0,
                hist.GetXaxis().GetXmax() / 1000.0
            )
            for i in range(1, hist.GetNbinsX()+1): 
                hTeV.SetBinContent(i, hist.GetBinContent(i))
                hTeV.SetBinError(i, hist.GetBinError(i))
            return hTeV
        else: #2d histogram
            if("TH2D" not in str(type(hist))):
                print("This functionality with choosing an axis to transform is only supported for TH2D. Do not pass axis option with TH1D")
                return hist
            xmin = hist.GetXaxis().GetXmin()
            xmax = hist.GetXaxis().GetXmax()
            ymin = hist.GetYaxis().GetXmin()
            ymax = hist.GetYaxis().GetXmax()
 
            nx = hist.GetNbinsX()
            ny = hist.GetNbinsY()
          
            if(axis == 0):
                xmin /= 1000.0
                xmax /= 1000.0
            elif(axis == 1):
                ymin /= 1000.0
                ymax /= 1000.0
            else:
                print("Axis",axis,"not supported for tev transformation")
                return hist
 
            hTeV = ROOT.TH2D(
                hist.GetName()+"tev",
                hist.GetTitle(),
                hist.GetNbinsX(),
                xmin,
                xmax,
                hist.GetNbinsY(),
                ymin,
                ymax
            )
            for ix in range(1, nx + 1):
                for iy in range(1, ny + 1):
                    hTeV.SetBinContent(ix, iy, hist.GetBinContent(ix, iy))
                    hTeV.SetBinError(ix, iy, hist.GetBinError(ix, iy))
            return hTeV

    #unroll 2D histogram
    def UnrollHist(self, h2, normalize = False, gevtotev = None):
        check = ("TH2" in str(type(h2)))
        if not check:
            print("This function only works for TH2 histograms. You have passed",type(h2),"Returning...")
            return
        if normalize:
            self.NormalizeHist(h2)
        xaxis = h2.GetXaxis()
        yaxis = h2.GetYaxis()
        
        nx = xaxis.GetNbins()
        ny = yaxis.GetNbins()
       
        xbins = [] 
        ybins = []
        totalbins = []
        for ix in range(1, nx + 1):
            xlow = xaxis.GetBinLowEdge(ix)
            xhi = xaxis.GetBinUpEdge(ix)
            xbins.append([xlow, xhi])

        for iy in range(1, ny + 1):
            ylow = yaxis.GetBinLowEdge(iy)
            yhi = yaxis.GetBinUpEdge(iy)
            if gevtotev == 1:
                ylow /= 1000
                yhi /= 1000
            ybins.append([ylow, yhi])
            totalxbins = []
            for ix in range(1, nx + 1):
                xlow = xaxis.GetBinLowEdge(ix)
                xhi = xaxis.GetBinUpEdge(ix)
                if gevtotev == 0:
                    xlow /= 1000
                    xhi /= 1000
                totalxbins.append([xlow, xhi])
            totalbins.append(([ylow, yhi], totalxbins))
        #create bin labels
        print("totalbins",totalbins, len(totalbins)) 
        
        nbins = nx*ny
        h1 = ROOT.TH1D(
            h2.GetName() + "_unrolled",
            h2.GetTitle(),
            nbins,
            0,
            nbins
        )
        h1.Sumw2()
        
        for iy in range(1, ny + 1):
            for ix in range(1, nx + 1):
                ibin = (iy - 1) * nx + ix
                h1.SetBinContent(ibin, h2.GetBinContent(ix, iy))
                h1.SetBinError(ibin,   h2.GetBinError(ix, iy))
        #set bin labels
        ax = h1.GetXaxis()
        ntotbin = 1
        for ylabel, xlabel in totalbins:
            for xl, xlabel in enumerate(xlabel):
                ax.SetBinLabel(ntotbin, interval_to_label(xlabel))
                ntotbin += 1
        return h1, totalbins 
        

    
    def MakeLegendLabels(self, inlabels, option=None):
        procs = set()
        chs = set()
        regs = set()
        label = inlabels[0]
        histtype = label[label.rfind("_")+1:] #remove histogram type tag
        for label in inlabels:
            print("label",label)
            obs, proc, ch, reg, histtype_l = label.split("_")
            if option is not None:
                reg += "_"+option
            if histtype_l != histtype:
                print("Error: combining hists with types",histtype_l,"and",histtype)
            if proc not in self._labels_dict.keys():
                procs.add(proc)
            else:
                procs.add(self._labels_dict[proc])
            if ch not in self._labels_dict.keys():
                chs.add(ch)
            else:
                chs.add(self._labels_dict[ch])
            if reg not in self._labels_dict.keys():
                regs.add(reg)
            else:
                regs.add(self._labels_dict[reg])
        list_types = [procs, chs, regs]
        final_labels = []
        group_label = []
        for ltype in list_types:
            if len(ltype) == 1:
                group_label.append(list(ltype)[0])
            elif len(ltype) > 1:
                if len(final_labels) == 0:
                    final_labels = [ll for ll in ltype]
                else:
                    new_final_labels = []
                    for flabel in final_labels:
                        for ll in ltype:
                            new_final_labels.append(flabel+" "+ll) 
                    final_labels = new_final_labels
            else: #len(l) == 0
                continue
        #remove repeated instances of labels
        group_label = np.unique(group_label).tolist() 
        print("final_labels",final_labels,"group_label",group_label)
        return final_labels, group_label

class PlotFormatter():
    def __init__(self, helper, lumi = 138, plot_format = ".png"):
        self._styler = StyleManager(lumi) 
        self._plotter1d = Plotter1D(self._styler)
        self._plotter2D = Plotter2D(self._styler)
        self._plotterDataMC = PlotterDataMC(self._styler)
        self._plotterDataMC._ensure_mc_colors()
        self._plot_format = plot_format
        self._lumi = lumi
        self._helper = helper
      

    def format_2d_hist(self, name, hist, sample_label, x_label, x_min, x_max, y_label, y_min, y_max, normalize = False, globallabel = "", sample_label_x_pos = 0.65):
        canvas = CMS.cmsCanvas(name, x_min, x_max, y_min, y_max, x_label, y_label, 
                              square=False, extraSpace=0.01, iPos=0, with_z_axis=True)
        axis_labels = {}
        axis_labels['x'] = x_label
        axis_labels['y'] = y_label
        final_state_label = ""
        print("sample_label",sample_label)
        if sample_label not in self._helper._labels_dict.values():
            if "SMS" in sample_label:
                sample_label = get_signal_label(sample_label)
                sample_label_x_pos = 0.32
            else:
                sample_label = self._helper._labels_dict[sample_label]
        if globallabel != "":
            final_state_label = transform_to_final_state(globallabel)
        self._plotter2D.plot_2d_baseFormat(hist, canvas, axis_labels, sample_label, final_state_label, sample_label_x_pos=sample_label_x_pos,normalize=normalize)
        return canvas
 
    def get_hist_process(self, hist): 
        procname = hist.GetName()
        procname = procname[procname.find("_")+1:]
        procname = procname[:procname.find("_")]
        return procname       
    
    def format_hist_1d_base(self, hist, color, style, normalize = False):
        hist.SetLineColor(color)
        hist.SetLineStyle(style)
        hist.SetLineWidth(2)
        hist.SetFillColor(color)
        hist.SetFillStyle(3003)
        hist.SetStats(0)
        if normalize and hist.Integral() > 0:
            hist.Scale(1.0 / hist.Integral())

    def format_1d_hist(self, name, sample_label, hist, var_label, color, style, xmin, xmax, normalize = False, logy=False, sample_label_x_pos = 0.65):
        self.format_hist_1d_base(hist, color, style, normalize=normalize)
        #do canvas
        can = self._plotter1d._initialize_canvas(name, xmin, xmax, var_label)	
        if logy:
            can.SetLogy()
        
        #setup hist axes
        x_label = var_label 
        self._plotter1d.setup_axes(hist, x_label, normalized=normalize) 
        hist.Draw("HIST")
        #add histograms after dereferencing
        proc = self.get_hist_process(hist)
        if "MET" not in proc:
        	prelim_str = "Simulation"
        else:
        	prelim_str = "Preliminary"
        self._styler.draw_cms_labels(prelim_str=prelim_str)#cms_x=0.16, cms_y=0.93, prelim_str="Preliminary", prelim_x=0.235, lumi_x=0.75, cms_text_size_mult=1.25)
        self._styler.draw_process_label(sample_label, x_pos=sample_label_x_pos, y_pos=0.88)
        return can

    def format_hists_1d(self, hists, colors, styles, normalize = False, logy = False):
        for hist, color, style in zip(hists, colors, styles):
            format_hist_1d_base(hist, color, style, normalize)
        #do canvas
        can = plotter._initialize_canvas(name, xmin, xmax, var_label)	
        if logy:
            can.SetLogy()

        #setup hist axes
        x_label = var_label 
        self._plotter1d.setup_axes(hists[0], x_label, normalized=normalize)
        hist[0].Draw("HIST")
        for hist in hists[1:]:
            hist.Draw("HIST SAME")

        #add histograms after dereferencing
        self._styler.draw_cms_labels(prelim_str="Preliminary")#cms_x=0.16, cms_y=0.93, prelim_str="Preliminary", prelim_x=0.235, lumi_x=0.75, cms_text_size_mult=1.25)
        self._styler.draw_process_label(sample_label, x_pos=sample_label_x_pos, y_pos=0.88)
        return can
    
    def add_histogram_to_canvas(self, canvas: ROOT.TCanvas, hist: ROOT.TH1D, 
                               draw_option: str = "hist", is_first: bool = True) -> None:
        """
        Add a histogram to the canvas with proper scaling.
        
        Args:
            canvas: Canvas to draw on
            hist: Histogram to add
            draw_option: ROOT draw option
            is_first: Whether this is the first histogram (sets axes)
        """
        canvas.cd()
        
        if is_first:
            # First histogram sets up axes
            hist.Draw(draw_option)
        else:
            # Subsequent histograms use "same"
            if "same" not in draw_option.lower():
                draw_option += " same"
            hist.Draw(draw_option)
        
        canvas.Update()

    def create_legend(self, entries: List[Dict], x1: float = 0.77, y1: float = 0.8,
                      x2: float = 1., y2: float = 0.88) -> ROOT.TLegend:
        """
        Create a legend for multiple histograms.
        
        Args:
            entries: List of dicts with 'object', 'label', 'option' keys
            x1, y1, x2, y2: Legend position in NDC coordinates
            
        Returns:
            ROOT.TLegend object
        """
        legend = ROOT.TLegend(x1, y1, x2, y2)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetMargin(0.2)  
        legend.SetEntrySeparation(0.1)
        
        # Use two columns if more than 3 entries
        #if len(entries) > 3:
        #legend.SetNColumns(2)
        #legend.SetColumnSeparation(-0.25)
        
        for entry in entries:
            legend.AddEntry(entry['object'], entry['label'], entry['option'])
        
        return legend
    def universal_cms_mark(self, cms_x: float, cms_y: float, text_size: float, 
                           preliminary_x: float = None, preliminary_y: float = None,
                           cms_text: str = "CMS", preliminary_text: str = "Preliminary", ) -> List[ROOT.TLatex]:
        """
        Universal CMS mark function with configurable positions and sizes.
        
        Args:
            cms_x: X position for CMS text (NDC coordinates)
            cms_y: Y position for CMS text (NDC coordinates)  
            text_size: Base text size (CMS text will be 1.3x larger)
            preliminary_x: X position for preliminary text (defaults to cms_x + 0.06)
            preliminary_y: Y position for preliminary text (defaults to cms_y)
            cms_text: CMS text (default: "CMS")
            preliminary_text: Preliminary text (default: "Preliminary")
            
        Returns:
            List of TLatex objects for memory management
        """
        if preliminary_x is None:
            preliminary_x = cms_x + 0.06
        if preliminary_y is None:
            preliminary_y = cms_y
            
        latex_objects = []
        
        # Draw CMS text (1.3x larger than base text size)
        cms_latex = ROOT.TLatex()
        cms_latex.SetNDC()
        cms_latex.SetTextAlign(11)  # Left bottom align
        cms_latex.SetTextFont(61)   # Bold font
        cms_latex.SetTextSize(text_size * 1.3)
        cms_latex.DrawLatex(cms_x, cms_y, cms_text)
        latex_objects.append(cms_latex)
        
        # Draw preliminary text (base text size)
        prelim_latex = ROOT.TLatex()
        prelim_latex.SetNDC() 
        prelim_latex.SetTextAlign(11)  # Left bottom align
        prelim_latex.SetTextFont(52)   # Italic font
        prelim_latex.SetTextSize(text_size)
        prelim_latex.DrawLatex(preliminary_x, preliminary_y, preliminary_text)
        latex_objects.append(prelim_latex)
        
        return latex_objects
    def add_cms_labels(self, canvas: ROOT.TCanvas,
                       x_location: float = 0.12,
                       y_location: float = 0.915,
                       text_size: float = 0.04,
                       lumi_location: float = 0.85) -> List[ROOT.TLatex]:
        """
        Add CMS preliminary mark and luminosity label.
        
        Args:
            canvas: Canvas to add labels to
            
        Returns:
            List of TLatex objects (for memory management)
        """
        # Draw on overlay pad
        overlay_pad = canvas.GetListOfPrimitives().FindObject("overlay")
        overlay_pad.cd()
        
        # Use universal CMS mark
        cms_objects = self.universal_cms_mark(x_location, y_location, text_size, preliminary_x=x_location+0.056)
        
        # Add luminosity label
        lumi_latex = ROOT.TLatex()
        lumi_latex.SetTextFont(42)
        lumi_latex.SetNDC()
        lumi_latex.SetTextSize(text_size)
        lumi_latex.SetTextAlign(31)  # Right align
        lumi_latex.DrawLatex(lumi_location, y_location, f"{self._lumi:.0f} fb^{{-1}} (13 TeV)")
        
        return cms_objects + [lumi_latex]

    def format_stack_mc_hist(self, hist, color):
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineStyle(1)
        hist.SetLineWidth(1)
        hist.SetFillColor(color) 
        hist.SetStats(0)

    def format_stack_data_hist(self, hist):
        hist.SetMarkerStyle(20)
        hist.SetLineColor(ROOT.kBlack)
        hist.SetMarkerSize(self._styler.data_marker_size)
        hist.SetLineWidth(self._styler.data_line_width)

    def format_stack_hists_datamc(self, canvas_name, data_histogram, mc_histograms, var_label, x_min, x_max, normalize, draw_mc_uncertainty = False, globallabel = ""):
        comp_maker = ComparisonMaker()
        # Use shared canvas setup
        canvas, pad1, pad2 = self._plotterDataMC._setup_comparison_canvas(canvas_name, x_min, x_max, var_label)
        # Setup main pad with standard data/MC grid
        pad1.cd() 
        pad1.SetGridx(True)
        pad1.SetGridy(True)
        pad1.SetLogy(True)
        pad1.SetLeftMargin(self._styler.margin_left+0.04)
        pad1.SetRightMargin(self._styler.margin_right_ratio)
        # Create and add MC uncertainty band using helper
        mc_histograms_labeled = []
        for hist in mc_histograms:
            print(hist.GetName(), hist.Integral())
            bgname = self.get_hist_process(hist)
            color = comp_maker._get_background_color_index(bgname)
            self.format_stack_mc_hist(hist, color)
            format_bg_name = comp_maker.label_mapping[bgname]
             
            mc_histograms_labeled.append((hist, format_bg_name))
        # Sort MC histograms by yield (ascending order)
        mc_histograms_labeled.sort(key=lambda x: x[0].Integral())
        
        # Create THStack for MC
        stack = ROOT.THStack("stack", "")
        for mc_hist, _ in mc_histograms_labeled:
            stack.Add(mc_hist)
            
        # Apply normalization after histograms are created but before drawing
        if normalize:
            # Get total MC integral for normalization
            total_mc_integral = 0
            for mc_hist, _ in mc_histograms_labeled:
                total_mc_integral += mc_hist.Integral()
            
            # Normalize each MC histogram by the total MC integral
            if total_mc_integral > 0:
                for mc_hist, _ in mc_histograms_labeled:
                    mc_hist.Scale(1.0 / total_mc_integral)
        
        
        if normalize:
            data_histogram.Scale(1/data_histogram.Integral())
        self.format_stack_data_hist(data_histogram)
        # Set axis ranges
        data_max = data_histogram.GetMaximum()
        stack_max = stack.GetMaximum()
        max_val = max(data_max, stack_max)
        
        #stack.GetXaxis().SetRangeUser(x_min,x_max)
        # Draw stack and data with grid behind histograms
        stack.Draw("HIST")
        #pad1.RedrawAxis("G")  # Draw grid lines behind everything  
        #stack.Draw("HIST SAME")  # Redraw histogram content on top of grid
        
        # Set axis ranges after all drawing is complete
        if normalize:
            # For normalized plots, use fixed range that works with log scale
            stack.SetMaximum(max_val * 5.)
            stack.SetMinimum(2e-4)
            stack.GetHistogram().GetYaxis().SetRangeUser(2e-4, max_val * 5.)
        else:
            # For regular plots, use the original scaling
            stack.SetMaximum(max_val * 10.)
            stack.SetMinimum(0.5)
            stack.GetHistogram().GetYaxis().SetRangeUser(0.5, max_val * 10)
        stack.GetXaxis().SetLabelSize(0)
        
        # Set y-axis title based on normalization
        if normalize:
            # Extract just the variable name without units for the normalized y-axis title
            import re
            clean_var = re.sub(r'\s*\[.*?\]', '', var_label).strip()
            y_axis_title = "normalized events"#f"#frac{{1}}{{N}}  #frac{{dN}}{{d({clean_var})}}"
        else:
            y_axis_title = "number of events"
        stack.GetYaxis().SetTitle(y_axis_title)
        stack.GetYaxis().SetTitleSize(0.06)
        stack.GetYaxis().SetTitleOffset(1.1)
        stack.GetYaxis().SetLabelSize(0.05)
        stack.GetYaxis().CenterTitle(True)
        
        mc_uncertainty = self._plotterDataMC._create_mc_uncertainty_band(mc_histograms_labeled)
        if mc_uncertainty and draw_mc_uncertainty:
            mc_uncertainty.Draw("E2 SAME")  # E2 = error band only (no markers)
        
        if data_histogram:
            data_histogram.Draw("PEX0 SAME")
        
        # Create legend using helper
        legend = self._plotterDataMC._create_standard_legend(data_histogram, mc_uncertainty, mc_histograms_labeled)
        legend.Draw()
        
        # Draw standard labels using helper
        self._styler.draw_cms_labels(cms_x=0.16, cms_y=0.93, prelim_str="Preliminary", prelim_x=0.235, lumi_x=0.75, cms_text_size_mult=1.25)
    
        if globallabel != "":
            final_state_label = transform_to_final_state(globallabel)
            self._plotterDataMC._draw_region_label(canvas, final_state_label, x_pos=0.45, y_pos=0.93, textsize=0.05, plot_type="datamc")
        # Draw bottom pad (ratio) only if data is available
        pad2.cd()
        pad2.SetGridx(True)
        pad2.SetGridy(True)
        pad2.SetLeftMargin(self._styler.margin_left+0.04)
        pad2.SetRightMargin(self._styler.margin_right_ratio)
        
        # Create total MC histogram for ratio
        total_mc_hist = mc_histograms_labeled[0][0].Clone("total_mc")
        total_mc_hist.Reset()
        for mc_hist, _ in mc_histograms_labeled:
            total_mc_hist.Add(mc_hist)
        
        # Create ratio histogram
        ratio_hist = data_histogram.Clone("ratio")
        ratio_hist.Divide(total_mc_hist)
        
        # Style ratio plot
        ratio_hist.SetLineColor(self._plotterDataMC.data_color)
        ratio_hist.SetLineWidth(self._styler.data_line_width)
        ratio_hist.SetMarkerColor(self._plotterDataMC.data_color)
        ratio_hist.SetMarkerStyle(20)
        ratio_hist.SetMarkerSize(self._styler.data_marker_size)

        ratio_hist.SetTitle("")
        ratio_hist.GetXaxis().SetTitle(var_label)
        ratio_hist.GetYaxis().SetTitle("#frac{data}{model}")
        ratio_hist.GetYaxis().SetRangeUser(0.5, 1.5)
        ratio_hist.GetXaxis().SetTitleSize(0.14)
        ratio_hist.GetYaxis().SetTitleSize(0.14)
        ratio_hist.GetXaxis().SetLabelSize(0.12)
        ratio_hist.GetYaxis().SetLabelSize(0.12)
        ratio_hist.GetYaxis().SetTitleOffset(0.45)
        ratio_hist.GetXaxis().SetTitleOffset(1.15)
        ratio_hist.GetYaxis().SetNdivisions(505)
        ratio_hist.GetXaxis().CenterTitle(True)
        ratio_hist.GetYaxis().CenterTitle(True)
        
        # Draw ratio histogram first to establish axis formatting
        ratio_hist.Draw("PEX0")
        
        # Create and draw MC uncertainty band for ratio plot
        mc_ratio_uncertainty = self._plotterDataMC._create_mc_ratio_uncertainty_band(total_mc_hist)
        if mc_ratio_uncertainty and draw_mc_uncertainty:
            mc_ratio_uncertainty.Draw("E2 SAME")
        
        # Reference line at 1
        x_min_ratio = ratio_hist.GetXaxis().GetXmin()
        x_max_ratio = ratio_hist.GetXaxis().GetXmax()
        line = ROOT.TLine(x_min_ratio, 1, x_max_ratio, 1)
        line.SetLineStyle(2)
        line.Draw()
        
        # Keep objects alive
        canvas.ratio_hist = ratio_hist
        canvas.total_mc_hist = total_mc_hist
        canvas.mc_ratio_uncertainty = mc_ratio_uncertainty
        canvas.line = line

        pad1.SetFixedAspectRatio()
        pad2.SetFixedAspectRatio()
            
        canvas.Modified()
        canvas.Update()
        # Keep objects alive
        canvas.pad1 = pad1
        canvas.pad2 = pad2
        canvas.stack = stack
        canvas.mc_histograms = mc_histograms_labeled
        canvas.mc_uncertainty = mc_uncertainty
        canvas.data_hist = data_histogram
        canvas.legend = legend
        
        return canvas


    def format_unrolled_hists(self, name, histograms, totalbins, xlabel, ylabel, labels, globallabel = ""):
        n_hists = len(histograms)
        # Create base canvas
        unroll_maker = UnrollMaker()
        canvas = unroll_maker.create_base_canvas(f"{name}_comparison_markers", use_grid=False)
        if n_hists != len(labels):
            print("Gave unroll maker",n_hists,"histograms but only",len(labels),"labels: ",labels)
            return canvas
        # Create overlay pad first, then draw grid lines, then redraw separator lines
        xbins_grouped = [xb[1] for xb in totalbins]
        xbins = [xxb for xb in xbins_grouped for xxb in xb] #flatten
        ybins = [yb[0] for yb in totalbins]
        nbins = len(xbins)
        binning_scheme = {}
        binning_scheme["separator_bins"] = []
        offset = 0
        for yybins, xxbins in totalbins[:-1]:
            binning_scheme["separator_bins"].append(len(xxbins) + offset)
            offset += len(xxbins)
        binning_scheme["total_bins"] = nbins
        binning_scheme["group_label_y_position"] = 0.87
        binning_scheme["group_widths"] = [len(xb) for xb in xbins_grouped]
        print("totalbins",totalbins,"separator_bins",binning_scheme["separator_bins"],"group_widths",binning_scheme["group_widths"])
       
        
        # Find maximum for proper scaling
        max_val = max(hist.GetMaximum() for hist in histograms)
        max_norm = max(hist.Integral() for hist in histograms)
        print("max_norm",max_norm)
        yaxistitle = "events"
        if(max_norm == 1):
            yaxistitle = "normalized events"

 
        # Create first histogram as invisible for axis setup
        axis_hist = histograms[0].Clone(f"{name}_axis_template")
        axis_hist.SetMaximum(max_val * 5.)  # 5x headroom for log scale
        axis_hist.SetMinimum(1e-2)  # 5x headroom for log scale
        axis_hist.SetLineColor(0)  # Invisible
        axis_hist.SetMarkerSize(0)  # No markers
        axis_hist.SetFillStyle(0)  # No fill
        #add individual xaxis label    
        axis_hist.GetXaxis().SetTitle(xlabel)
        axis_hist.GetYaxis().SetTitle(yaxistitle)
        axis_hist.GetXaxis().SetTitleSize(0.045)
        axis_hist.GetYaxis().SetTitleSize(0.045)
        axis_hist.GetXaxis().SetLabelSize(0.045)
        axis_hist.GetYaxis().SetLabelSize(0.045)
        axis_hist.GetYaxis().SetTitleOffset(0.86)
        axis_hist.GetXaxis().SetTitleOffset(1.04)
        axis_hist.GetXaxis().CenterTitle()
        axis_hist.GetYaxis().CenterTitle()
        self.add_histogram_to_canvas(canvas, axis_hist, "hist", is_first=True)
        
        # Ensure grid is visible after drawing axis
        #canvas.SetGridx()
        #canvas.SetGridy()
        #canvas.Update()

        # Sort histograms by total yield (highest to lowest for display order)
        hist_with_yields = []
        for i, hist in enumerate(histograms):
            # Use original yields if provided, otherwise fall back to histogram integral
            #if original_yields and i < len(original_yields):
            #    yield_total = original_yields[i]
            #else:
            #    yield_total = hist.Integral()
            yield_total = hist.Integral()
            hist_with_yields.append((yield_total, hist, labels[i], i))
        
        # Sort by yield (highest first)
        hist_with_yields.sort(key=lambda x: x[0], reverse=True)
        
        
        # Convert histograms to offset graphs and draw them
        graphs = []
        legend_entries = []
       
        for display_idx, (yield_total, hist, label, orig_idx) in enumerate(hist_with_yields):
            # Calculate offset: spread evenly across bin width
            # For n histograms, offsets go from -0.4 to +0.4 (leaving 20% margin on each side)
            if n_hists == 1:
                offset = 0.0
            else:
                offset = -0.4 + (0.8 * display_idx / (n_hists - 1))
            
            # Get histogram color for original index
            color = unroll_maker.comparison_colors[orig_idx % len(unroll_maker.comparison_colors)]
            marker_style = unroll_maker.marker_styles[orig_idx % len(unroll_maker.marker_styles)]
            
            # Create offset graph
            graph = unroll_maker.create_offset_graph(hist, offset, marker_style=marker_style, 
                                           marker_size=1.4, marker_color=color)
            
            graphs.append(graph)
            legend_entries.append({'object': graph, 'label': label, 'option': 'p'})
        
        # Draw all graphs
        canvas.cd()
        for graph in graphs:
            graph.Draw("P SAME")  # P = markers, SAME = on same canvas
        
        canvas.Update()
        
        # Create legend
        legend = self.create_legend(legend_entries, 0.71, 0.7, 1.01, 0.91)
        canvas.cd()
        legend.Draw()
        canvas.legend = legend
 
        separator_lines = unroll_maker.add_separator_lines(canvas, histograms[0], binning_scheme=binning_scheme)
        
        # Draw vertical grid lines on overlay pad
        overlay_pad = canvas.GetListOfPrimitives().FindObject("overlay")
        if overlay_pad:
            overlay_pad.cd()
            
            # Get axis ranges and convert to NDC coordinates for overlay
            hist = axis_hist
            x_min = hist.GetXaxis().GetXmin()
            x_max = hist.GetXaxis().GetXmax()
            
            # Get main pad margins to map to overlay NDC coordinates
            main_pad = canvas.GetPad(0)
            left_margin = main_pad.GetLeftMargin()
            right_margin = main_pad.GetRightMargin()
            bottom_margin = main_pad.GetBottomMargin()
            top_margin = main_pad.GetTopMargin()
            
            # Calculate plotting area in NDC
            plot_left = left_margin
            plot_right = 1.0 - right_margin
            plot_bottom = bottom_margin
            plot_top = 1.0 - top_margin
            
            # Get separator line positions to avoid drawing grid lines there
            separator_bins = binning_scheme["separator_bins"]
            
            # Draw vertical grid lines at bin edges, EXCEPT at plot edges and separator positions
            n_bins = hist.GetNbinsX()
            grid_lines = []
            for i in range(2, n_bins + 1):  # Start from bin 2, skip first and last edges
                if i == n_bins + 1:  # Skip last edge
                    continue
                if (i - 1) in separator_bins:  # Skip separator positions (use i-1 for indexing)
                    continue
                x_pos = hist.GetXaxis().GetBinLowEdge(i)
                
                # Convert x position to NDC fraction
                x_frac = (x_pos - x_min) / (x_max - x_min)
                x_ndc = plot_left + x_frac * (plot_right - plot_left)
                
                line = ROOT.TLine(x_ndc, plot_bottom, x_ndc, plot_top)
                line.SetLineStyle(3)  # Dotted
                line.SetLineColor(ROOT.kBlack)
                line.SetLineWidth(1)
                line.SetNDC(True)  # Use NDC coordinates
                line.Draw()
                grid_lines.append(line)
            
            canvas.grid_lines = grid_lines  # Store to prevent garbage collection
        
        # Now add other decorations on top of grid lines
        #draw grouped (yaxis) labels
        group_labels = []
        for l in ybins:
            print("ylabel",l)
            label = ylabel + " #in " +interval_to_label(l)
            group_labels.append(label)
        text_objects = unroll_maker.add_group_labels(canvas, group_labels, binning_scheme=binning_scheme)

        if globallabel != "":
            unroll_maker._add_global_label(overlay_pad, globallabel)
    
        # Add individual labels for merged schemes
        #if binning_scheme in ['merged_rs', 'merged_ms'] and histogram_maker:
        #    individual_labels = histogram_maker.data_processor.get_individual_labels(binning_scheme)
        #    if individual_labels:
        #        individual_text_objects = self.add_individual_labels(canvas, individual_labels, binning_scheme=binning_scheme)
        #        text_objects.extend(individual_text_objects)
        
        # Add centered labels directly with TLatex for merged schemes
        #if binning_scheme in ['merged_rs', 'merged_ms', 'reversed_rs', 'reversed_ms']:
        #    try:
        #        centered_text_objects = self.add_merged_centered_labels(canvas, binning_scheme)
        #        text_objects.extend(centered_text_objects)
        #    except Exception as e:
        #        print(f"DEBUG: Error adding centered labels: {e}")
        #        import traceback
        #        traceback.print_exc()
        
        cms_objects = self.add_cms_labels(canvas,x_location = 0.07, lumi_location = 0.72)
        
        # Store objects
        canvas.lines = separator_lines
        #canvas.text_objects = text_objects
        canvas.cms_objects = cms_objects
        canvas.graphs = graphs
        canvas.axis_hist = axis_hist
        
        # Keep horizontal grid from ROOT
        canvas.cd()
        canvas.SetGridy(1)
        canvas.Modified()
        canvas.Update()
        return canvas


