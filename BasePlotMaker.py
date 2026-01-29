import ROOT
from typing import Dict, List, Tuple, Optional
import numpy as np
import cmsstyle as CMS

class UnrollMaker:
    def __init__(self):
        # Group label configuration
        self.label_config = {
            'group_label_size': 0.038,
            'group_label_y_position': 0.87,  # Near top of canvas
            'separator_line_width': 2,
            'separator_line_color': ROOT.kBlack
        }
        self.comparison_colors = [
            ROOT.kBlue + 1,
            ROOT.kGreen + 2,
            ROOT.kRed - 4,
            ROOT.kOrange - 3,
            ROOT.kViolet - 4,
            ROOT.kBlack,
            ROOT.kAzure + 2,
            ROOT.kGray
        ]
        self.marker_styles = [20, 21, 22, 23, 47, 25, 26, 32, 33]  # Different marker styles
        # Canvas configuration
        self.unrolled_canvas_config = {
            'width': 1200,
            'height': 600,
            'bottom_margin': 0.15,
            'left_margin': 0.08,
            'right_margin': 0.28,
            'top_margin': 0.1
        }

    def _add_global_label(self, overlay_pad: ROOT.TPad, global_label: str, x_pos: float = 0.55, y_pos: float = 0.915):
        """Add SV label for data/MC canvas (separate from CMS labels)."""
        sv_text = ROOT.TLatex()
        sv_text.SetTextFont(42)
        sv_text.SetNDC()
        sv_text.SetTextSize(0.04)
        sv_text.SetTextAlign(31)  # Right align
        sv_text.DrawLatex(x_pos, y_pos, global_label)

    def create_base_canvas(self, name: str, title: str = "", 
                          use_log_y: bool = True, use_grid: bool = True) -> ROOT.TCanvas:
        """
        Create a base canvas with proper configuration.
        
        Args:
            name: Canvas name
            title: Canvas title (usually empty)
            use_log_y: Enable logarithmic y-axis
            use_grid: Enable grid lines
            
        Returns:
            Configured ROOT.TCanvas
        """
        # Create canvas
        canvas = ROOT.TCanvas(name, title, 
                            self.unrolled_canvas_config['width'], 
                            self.unrolled_canvas_config['height'])
        
        # Set margins (bottom, left, and right)
        canvas.SetTopMargin(self.unrolled_canvas_config['top_margin'])
        canvas.SetBottomMargin(self.unrolled_canvas_config['bottom_margin'])
        canvas.SetLeftMargin(self.unrolled_canvas_config['left_margin'])
        canvas.SetRightMargin(self.unrolled_canvas_config['right_margin'])
        
        # Set log scale and grid
        if use_log_y:
            canvas.SetLogy()
        if use_grid:
            canvas.SetGridx()
            canvas.SetGridy()
        
        return canvas
    def create_offset_graph(self, hist: ROOT.TH1D, offset_fraction: float, 
                           marker_style: int = 20, marker_size: float = 1.0, 
                           marker_color: int = ROOT.kBlack) -> ROOT.TGraphErrors:
        """
        Convert histogram to TGraph with x-offset markers for jittered plotting.
        
        Args:
            hist: Input histogram
            offset_fraction: Fraction of bin width to offset (between -0.5 and 0.5)
            marker_style: ROOT marker style
            marker_size: Marker size
            marker_color: Marker color
            
        Returns:
            TGraphErrors with offset x positions
        """
        n_bins = hist.GetNbinsX()
        x_vals = []
        y_vals = []
        x_errors = []
        y_errors = []
        
        for i in range(1, n_bins + 1):
            bin_center = hist.GetBinCenter(i)
            bin_content = hist.GetBinContent(i)
            bin_error = hist.GetBinError(i)
            bin_width = hist.GetBinWidth(i)
            
            # Apply offset as fraction of bin width
            x_offset = offset_fraction * bin_width
            x_position = bin_center + x_offset
            
            x_vals.append(x_position)
            y_vals.append(bin_content)
            x_errors.append(0.0)  # No x error bars
            y_errors.append(bin_error)
        
        # Create TGraphErrors
        graph = ROOT.TGraphErrors(len(x_vals), 
                                np.array(x_vals, dtype=float),
                                np.array(y_vals, dtype=float),
                                np.array(x_errors, dtype=float), 
                                np.array(y_errors, dtype=float))
        
        # Style the graph
        graph.SetMarkerStyle(marker_style)
        graph.SetMarkerSize(marker_size)
        graph.SetMarkerColor(marker_color)
        graph.SetLineColor(marker_color)
        graph.SetLineWidth(2)
        
        return graph
    def add_separator_lines(self, canvas: ROOT.TCanvas, hist: ROOT.TH1D, binning_scheme: dict):
        """
         Add separator lines between groups.
         Args:
            canvas: Canvas to add lines to
            hist: Histogram for range information
            binning_scheme: The name of the binning scheme ('legacy_9bin' or 'merged_6bin')
        Returns:
            List of line objects (for memory management)
         """
        separator_bins = binning_scheme['separator_bins']
        total_bins = binning_scheme['total_bins']
                                  
        canvas.cd()
        canvas.Update()
        
        # === Create an overlay pad (fully visual coordinate system, NDC) ===
        overlay = ROOT.TPad("overlay","overlay",0,0,1,1)
        overlay.SetFillStyle(0)
        overlay.SetFrameFillStyle(0)
        overlay.SetBorderSize(0)
        overlay.SetBorderMode(0)
        overlay.SetMargin(0,0,0,0)
        overlay.SetBit(ROOT.kCannotPick)  # Make transparent to mouse events
        overlay.Draw()
        overlay.cd()
        
        # === Convert x positions to NDC within the *main* pad ===
        main = canvas.GetPad(0)
        main.Update()
        
        x_axis = hist.GetXaxis()
        x_min = x_axis.GetXmin()
        x_max = x_axis.GetXmax()
        
        # Pad margins
        left_ndc = main.GetLeftMargin()
        right_ndc = 1.0 - main.GetRightMargin()
        data_ndc_width = right_ndc - left_ndc
        
        # If the histogram has a custom binning (e.g. from GetXaxis().SetBinLabel),
        # its internal range might still be 0 to N-1.
        # We need to map this to the visual range for proper NDC conversion.
        # Assuming histogram is set up with 0 to total_bins-1, and each bin is unit width.
        hist_axis_range = total_bins
        
        is_logx = bool(main.GetLogx())
        import math

        def normalized_x(xval):
            # For unrolled plots, xval represents bin index from 0 to total_bins-1
            # We want to place separators at edges *between* bins
            # So, normalize based on total_bins
            if is_logx: # Logarithmic X axis is unlikely for categorical unrolled plots, but keep for robustness
                # This case might need specific handling if categorical bins were log-scaled
                return (math.log(xval) - math.log(x_min)) / (math.log(x_max)-math.log(x_min))
            else:
                return xval / hist_axis_range

        def x_to_ndc_from_bin_edge(bin_edge_index):
            # The bin_edge_index refers to the numerical position (e.g., 3.0 for between bin 2 and 3)
            norm = normalized_x(bin_edge_index)
            return left_ndc + norm * data_ndc_width

        # === Now draw the vertical lines in PURE NDC ===
        y_bottom = 0.1   # constant, visual extent ONLY
        y_top    = 0.9
        
        lines = []
        for b_edge in separator_bins:
            x_ndc = x_to_ndc_from_bin_edge(b_edge)
            line = ROOT.TLine()
            line.SetNDC(True)
            line.SetLineColor(self.label_config['separator_line_color'])
            line.SetLineWidth(self.label_config['separator_line_width'])
            line.SetLineStyle(1)
            line.DrawLine(x_ndc, y_bottom, x_ndc, y_top)
            lines.append(line)

        canvas.Modified()
        canvas.Update()
        return lines
    
    def add_group_labels(self, canvas: ROOT.TCanvas, group_labels: List[str], binning_scheme: dict) -> List[ROOT.TLatex]:
        """
        Add group labels at the top of the plot.
        
        Args:
            canvas: Canvas to add labels to
            group_labels: List of group label strings
            binning_scheme: The name of the binning scheme ('legacy_9bin' or 'merged_6bin')
            
        Returns:
            List of TLatex objects (for memory management)
        """
        group_widths = binning_scheme['group_widths']
        total_bins = binning_scheme['total_bins']

        # Draw on overlay pad 
        overlay_pad = canvas.GetListOfPrimitives().FindObject("overlay")
        overlay_pad.cd()
        
        # Get actual pad boundaries in NDC from main pad
        main_pad = canvas.GetPad(0)
        pad_left = main_pad.GetLeftMargin()
        pad_right = 1.0 - main_pad.GetRightMargin()
        plot_area_ndc_width = pad_right - pad_left
        
        text_objects = []
        latex = ROOT.TLatex()
        latex.SetTextAlign(22)  # Center alignment
        latex.SetTextSize(self.label_config['group_label_size'])
        latex.SetTextFont(42)  # Helvetica (normal, not bold)
        latex.SetNDC(True)  # Use NDC coordinates
        
        current_bin_edge = 0
        for i, group_name in enumerate(group_labels):
            if i >= len(group_widths):
                # This might happen if group_labels has more elements than group_widths (e.g. for legacy scheme with ms/rs grouping)
                # In that case, we fall back to a simple even distribution for the remaining labels
                # Or it indicates an inconsistency between group_labels and group_widths.
                # For robustness, we'll try to estimate remaining group width, but a warning might be useful.
                print(f"Warning: More group_labels than defined group_widths for scheme {binning_scheme}. Distributing evenly.")
                remaining_bins = total_bins - current_bin_edge
                if (len(group_labels) - i) > 0:
                    group_actual_width = remaining_bins / (len(group_labels) - i)
                else:
                    group_actual_width = 0 # Should not happen if check above is true
            else:
                group_actual_width = group_widths[i]

            # Calculate center of each group in NDC relative to the plot area
            group_center_relative_to_plot_area = (current_bin_edge + group_actual_width / 2.0) / total_bins
            group_center_ndc = pad_left + group_center_relative_to_plot_area * plot_area_ndc_width
            
            y_ndc = binning_scheme['group_label_y_position']
            
            text_obj = latex.DrawLatex(group_center_ndc, y_ndc, group_name)
            text_objects.append(text_obj)

            current_bin_edge += group_actual_width # Advance the bin edge for the next group
        
        return text_objects

