import re

class EfficiencyParser:
    def __init__(self):
    	return
    
    
    def add_dollar_to_inequalities(self, s):
        """
        Finds all simple inequalities in a string and wraps them in $...$.
        Example: "rjrPTS[0] < 150" -> "rjrPTS[0] $< 150$"
        """
        # pattern: operator optionally preceded by whitespace, then a number or variable
        pattern = r"(\s*(<=|>=|<|>|==|!=)\s*[^&|]+)"
        
        # replace matches with $...$
        result = re.sub(pattern, lambda m: f" ${m.group(1).strip()}$ ", s)
        
        # clean extra spaces
        result = re.sub(r"\s+", " ", result).strip()
        
        return result
    
    
    
    def write_latex_table(self, output_path, selected_data, lumi, total_eff = None):
        with open(output_path, "w") as f:
            for infile, data in selected_data.items():
                f.write("\\begin{table}\n")
                f.write("\\centering\n")
                f.write("\\caption{"+infile+" weighted to $\\mathcal{L} = "+str(lumi)+" \\text{fb}^{-1}$}\n")
                f.write("\\begin{tabular}{l c c}\n")
                f.write("\\hline\n")
                f.write("Cut & Weighted Entries & Efficiency (\\%)\\\\\n")
                f.write("\\hline\n")
                for label, entries, eff, wt_pass in data:
                    if ">" in label or "<" in label:
                        label = self.add_dollar_to_inequalities(label)
                    writestr = f"{label}"
                    if label != "presel" and label != "ge1KUBaseLinePhoton" and "Early" not in label and "Late" not in label and "Prompt" not in label:
                        writestr = "\\textbf{"+writestr+"}"
                    #writestr += f" & {int(entries)} & {wt_pass:.3f} & {eff:.3f} \\\\\n"
                    writestr += f" & {wt_pass:.3f} & {eff:.3f} \\\\\n"
                    f.write(writestr)
                if(total_eff is not None):
                	f.write("\\hline\n")
                	f.write(f"Total efficiency over all regions &  & {total_eff:.3f} \\\\\n")
                f.write("\\hline\n")
                f.write("\\end{tabular}\n")
                f.write("\\end{table}\n\n")
    
    def report2str(self, report, wt_counts_map = None):
        begin = report.begin()
        if begin == report.end(): return ""
        allEntries = begin.GetAll()
        result = []
        for ci in report:
            name = ci.GetName().strip()
            pass_val = ci.GetPass()
            all = ci.GetAll()
            eff = ci.GetEff()
            cumulativeEff = 100.0 * float(pass_val) / float(allEntries) if allEntries > 0 else 0.0
            wt_pass = -1
            if wt_counts_map is not None:
                if name in wt_counts_map.keys():
                    wt_pass = wt_counts_map[name].GetValue()
            result+=[f"{name:10}: pass={pass_val:<10} all={all:<10} weighted={wt_pass:<10} -- eff={eff:.2f} % cumulative eff={cumulativeEff:.2f} %"]
        return result
    
    def get_denom_line(self, lines, name):
        for line in lines:
            linename, info = line.split(":",1)
            linename = linename.strip()
            if name == linename:
                info = info.strip().split()
                return info 
    
    def parse_eff_line(self, line,denom = None):
        # line is like: "cutName: 1234 (0.567)"
        # split by ':' first
        if ':' not in line:
            return None
        name, rest = line.split(':', 1)
        name = name.strip()
        # rest has entries and efficiency
        parts = rest.strip().split()
        if len(parts) < 2:
            return None
        pass_entries = parts[0]
        pass_entries = pass_entries[pass_entries.find("=")+1:]
        n_entries = int(pass_entries)
        wt_pass = parts[2]
        wt_pass = wt_pass[wt_pass.find("=")+1:]
        wt_pass = float(wt_pass)
        if denom is None:
            eff_perc = parts[4]
            eff_perc = eff_perc[eff_perc.find("=")+1:]
            eff = float(eff_perc)
        else: #denom is a line already processed
            denom = denom[0]
            denom = float(denom[denom.find("=")+1:])
            if denom < n_entries: #don't parse cuts that happened *before* denom cut
                return None
            eff = (float(n_entries) / denom)*100
 
        #for latex 
        name = name.replace("_","\_")
        return name, n_entries, eff, wt_pass
