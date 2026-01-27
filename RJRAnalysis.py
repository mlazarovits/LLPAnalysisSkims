from ROOT import RDataFrame, TChain, TFile, TH1, TH2, gInterpreter
import sys
import argparse
import os
from FileProcessor import FileProcessor
import re
# ROOT.EnableImplicitMT()  # optional

# -------------------------
# Utilities
# -------------------------

from array import array
TH1.SetDefaultSumw2(True)
TH2.SetDefaultSumw2(True)



def add_dollar_to_inequalities(s):
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



def write_latex_table(output_path, selected_data):
    with open(output_path, "w") as f:
        for infile, data in selected_data.items():
            f.write("\\begin{table}\n")
            f.write("\\centering\n")
            f.write("\\caption{"+infile+"}\n")
            f.write("\\begin{tabular}{l c c}\n")
            f.write("\\hline\n")
            f.write("Cut & Entries & Efficiency (\\%)\\\\\n")
            f.write("\\hline\n")
            for label, entries, eff in data:
                if ">" in label or "<" in label:
                    label = add_dollar_to_inequalities(label)
                f.write(f"{label} & {int(entries)} & {eff:.3f} \\\\\n")
            f.write("\\hline\n")
            f.write("\\end{tabular}\n")
            f.write("\\end{table}\n\n")

def report2str(report):
    begin = report.begin()
    if begin == report.end(): return ""
    allEntries = begin.GetAll()
    result = []
    for ci in report:
        name = ci.GetName()
        pass_val = ci.GetPass()
        all = ci.GetAll()
        eff = ci.GetEff()
        cumulativeEff = 100.0 * float(pass_val) / float(allEntries) if allEntries > 0 else 0.0
        result+=[f"{name:10}: pass={pass_val:<10} all={all:<10} -- eff={eff:.2f} % cumulative eff={cumulativeEff:.2f} %"]
    print(result)
    return result

def parse_eff_line(line):
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
    eff_perc = parts[3]
    eff_perc = eff_perc[eff_perc.find("=")+1:]
    eff = float(eff_perc)
  
    #for latex 
    name = name.replace("_","\_")
    return name, n_entries, eff

class RJRAnalysis:
    def __init__(self):
        self._branches = [
            "rjr_Rs","rjr_Ms","selCMet","rjrPTS",
            "evtFillWgt",
            "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
            "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight",
            "Trigger_PFMET120_PFMHT120_IDTight_PFHT60",
            "Trigger_PFMET120_PFMHT120_IDTight",
            "Flag_BadChargedCandidateFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_BadPFMuonFilter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter",
            "Flag_ecalBadCalibFilter",
            "Flag_eeBadScFilter",
            "Flag_goodVertices",
            "Flag_hfNoisyHitsFilter",
            "nSelPhotons",
            "selPhoWTime",
            "selPho_beamHaloCNNScore",
            "selPho_physBkgCNNScore",
            "selPhoWTimeSig",
            "selPho_nonIsoANNScore",
            "selPho_isoANNScore",
            "selPhoEta",
            "selPhoWTimeSig"
        ]
        self._threshs = {
            "bh": "0.917252",
            "pb": "0.81476355",
            "EE_nonIso": "0.9290591",
            "EE_veryNonIso": "0.9939665",
            "EE_iso" : "0.9994431", #80% efficiency, 5% bkg contanimination from SMS-GlGl ROC 
            "early_time": "-2",
            "late_time": "2",
            "prompt_time": "1",
            "early_timesig": "-2.5",
            "late_timesig": "2.5",
            "prompt_timesig": "2.5",
        }
        # --- Preselection Cuts ---
        self._metcut = "(selCMet > 150)"
        self._ptscut = "(rjrPTS[0] < 150)"
        self._triggers = (
            "(Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || "
            "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || "
            "Trigger_PFMET120_PFMHT120_IDTight_PFHT60 == 1 || "
            "Trigger_PFMET120_PFMHT120_IDTight == 1)"
        )
        self._met_filters = (
            "(Flag_BadPFMuonFilter == 1 && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && "
            "Flag_HBHENoiseFilter == 1 && "
            "Flag_HBHENoiseIsoFilter == 1 && "
            "Flag_ecalBadCalibFilter == 1 && "
            "Flag_eeBadScFilter == 1 && "
            "Flag_goodVertices == 1)"
        )
        #beam halo CR bins - limited stats regime
        self._msrs_bins = {}
        self._msrs_bins["BHCR"] = {"ms" : array("d", [700, 5000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._msrs_bins["PB"] = {"ms" : array("d", [700, 5000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._msrs_bins["nonIsoEE"] = {"ms": array("d", [500, 1000, 2000, 5000]),"rs" :array("d", [0.0, 0.15, 0.4, 1.0]) }
        self._msrs_bins["isoEE"] = {"ms": array("d", [500, 1000, 2000, 5000]),"rs" :array("d", [0.0, 0.15, 0.4, 1.0]) }

        
    def define_photon_quantities(self, df, thr):
        return (
            df
            .Define(
                "EEnonIsoScore",
                "selPho_nonIsoANNScore[(selPhoEta < -1.479) || (selPhoEta > 1.479)]"
            )
            .Define(
                "nEndcapPhotons",
                "selPhoEta[selPhoEta > 1.479 || selPhoEta < -1.479].size()"
            )
            .Define(
                "nEndcapPhotons_nonIso",
                f"selPho_nonIsoANNScore[(selPhoEta > 1.479 || selPhoEta < -1.479) && "
                f"(selPho_nonIsoANNScore > {self._threshs['EE_nonIso']})].size()"
            )
            .Define(
                "nEndcapPhotons_veryNonIso",
                f"EEnonIsoScore[EEnonIsoScore > {self._threshs['EE_veryNonIso']}].size()"
            )
            .Define(
                "nBHPhotons",
                f"selPho_beamHaloCNNScore[selPho_beamHaloCNNScore > {self._threshs['bh']}].size()"
            )
        )
    
    def apply_preselection(self, df):
        baseline = f"{self._metcut} && {self._ptscut} && {self._triggers} && {self._met_filters}"
        return df.Filter(baseline, "baseline")
    
    
    def print_photon_efficiencies(self, df):
        npho = df.Sum("nSelPhotons").GetValue()
        npho_ee = df.Sum("nEndcapPhotons").GetValue()
        npho_ee_noniso = df.Sum("nEndcapPhotons_nonIso").GetValue()
        npho_ee_vnoniso = df.Sum("nEndcapPhotons_veryNonIso").GetValue()
        npho_bh = df.Sum("nBHPhotons").GetValue()
    
        print("sum nSelPhotons:", npho)
        print("sum nEndcapPhotons:", npho_ee)
    
        if npho_ee > 0:
            print(
                "EE nonIso eff:",
                100 * npho_ee_noniso / npho_ee,
                "%  total:",
                100 * npho_ee_noniso / npho,
                "%"
            )
    
            print(
                "EE veryNonIso eff:",
                100 * npho_ee_vnoniso / npho_ee,
                "%  total:",
                100 * npho_ee_vnoniso / npho,
                "%"
            )
    
        print("BH eff:", 100 * npho_bh / npho, "%")
    
    def define_channels(self, df):
        return {
            "1pho": df.Filter("nSelPhotons == 1", "1nSelPho"),
            "ge2pho": df.Filter("nSelPhotons >= 2", "ge2nSelPho"),
        }
        '''
        df_1pho = df.Filter("nSelPhotons == 1", "1nSelPho")
        df_ge2pho = df.Filter("nSelPhotons >= 2", "ge2nSelPho")
      
        df_1pho = self.define_lead_photon_vars(df_1pho)
        df_ge2pho = self.define_lead_photon_vars(df_ge2pho)
 
        #define sublead branches for ge2pho channel 
        df_ge2pho = (df_ge2pho.Define("subleadPhotonBHScore", "selPho_beamHaloCNNScore[1]")
            .Define("subleadPhotonPBScore", "selPho_physBkgCNNScore[1]")
            .Define("subleadPhotonTime", "selPhoWTimeSig[1]")
            .Define("subleadPhotonTimeSig", "selPhoWTime[1]")
            .Define("subleadPhotonNonIsoScore", "selPho_nonIsoANNScore[1]")
            .Define("subleadPhotonEEnonIsoScore", "EEnonIsoScore[1]"))
        return {
            "1pho": df_1pho,
            "ge2pho": df_ge2pho,
        }
        ''' 
    
    def define_lead_photon_vars(self, df):
        return (
            df
            .Define("leadPhotonBHScore", "selPho_beamHaloCNNScore[0]")
            .Define("leadPhotonPBScore", "selPho_physBkgCNNScore[0]")
            .Define("leadPhotonTime", "selPhoWTimeSig[0]")
            .Define("leadPhotonTimeSig", "selPhoWTime[0]")
            .Define("leadPhotonNonIsoScore", "selPho_nonIsoANNScore[0]")
            .Define("leadPhotonEEnonIsoScore", "EEnonIsoScore[0]")
        )
   
    def define_regions(self, df, ch_name, mc):
        pho_late = f"(leadPhotonTimeSig > {self._threshs['late_timesig']})"
        pho_early = f"(leadPhotonTimeSig < {self._threshs['early_timesig']})"
        pho_prompt = (
            f"(leadPhotonTimeSig < {self._threshs['prompt_timesig']} && "
            f"leadPhotonTimeSig > -{self._threshs['prompt_timesig']})"
        )
    
        regions = {
            "BHCR": df.Filter(
                f"leadPhotonBHScore > {self._threshs['bh']}",
                f"leadPhoBHTag_{ch_name}"
            ),
            "earlyBHCR": df.Filter(
                f"leadPhotonBHScore > {self._threshs['bh']} && {pho_early}",
                f"leadPhoBHTag_early_{ch_name}"
            ),
            "lateBHCR": df.Filter(
                f"leadPhotonBHScore > {self._threshs['bh']} && {pho_late}",
                f"leadPhoBHTag_late_{ch_name}"
            ),
        }
    
        # EE non-iso CRs
        if "ge2pho" not in ch_name:
            eeNonIso = (
                #f"(selPho_nonIsoANNScore[0] > {self._threshs['EE_nonIso']}) && "
                f"(selPho_nonIsoANNScore[0] > {self._threshs['EE_nonIso']} && selPho_nonIsoANNScore[0] < {self._threshs['EE_veryNonIso']}) && "
                f"(selPhoEta[0] < -1.479 || selPhoEta[0] > 1.479) && "
                f"{pho_prompt}"
            )
            eeVeryNonIso = (
                f"(selPho_nonIsoANNScore[0] > {self._threshs['EE_veryNonIso']}) && "
                f"(selPhoEta[0] < -1.479 || selPhoEta[0] > 1.479) && "
                f"{pho_prompt}"
            )
        else:
            sub_prompt = (
                f"(selPhoWTimeSig[1] > -{self._threshs['prompt_timesig']} && "
                f"selPhoWTimeSig[1] < {self._threshs['prompt_timesig']})"
            )
    
            lead_nonIso = (
                f"((selPho_nonIsoANNScore[0] > {self._threshs['EE_nonIso']}) && "
                f"(selPhoEta[0] < -1.479 || selPhoEta[0] > 1.479) && {pho_prompt})"
            )
    
            sub_nonIso = (
                f"((selPho_nonIsoANNScore[1] > {self._threshs['EE_nonIso']}) && "
                f"(selPhoEta[1] < -1.479 || selPhoEta[1] > 1.479) && {sub_prompt})"
            )
            sub_notvnonIso = ( 
                f"((selPho_nonIsoANNScore[1] < {self._threshs['EE_veryNonIso']}) && "
                f"(selPhoEta[1] < -1.479 || selPhoEta[1] > 1.479) && {sub_prompt})"
            )
    
            lead_vnonIso = lead_nonIso.replace(self._threshs['EE_nonIso'], self._threshs['EE_veryNonIso'])
            sub_vnonIso = sub_nonIso.replace(self._threshs['EE_nonIso'], self._threshs['EE_veryNonIso'])
            #if wanting to include sublead photon, make sure the definitions are orthogonal bw nonIso and VeryNonIso
            #eeNonIso = f"{lead_nonIso} && {sub_notvnonIso}"
            #eeVeryNonIso = f"{lead_vnonIso} || {sub_vnonIso}"
            eeNonIso = f"{lead_nonIso}"# || {sub_nonIso}"
            eeVeryNonIso = f"{lead_vnonIso}"# || {sub_vnonIso}"
        #print("channel:",ch_name,"\n  eeNonIso",eeNonIso,"\n  eeVeryNonIso",eeVeryNonIso)
        regions["nonIsoEECR"] = df.Filter(
            eeNonIso, f"leadPhoEEnonIsoTag_{ch_name}"
        )
        regions["verynonIsoEECR"] = df.Filter(
            eeVeryNonIso, f"leadPhoEEVeryNonIsoTag_{ch_name}"
        )
    
        #if MC, define PB early/late regions and iso SRs
        if mc:
            pbEarly = f"leadPhotonPBScore > {self._threshs['pb']} && {pho_early}"
            pbLate = f"leadPhotonPBScore > {self._threshs['pb']} && {pho_late}"
            regions["earlyPBCR"] = df.Filter(pbEarly, f"leadPhoPBTag_early_{ch_name}")
            regions["latePBSR"] = df.Filter(pbLate, f"leadPhoPBTag_late_{ch_name}")
    
            lead_eeIso = (
                f"((selPho_isoANNScore[0] > {self._threshs['EE_iso']}) && "
                f"(selPhoEta[0] < -1.479 || selPhoEta[0] > 1.479) && {pho_prompt})"
            )
            regions["isoEESR"] = df.Filter(lead_eeIso,f"leadPhoEEIso_prompt_{ch_name}")
    
        return regions
    
    
    def fill_region_hists(self, df, proc_name, reg_name, ch_name, h1d, h2d):
        reg_name_key = next(k for k in self._msrs_bins if k in reg_name)
        msbins = self._msrs_bins[reg_name_key]['ms']
        rsbins = self._msrs_bins[reg_name_key]['rs']
        
        n_msbins = len(msbins) - 1
        n_rsbins = len(rsbins) - 1
   
        h2d.append(
            df.Histo2D(
	            (f"yields_{proc_name}_{ch_name}_{reg_name}", ";Ms;Rs",
                 n_msbins, msbins, n_rsbins, rsbins),
                "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
            )
        )
      
        h2d.append(
            df.Histo2D(
	            (f"MsRs_{proc_name}_{ch_name}_{reg_name}", ";Ms;Rs",
                 50, 0, 5000, 50, 0, 1.01),
                "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
            )
        )
   
        h1d.append(
            df.Histo1D(
	            (f"Rs_{proc_name}_{ch_name}_{reg_name}", "", 50, 0, 1.01),
                "rjr_Rs0", "evtFillWgt"
            )
        )
   	
        h1d.append(
            df.Histo1D(
	            (f"Ms_{proc_name}_{ch_name}_{reg_name}", "", 50, 0, 5000),
                "rjr_Ms0", "evtFillWgt"
            )
        )
   
        h1d.append(
            df.Histo1D(
	            (f"timeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -20, 20),
                "selPhoWTimeSig", "evtFillWgt"
            )
        )






    def doSignalEfficiencies(self, args):
        mGl = args.mGl
        mN2 = args.mN2
        mN1 = args.mN1
        ctau = args.ctau
      
        fileprocessor = FileProcessor()
        for proc in args.proc:
            selected_data = {}
            procstr = self.GetProcessName(proc, mGl, mN2, mN1, ctau)
            if "SMS" not in proc:
                print("Can only do efficiencies for signal, doing efficiencies for",proc)
                exit()
            mc = True
            if "PD" in proc:
                mc = False
            files = fileprocessor.GetFiles(proc, mGl, mN2, mN1, ctau)
            for file in files:
                print("file",file)
                infilename = file
                infilename = infilename[infilename.rfind("/")+1:infilename.rfind("_")]
                infilename = infilename.replace("_","\_")
                selected_data[infilename] = []

                df00 = RDataFrame("kuSkimTree", file, self._branches)
                df0 = df00.Filter("rjr_Rs.size() > 0 && rjr_Ms.size() > 0") #in case these have size 0, can lead to undefined behavior
                df = (
                    df0.Define("rjr_Rs0", "rjr_Rs[0]")
                      .Define("rjr_Ms0", "rjr_Ms[0]")
                )
    
                df1 = self.define_photon_quantities(df, self._threshs)
                #do individual presel cuts here so they are printed out
                df_metcut = df1.Filter(self._metcut,self._metcut)
                df_pts = df1.Filter(self._ptscut, self._ptscut)
                df_triggers = df1.Filter(self._triggers,"triggers")
                df_filters = df1.Filter(self._met_filters,"met_filters")
    
                #do all presel cuts
                df_presel = self.apply_preselection(df1)
    
                channels = self.define_channels(df_presel)
   
                df_list = [df_presel]
                region_list = []
                for ch_name, df_ch in channels.items():
                    print(" doing channel",ch_name)
                    df_ch = self.define_lead_photon_vars(df_ch) #done in define channels
                    regions = self.define_regions(df_ch, ch_name, mc)
                    region_list.append(regions)

                report = df00.Report()
                # select cuts
                lines = report2str(report)
                for line in lines:
                    parsed_eff = parse_eff_line(line)
                    if parsed_eff is None:
                        continue
                    selected_data[infilename].append(parsed_eff) 
            # write LaTeX table
            outfile = f"{procstr}_eff_table.tex"
            write_latex_table(outfile, selected_data)
            print("Wrote efficiencies for process",proc,"to",outfile) 

    def GetProcessName(self, proc, mGl = None, mN2 = None, mN1 = None, ctau = None):
        procstr = proc
        if proc == "METPD":
            procstr = "METFullRunII"   
        if "SMS" in proc: #assume only one mass point given
            if mGl is not None:
                procstr += "mGl"+mGl
            else:
                procstr += "mGlAll"
            if mN2 is not None:
                procstr += "mN2"+mN2
            else:
                procstr += "mN2All"
            if mN1 is not None:
                procstr += "mN1"+mN1
            else:
                procstr += "mN1All"
            if ctau is not None:
                procstr += "ct"+ctau
            else:
                procstr += "ctAll"
        return procstr

# -------------------------
# Main analysis
# -------------------------

    def runRJRAnalysis(self, args, ofilename_extra: str = ""):
        #self.define_region_selection_root_functions() 
        procs = args.proc
        mGl = args.mGl
        mN2 = args.mN2
        mN1 = args.mN1
        ctau = args.ctau
      
        fileprocessor = FileProcessor()
    
        checkForBadWgts = False 
       
        ofilename = ""
        hists1d, hists2d = [], []
        for proc in procs:
            procstr = self.GetProcessName(proc, mGl, mN2, mN1, ctau)
            files = []
            if "SMS" in proc: #assume only one mass point given
                files += fileprocessor.GetFiles(proc, mGl, mN2, mN1, ctau)
                if "gogoGZ" in files[0] and "skims_v46" in files[0]: #only do for skims_v46 gogoGZ signal samples
                    checkForBadWgts = True
            else: 
                files += fileprocessor.GetFiles(proc)
            if len(files) == 0:
                print(f"No files found for proc {proc} with ctau={args.ctau}, mGl={args.mGl}, mN2={args.mN2}, mN1={args.mN1}. Going to next process")
                continue
            #remove underscores from procstr
            procstr = procstr.replace("_","")
            if ofilename == "":
                ofilename = procstr
            else:
                ofilename += "_"+procstr
            print("Doing process",procstr)
            
            mc = True
            if "MET" in proc:
                mc = False
         
            chain = TChain("kuSkimTree")
            for f in files:
                chain.Add(f)
    
    
            df00 = RDataFrame(chain, self._branches)
            print("There are",df00.Filter("rjr_Ms.size() == 0 && nSelPhotons > 1").Count().GetValue(),"events without RJR info")
            df = df00.Filter("rjr_Rs.size() > 0 && rjr_Ms.size() > 0") #in case these have size 0, can lead to undefined behavior
            if checkForBadWgts:
                print("checking for bad event weights - will cost one (1) event loop")
                # Count rows where the column is NaN or Inf
                n_bad = (
                    df
                    .Filter("std::isnan(evtFillWgt) || std::isinf(evtFillWgt)")
                    .Count()
                    .GetValue()
                )
                has_bad = n_bad > 0
                if has_bad:
                    print("Bad weights - setting all weights to 1")
                    df = df.Redefine("evtFillWgt", "1")

            df0 = (
                df.Define("rjr_Rs0", "rjr_Rs[0]")
                  .Define("rjr_Ms0", "rjr_Ms[0]")
            )
            df1 = self.define_photon_quantities(df0, self._threshs)
            #do individual presel cuts here so they are printed out
            df_metcut = df1.Filter(self._metcut,self._metcut)
            df_pts = df1.Filter(self._ptscut, self._ptscut)
            df_triggers = df1.Filter(self._triggers,"triggers")
            df_filters = df1.Filter(self._met_filters,"met_filters")
    
            #do all presel cuts
            df_presel = self.apply_preselection(df1)
    
            #print_photon_efficiencies(df_presel)
    
            channels = self.define_channels(df_presel)
   
            df_list = [df_presel]
            for ch_name, df_ch in channels.items():
                #print(" doing channel",ch_name)
                df_ch = self.define_lead_photon_vars(df_ch) #done in define channels
    
                regions = self.define_regions(df_ch, ch_name, mc)
                hists1d.append(
                    df_ch.Histo1D(
                        (f"leadPhotonBHdiscrScore_{procstr}_{ch_name}", "", 50, 0, 1.01),
                        "leadPhotonBHScore"
                    )
                )
                
                hists1d.append(
                    df_ch.Histo1D(
                        (f"leadPhotonPBdiscrScore_{procstr}_{ch_name}", "", 50, 0, 1.01),
                        "leadPhotonPBScore"
                    )
                )
    
                hists1d.append(
                    df_ch.Histo1D(
                        (f"leadPhotonTimeSig_{procstr}_{ch_name}", "", 50, -5, 5),
                        "leadPhotonTimeSig"
                    )
                )
                for reg_name, df_reg in regions.items():
                    #print("  doing region", reg_name)
                    self.fill_region_hists(
                        df_reg, proc, reg_name, ch_name,
                        hists1d, hists2d
                    )
                    df_list.append(df_reg)
            df.Report().Print() #triggers event loop with Print() call dereferencing 
            exit() 
        #process loop end
    
    
        if ofilename_extra:
            ofilename += f"_{ofilename_extra}"
        ofilename += "_rjrObs.root"
        
        fout = TFile(ofilename, "RECREATE")
        for h in hists1d:
            h.Write()
        for h in hists2d:
            h.Write()
        fout.Close()
        print("Writing output to", ofilename)

if __name__ == "__main__":
    #import kerebos credentials to conda env if not already there
    kerb = os.getenv("KRB5CCNAME")
    if(kerb is None):
        print("Setting kerebos credentials")
        os.environ["KRB5CCNAME"] = "API:"
    parser = argparse.ArgumentParser(
        description="Run RJR analysis"
    )

    #TODO - set integrated luminosity
    #parser.add_argument(
    #    "--lumi",
    #    type=float,
    #    default=138,
    #    help="Integrated lumi for scaling MC samples"
    #)

    parser.add_argument(
        "--proc",
        "-p",
        nargs="+",
        type=str,
        required=True,
        help="Process name (e.g. MET(18), GJets, SMS-GlGl)"
    )

    parser.add_argument(
        "--mGl",
        default=None,
        help="gluino mass for SMS process"
    )
    
    parser.add_argument(
        "--mN2",
        default=None,
        help="N2 mass for SMS process"
    )
    
    parser.add_argument(
        "--mN1",
        default=None,
        help="N1 mass for SMS process"
    )
    
    parser.add_argument(
        "--ctau",
        default=None,
        help="ctau for SMS process"
    )

    # optional argument with default
    parser.add_argument(
        "--ofilename-extra",
        dest="ofilename_extra",
        type=str,
        default="",
        help="Extra string appended to output filename"
    )

   
    parser.add_argument(
        '--doPerFileEffs',
        help='show efficiencies per file',
        action='store_true',
        default=False
    )
 
    args = parser.parse_args()

    #event histograms (ie yields, Ms, Rs, time sig, etc) are weighted
    #object histograms (ie photon scores, etc) are not
    rjrana = RJRAnalysis()

    if args.doPerFileEffs:
        rjrana.doSignalEfficiencies(args)
        exit()

    rjrana.runRJRAnalysis(
        args,
        ofilename_extra=args.ofilename_extra
    )
