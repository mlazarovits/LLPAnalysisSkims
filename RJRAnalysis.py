from ROOT import RDataFrame, TChain, TFile, TH1, TH2, gInterpreter, std
from tools import EfficiencyParser
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


class RJRAnalysis:
    def __init__(self):
        self._eff_parser = EfficiencyParser()
        self._branches = [
            "rjr_Rs","rjr_Ms","selCMet","rjrPTS","rjrIsr_PtIsr","rjrIsr_Rs",
            "HadronicSV_dxySig",
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
            "SV_nLeptonic",
            "SV_nHadronic",
            "selPhoWTime",
            "selPho_beamHaloCNNScore",
            "selPho_physBkgCNNScore",
            "selPhoWTimeSig",
            "selPho_nonIsoANNScore",
            "selPho_isoANNScore",
            "selPhoEta",
            "selPhoPt",
            "selPhoEcalRHSumEtConeDR04",
            "selPhoHadTowOverEM",
            "selPhoTrkSumPtSolidConeDR04",
        ]
        self._threshs = {
            "bh": "0.917252",
            "pb": "0.81476355",
            "EE_nonIso": "0.9290591",
            "EB_nonIso": "0.99661630", #
            #"EE_veryNonIso": "0.9939665",
            "EE_veryNonIso": "0.99",
            "EE_iso" : "0.9994431", #80% efficiency, 5% bkg contanimination from SMS-GlGl ROC 
            "EB_iso" : "0.003383696", #1-nonIso threshold 
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
        #10000 is upperlimit (ie inclusive)
        self._msrs_bins = {}
        self._msrs_bins["1pho"] =   {"ms" : array("d", [0, 10000]), "rs": array("d", [0, 1.0])}
        self._msrs_bins["1HadSV"] =   {"ms" : array("d", [0, 10000]), "rs": array("d", [0., 1.0])}
        self._msrs_bins["ge2pho"] = {"ms" : array("d", [0, 10000]), "rs": array("d", [0., 1.0])}
        self._msrs_bins["BHCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._msrs_bins["PBCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._msrs_bins["PBSR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._msrs_bins["nonIsoEECR"] = {"ms": array("d", [700, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["isoEESR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["nonIsoEBCR"] = {"ms": array("d", [700, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["isoEBSR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["MsCR"] = {"ms" : array("d", [0, self._msrs_bins["BHCR"]['ms'][0]]), "rs": array("d", [0.15, 1.0])}
        self._msrs_bins["RsCR"] = {"ms" : array("d", [0, 10000]), "rs": array("d", [self._msrs_bins["BHCR"]['rs'][0], 1.0])}
        self._msrs_bins["dxySigCR"] = {"ms" : array("d", [1000, 10000]), "rs": array("d", [0.15, 1.0])}

        #declare functions for defining photon CRs
        #only works for 1 + >=2 photon channels rn
        #returns index of tagged photon
        #if ==1 photon: tag is based on lead selected photon
        #if >1 photon: tag is based on lead photon of (lead, sublead) selected photons
        gInterpreter.Declare(
            """
            using ROOT::RVecF;
            int getTagIdx(const RVecF& scores, const float score_thresh){
                if(scores.size() < 1){
                    return -1;
                }
                else if(scores.size() == 1){
                    if(scores[0] > score_thresh)
                        return 0;
                    else
                        return -1; 
                }
                else{
                    if(scores[0] > score_thresh)
                        return 0;
                    else{
                        if(scores[1] > score_thresh)
                            return 1;
                        else
                            return -1;

                    }
                }
            }
        """
        )
        gInterpreter.Declare(
            """
            using ROOT::RVecF;
            using ROOT::VecOps::RVec;
            //iso + noniso = 1; bh + pb = 1
            int getRegionIdx(const RVecF& timesigs, const RVecF& etas, const RVecF& bh_scores, const RVecF& noniso_scores, const float prompt_timesig, const float early_timesig, const float late_timesig, const float nonisoEE_scorethresh, const float isoEE_scorethresh, const float nonisoEB_scorethresh, const float isoEB_scorethresh, const float bh_scorethresh, const float pb_scorethresh){
                bool lead_prompt = (timesigs[0] < prompt_timesig) && (timesigs[0] > -prompt_timesig);
                bool lead_endcap = (etas[0] > 1.479) || (etas[0] < -1.479);
                bool lead_early = timesigs[0] < early_timesig;
                bool lead_late = timesigs[0] > late_timesig;
                float lead_nonIsoScore = noniso_scores[0];
                float lead_isoScore = (1 - noniso_scores[0]);
                float lead_BHscore = bh_scores[0];
                float lead_PBscore = (1 - bh_scores[0]);

                //loosest isolation cut - photons not passing this do not make it into the analysis
                float EE_veryNonIso = 0.9939665;


                int lead_fail = 0;
                if(lead_prompt){
                    if(lead_endcap){
                        if(lead_nonIsoScore >= EE_veryNonIso){
                            lead_fail = -11;
                        }
                        else if(lead_nonIsoScore >= nonisoEE_scorethresh && lead_nonIsoScore < EE_veryNonIso){
                            return 1;
                        }
                        else if(lead_isoScore >= isoEE_scorethresh){
                            return 2;
                        }
                        else{
                            lead_fail = -1;
                        }
                    }
                    else{ //(barrel)
                        if(lead_nonIsoScore >= nonisoEB_scorethresh){
                            return 3;
                        }
                        if(lead_isoScore >= isoEB_scorethresh){
                            return 4;
                        }
                        else{
                            lead_fail = -2;
                        }
                    }
                }
                else{ //(lead nonprompt)
                    if(lead_BHscore >= bh_scorethresh){
                        if(lead_early){
                            return 5;
                        }
                        else if(lead_late){ //late
                            return 6;
                        }
                        else{
                            lead_fail = -3;
                        }
                    }
                    else if(lead_PBscore >= pb_scorethresh){
                        if(lead_early){
                            return 7;
                        }
                        else if(lead_late){ //late
                            return 8;
                        }
                        else{
                            lead_fail = -4;
                        }
                    }
                    else{
                        lead_fail = -5;
                    }
                }
                //do only one photon case too!
                if(timesigs.size() < 2 && lead_fail < 0)
                    return lead_fail;


                //check sublead
                bool sublead_prompt = (timesigs[1] < prompt_timesig) && (timesigs[1] > -prompt_timesig);
                bool sublead_endcap = (etas[1] > 1.479) || (etas[1] < -1.479);
                bool sublead_early = timesigs[1] < early_timesig;
                bool sublead_late = timesigs[1] > late_timesig;
                
                float sublead_nonIsoScore = noniso_scores[0];
                float sublead_isoScore = (1 - noniso_scores[0]);
                float sublead_BHscore = bh_scores[0];
                float sublead_PBscore = (1 - bh_scores[0]);
                if(sublead_prompt){
                    if(sublead_endcap){
                        if(sublead_nonIsoScore >= EE_veryNonIso){
                            return -12;
                        }
                        else if(sublead_nonIsoScore >= nonisoEE_scorethresh && sublead_nonIsoScore < EE_veryNonIso){
                            return 1;
                        }
                        else if(sublead_isoScore >= isoEE_scorethresh){
                            return 2;
                        }
                        else{
                            return -6;
                        }
                    }
                    else{ //(barrel)
                        if(sublead_nonIsoScore >= nonisoEB_scorethresh){
                            return 3;
                        }
                        else if(sublead_isoScore >= isoEB_scorethresh){
                            return 4;
                        }
                        else{
                            return -7;
                        }
                    }
                }
                else{ //(sublead nonprompt)
                    if(sublead_BHscore >= bh_scorethresh){
                        if(sublead_early){
                            return 5;
                        }
                        else if(sublead_late){
                            return 6;
                        }
                        else{
                            return -8;
                        }
                    }
                    else if(sublead_PBscore >= pb_scorethresh){
                        if(sublead_early){
                            return 7;
                        }
                        else if(sublead_late){
                            return 8;
                        }
                        else{
cout << "sublead timesig " << timesigs[1] << endl;
                            return -9;
                        }
                    }
                    else{
                        return -10;
                    }
                }
                
            } 
            """
        )
        gInterpreter.Declare(
            """
            using ROOT::RVecF;
            int getTagIdx(const RVecF& scores, const RVecF& eta, const ROOT::VecOps::RVec<bool>& iso_presel, const map<string,float> score_thresh){
                float barrel_score_thresh = score_thresh.at("EB");
                float endcap_score_thresh = score_thresh.at("EE");
                if(scores.size() < 1){
                    return -1;
                }
                else if(scores.size() == 1){
                    float lead_eta = eta[0];
                    if(fabs(lead_eta) < 1.479){
                        if(iso_presel[0]){
                            if(scores[0] > barrel_score_thresh)
                                return 0;
                            else
                                return -1;
                        }
                        else
                            return -1; 
                    }
                    else{
                        if(scores[0] > endcap_score_thresh)
                            return 0;
                        else
                            return -1; 
                    }
                }
                else{
                    float lead_eta = eta[0];
                    float sublead_eta = eta[1];
                    if(fabs(lead_eta) < 1.479){
                        if(scores[0] > barrel_score_thresh && iso_presel[0])
                                return 0;
                        else{
                            if(fabs(sublead_eta) < 1.479){ 
                                if(scores[1] > barrel_score_thresh && iso_presel[1])
                                    return 1;
                                else
                                    return -1;
                            	
                            }
                            else{
                                if(scores[1] > endcap_score_thresh)
                                    return 1;
                                else
                                    return -1;
                            }
                        }
                    }
                    else{
                        if(scores[0] > endcap_score_thresh)
                            return 0;
                        else{
                            if(fabs(sublead_eta) < 1.479){ 
                                if(scores[1] > barrel_score_thresh && iso_presel[1])
                                    return 1;
                                else
                                    return -1;
                            	
                            }
                            else{
                                if(scores[1] > endcap_score_thresh)
                                    return 1;
                                else
                                    return -1;
                            }
                        }

                    }
                }
            }
        """
        )
        
    #includes barrel photon iso presel 
    def apply_preselection(self, df):
        baseline = f"{self._metcut} && {self._ptscut} && {self._triggers} && {self._met_filters}"
        return df.Filter(baseline, "baseline")
       
    def do_barrel_iso_presel(self, df): 
        #do barrel photon filtering
        df = df.Define("passBarrelIsoPresel", "((selPhoEcalRHSumEtConeDR04 < 10.) & (selPhoHadTowOverEM < 0.02) & (selPhoTrkSumPtSolidConeDR04 < 6.)) | ((selPhoEta > 1.479) | (selPhoEta < -1.479))")
        for branch in self._branches:
            if "selPho" not in branch:
                continue
            df = df.Redefine(branch,f"{branch}[passBarrelIsoPresel]")
        df = df.Redefine("nSelPhotons","selPhoPt.size()")
        return df 
    
    def define_channels(self, df):
        return {
            "1pho": df.Filter("nSelPhotons == 1 && SV_nLeptonic == 0 && SV_nHadronic == 0", "1nSelPho"),
            "ge2pho": df.Filter("nSelPhotons >= 2 && SV_nLeptonic == 0 && SV_nHadronic == 0", "ge2nSelPho"),
            #"1pho1HadSV" : df.Filter("nSelPhotons == 1 && SV_nHadronic == 1 && SV_nLeptonic == 0","1nSelPho1HadSV"),
            #"1HadSV" : df.Filter("nSelPhotons == 0 && SV_nHadronic == 1 && SV_nLeptonic == 0","1HadSV")
        }
   
    def define_regions(self, df, ch_name, mc):

        df_regs = df.Define("regionIdx",f"getRegionIdx(selPhoWTimeSig, selPhoEta, selPho_beamHaloCNNScore, selPho_nonIsoANNScore, {self._threshs['prompt_timesig']}, {self._threshs['early_timesig']}, {self._threshs['late_timesig']}, {self._threshs['EE_nonIso']}, {self._threshs['EE_iso']}, {self._threshs['EB_nonIso']}, {self._threshs['EB_iso']}, {self._threshs['bh']}, {self._threshs['pb']})")

        regions = {
            "earlyBHCR": df_regs.Filter(
                "regionIdx == 5",
                f"earlyBHCR_{ch_name}"
            ),
            "lateBHCR": df_regs.Filter(
                "regionIdx == 6",
                f"lateBHCR_{ch_name}"
            ),
            "loosenonIsoEECR": df_regs.Filter( #score passes EE_nonIso thresh
                "regionIdx == 1",
                f"loosenonIsoEECR_{ch_name}"
            ),
            "loosenonIsoEBCR": df_regs.Filter( #score passes EE_nonIso thresh
                "regionIdx == 3",
                f"loosenonIsoEBCR_{ch_name}"
            ),
            #"looseNotTightnonIsoEECR": df_regidxs.Filter( #score passes EE_nonIso thresh, but less than EE_veryNonIso thresh
            #    ee_loosenonIsoCut,
            #    f"looseNotTightIsoEECR_{ch_name}"
            #),
            #"tightnonIsoEECR": df_regidxs.Filter( #score needs to pass EE_veryNonIso thresh
            #    ee_tightnonIsoCut,
            #    f"tightnonIsoEECR_{ch_name}"
            #),
            #"loosenonIsoEBCR": df_regidxs.Filter( #score passes EB_nonIso thresh
            #    eb_loosenonIsoCut,
            #    f"loosenonIsoEBCR_{ch_name}"
            #),
        }
        #define kinematic sideband regions
        regions[f"MsCR"] = df_regs.Filter(f"rjr_Ms0 < {self._msrs_bins['BHCR']['ms'][0]}")
        regions[f"RsCR"] = df_regs.Filter(f"rjr_Rs0 < {self._msrs_bins['BHCR']['rs'][0]}")
        regions[f"dxySigCR"] = df_regs.Filter("HadronicSV_dxySig[0] < 1000")       

        #if MC, de#fine PB early/late regions and iso SRs
        if mc:
            regions["earlyPBCR"] = df_regs.Filter("regionIdx == 7", f"earlyPBCR_{ch_name}")
            regions["latePBSR"] = df_regs.Filter("regionIdx == 8", f"latePBSR_{ch_name}")

            regions["isoEESR"] = df_regs.Filter("regionIdx == 2",f"isoEESR_{ch_name}")
            regions["isoEBSR"] = df_regs.Filter("regionIdx == 4",f"isoEBSR_{ch_name}")
            #do inclusive object-multiplicity-defined region
            regions[ch_name] = df_regs
        regions["fail"] = df_regs.Filter("regionIdx < 0",f"fail_{ch_name}")
        #10 ways to fail
        if regions["fail"].Count().GetValue() > 10:
            for i in range(1,13):
                regions[f"failMode_{i}"] = regions["fail"].Filter(f"regionIdx == -{i}",f"fail_mode{i}_{ch_name}")
        return regions
   
    def do_ms_rs_cuts(self, df, reg_name):
        #in kinematic sidebands only do cuts on other obs
        if 'Ms' in reg_name:
            return df.Filter(f"rjr_Rs0 > {self._msrs_bins['BHCR']['rs'][0]}")
        elif 'Rs' in reg_name:
            return df.Filter(f"rjr_Ms0 > {self._msrs_bins['BHCR']['ms'][0]}")
        else:
            #if reg_name not in any(self._msrs_bins):
            #    print("Region",reg_name,"does not have MsRs yields binnings defined. Please add.")
            #    exit()
            reg_name_key = next(k for k in self._msrs_bins if k in reg_name)
            return df.Filter(f"rjr_Rs0 > {self._msrs_bins[reg_name_key]['rs'][0]}").Filter(f"rjr_Ms0 > {self._msrs_bins[reg_name_key]['ms'][0]}")
 
    def fill_region_hists(self, df, proc_name, reg_name, ch_name, h1d, h2d):
        reg_name_key = next(k for k in self._msrs_bins if k in reg_name)
        msbins = self._msrs_bins[reg_name_key]['ms']
        rsbins = self._msrs_bins[reg_name_key]['rs']
       
        if "CR" in reg_name:
            msmax = 4000
        else:
            msmax = 5000 
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
                 25, 1000, msmax, 25, rsbins[0], 1.0),
                "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
            )
        )
        #compressed kinematics
        h2d.append(
            df.Filter("rjrIsr_PtIsr != -1").Histo2D( #remove events with ill-defined ISR trees
	            (f"RISRPtISR_{proc_name}_{ch_name}_{reg_name}", ";RISR;ptISR",
                 25, 0, rsbins[-1]+0.01, 25, 50, 1500),
                "rjrIsr_RIsr", "rjrIsr_PtIsr", "evtFillWgt"
            )
        )
   
        h1d.append(
            df.Histo1D(
	            (f"Rs_{proc_name}_{ch_name}_{reg_name}", "", 25, rsbins[0], rsbins[-1]+0.01),
                "rjr_Rs0", "evtFillWgt"
            )
        )
   	
        h1d.append(
            df.Histo1D(
	            (f"Ms_{proc_name}_{ch_name}_{reg_name}", "", 25, msbins[0], msbins[-1]),
                "rjr_Ms0", "evtFillWgt"
            )
        )
   
        h1d.append(
            df.Histo1D(
	            (f"timeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -20, 20),
                "selPhoWTimeSig", "evtFillWgt"
            )
        )






    def doSignalEfficiencies(self, args, show_output=False):
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

                df00 = RDataFrame("kuSkimTree", file)
                df0 = df00.Filter("rjr_Rs.size() > 0 && rjr_Ms.size() > 0") #in case these have size 0, can lead to undefined behavior
                df1 = (
                    df0.Define("rjr_Rs0", "rjr_Rs[0]")
                      .Define("rjr_Ms0", "rjr_Ms[0]")
                )
    
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
                    regions = self.define_regions(df_ch, ch_name, mc)
                    region_list.append(regions)

                report = df00.Report()
                # select cuts
                lines = self._eff_parser.report2str(report)
                denom_info = self._eff_parser.get_denom_line(lines, 'baseline')
                for line in lines:
                    parsed_eff = self._eff_parser.parse_eff_line(line, denom_info)
                    if parsed_eff is None:
                        continue
                    else:
                        if(show_output):
                            print(line)
                            print("parsed_eff",parsed_eff)
                    selected_data[infilename].append(parsed_eff) 
            # write LaTeX table
            outfile = f"{procstr}_eff_table"
            if args.ofilename_extra is not None:
                outfile += f"_{args.ofilename_extra}"
            outfile += ".tex" 
            eff_parser.write_latex_table(outfile, selected_data)
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
                if(len(files) < 1):
                    print("No files found for ",proc,"with mGl",mGl,"mN2",mN2,"mN1",mN1,"ctau",ctau)
                    continue
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
    
    
            df00 = RDataFrame(chain)
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

            df1 = (
                df.Define("rjr_Rs0", "rjr_Rs[0]")
                  .Define("rjr_Ms0", "rjr_Ms[0]")
            )
            #do individual presel cuts here so they are printed out
            df_metcut = df1.Filter(self._metcut,self._metcut)
            df_pts = df1.Filter(self._ptscut, self._ptscut)
            df_triggers = df1.Filter(self._triggers,"triggers")
            df_filters = df1.Filter(self._met_filters,"met_filters")
    
            #do all presel cuts
            df_presel = self.apply_preselection(df1)
            #do barrel iso presel
            df_isopresel = self.do_barrel_iso_presel(df_presel)
    
            channels = self.define_channels(df_isopresel)
 
            df_list = [df_presel,df_isopresel]
            for ch_name, df_ch in channels.items():
                print(" doing channel",ch_name)
                regions = self.define_regions(df_ch, ch_name, mc)
                for reg_name, df_reg in regions.items():
                    #print("  doing region", reg_name)
                    #define columns with lowest Rs/Ms cuts
                    if "fail" in reg_name:
                        df_list.append(df_reg)
                        continue
                    df_reg = self.do_ms_rs_cuts(df_reg, reg_name)
                    self.fill_region_hists(
                        df_reg, procstr, reg_name, ch_name,
                        hists1d, hists2d
                    )
                    df_list.append(df_reg)
            report = df00.Report() 
            if(args.showOutput):
                # select cuts
                lines = self._eff_parser.report2str(report)
                for line in lines:
                    print(line)
            else:
                report.Print() #triggers event loop with Print() call dereferencing
        #process loop end
    
        if procstr == "":
            return 
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
        "--ofilename-extra","-e",
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
 
    parser.add_argument(
        '--showOutput',
        help='print efficiencies per file',
        action='store_true',
        default=False
    )

    args = parser.parse_args()

    #event histograms (ie yields, Ms, Rs, time sig, etc) are weighted
    #object histograms (ie photon scores, etc) are not
    rjrana = RJRAnalysis()

    if args.doPerFileEffs:
        rjrana.doSignalEfficiencies(args,args.showOutput)
        exit()

    rjrana.runRJRAnalysis(
        args,
        ofilename_extra=args.ofilename_extra
    )
