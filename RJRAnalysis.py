from ROOT import RDataFrame, TChain, TFile, TH1, TH2, gInterpreter
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
            "selPhoEcalRHSumEtConeDR04",
            "selPhoHadTowOverEM",
            "selPhoTrkSumPtSolidConeDR04"
            "selPhoWTimeSig"
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
            int getTagIdx(const RVecF& scores, const RVecF& eta, const float barrel_score_thresh, const float endcap_score_thresh){
                if(scores.size() < 1){
                    return -1;
                }
                else if(scores.size() == 1){
                    float lead_eta = eta[0];
                    if(fabs(lead_eta) < 1.479){
                        if(scores[0] > barrel_score_thresh)
                            return 0;
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
                        if(scores[0] > barrel_score_thresh)
                            return 0;
                        else{
                            if(fabs(sublead_eta) < 1.479){ 
                                if(scores[1] > barrel_score_thresh)
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
                                if(scores[1] > barrel_score_thresh)
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
        
    
    def apply_preselection(self, df):
        baseline = f"{self._metcut} && {self._ptscut} && {self._triggers} && {self._met_filters}"
        return df.Filter(baseline, "baseline")
    
    
    def define_channels(self, df):
        return {
            #"1pho": df.Filter("nSelPhotons == 1", "1nSelPho"),
            "1pho": df.Filter("nSelPhotons == 1 && SV_nLeptonic == 0 && SV_nHadronic == 0", "1nSelPho"),
            #"ge2pho": df.Filter("nSelPhotons >= 2", "ge2nSelPho"),
            "ge2pho": df.Filter("nSelPhotons >= 2 && SV_nLeptonic == 0 && SV_nHadronic == 0", "ge2nSelPho"),
            #"1pho1HadSV" : df.Filter("nSelPhotons == 1 && SV_nHadronic == 1 && SV_nLeptonic == 0","1nSelPho1HadSV"),
            #"1HadSV" : df.Filter("nSelPhotons == 0 && SV_nHadronic == 1 && SV_nLeptonic == 0","1HadSV")
        }
   
    #TODO - make sure that barrel photons have isolation preselection applied as part of the CR/SR assignment 
    def define_regions(self, df, ch_name, mc):
        #test custom filters
        df_regidxs = (df
                        .Define("selPhoIdx_BHTag",f"getTagIdx(selPho_beamHaloCNNScore,{self._threshs['bh']})")
                        .Define("selPhoIdx_PBTag",f"getTagIdx(selPho_physBkgCNNScore,{self._threshs['pb']})")
                        .Define("selPhoIdx_EEnonIsoTag",f"getTagIdx(selPho_nonIsoANNScore,{self._threshs['EE_nonIso']})")
                        .Define("selPhoIdx_EEIsoTag",f"getTagIdx(selPho_isoANNScore,{self._threshs['EE_iso']})")
                        .Define("selPhoIdx_EBnonIsoTag",f"getTagIdx(selPho_nonIsoANNScore,{self._threshs['EB_nonIso']})")
                        .Define("selPhoIdx_EBIsoTag",f"getTagIdx(selPho_isoANNScore,{self._threshs['EB_iso']})")
                        .Define("selPhoIdx_nonIsoTag",f"getTagIdx(selPho_nonIsoANNScore,selPhoEta,{self._threshs['EB_nonIso']}, {self._threshs['EE_nonIso']})") #eta-inclusive region
                    )

        #df_regidxs.Filter("selPhoIdx_EEIsoTag == 1").Display(["selPhoIdx_EEIsoTag","selPho_isoANNScore","selPhoPt","selPhoWTimeSig"],45).Print()

        #nonisotag cuts
        eecut_noniso = "(selPhoEta[selPhoIdx_EEnonIsoTag] > 1.479 || selPhoEta[selPhoIdx_EEnonIsoTag] < -1.479)"
        ebcut_noniso = "(selPhoEta[selPhoIdx_EBnonIsoTag] < 1.479 && selPhoEta[selPhoIdx_EBnonIsoTag] > -1.479)"
        ee_promptcut_noniso = f"(selPhoWTimeSig[selPhoIdx_EEnonIsoTag] > -{self._threshs['prompt_timesig']} && selPhoWTimeSig[selPhoIdx_EEnonIsoTag] < {self._threshs['prompt_timesig']})" 
        eb_promptcut_noniso = f"(selPhoWTimeSig[selPhoIdx_EBnonIsoTag] > -{self._threshs['prompt_timesig']} && selPhoWTimeSig[selPhoIdx_EBnonIsoTag] < {self._threshs['prompt_timesig']})"
        #barrel only
        isopresel_cut_ebnoniso = f"selPhoEcalRHSumEtConeDR04[selPhoIdx_EBnonIsoTag] < 10. && selPhoHadTowOverEM[selPhoIdx_EBnonIsoTag] < 0.02 && selPhoTrkSumPtSolidConeDR04[selPhoIdx_EBnonIsoTag] < 6."

        #isotag cuts
        eecut_iso = "(selPhoEta[selPhoIdx_EEIsoTag] > 1.479 || selPhoEta[selPhoIdx_EEIsoTag] < -1.479)"
        ebcut_iso = "(selPhoEta[selPhoIdx_EBIsoTag] < 1.479 && selPhoEta[selPhoIdx_EBIsoTag] > -1.479)"
        ee_promptcut_iso = f"(selPhoWTimeSig[selPhoIdx_EEIsoTag] > -{self._threshs['prompt_timesig']} && selPhoWTimeSig[selPhoIdx_EEIsoTag] < {self._threshs['prompt_timesig']})" 
        eb_promptcut_iso = f"(selPhoWTimeSig[selPhoIdx_EBIsoTag] > -{self._threshs['prompt_timesig']} && selPhoWTimeSig[selPhoIdx_EBIsoTag] < {self._threshs['prompt_timesig']})"
        #barrel only
        isopresel_cut_ebiso = f"selPhoEcalRHSumEtConeDR04[selPhoIdx_EBIsoTag] < 10. && selPhoHadTowOverEM[selPhoIdx_EBIsoTag] < 0.02 && selPhoTrkSumPtSolidConeDR04[selPhoIdx_EBIsoTag] < 6."

 
        ee_loosenonIsoCut = f"selPhoIdx_EEnonIsoTag != -1 && "+ee_promptcut_noniso+f" && selPho_nonIsoANNScore[selPhoIdx_EEnonIsoTag] >= {self._threshs['EE_nonIso']} && "+eecut_noniso 
        ee_looseNotTightnonIsoCut = f"selPhoIdx_EEnonIsoTag != -1 && "+ee_promptcut_noniso+f" && selPho_nonIsoANNScore[selPhoIdx_EEnonIsoTag] >= {self._threshs['EE_nonIso']} && selPho_nonIsoANNScore[selPhoIdx_EEnonIsoTag] < {self._threshs['EE_veryNonIso']} && "+eecut_noniso 
        ee_tightnonIsoCut = f"selPhoIdx_EEnonIsoTag != -1 && "+ee_promptcut_noniso+f" && selPho_nonIsoANNScore[selPhoIdx_EEnonIsoTag] >= {self._threshs['EE_veryNonIso']} && "+eecut_noniso 
        print(f"loose EE_nonIso cut: {ee_loosenonIsoCut}")
        #print(f"loose!tight EE_nonIso cut: {ee_looseNotTightnonIsoCut}")
        #print(f"tight EE_NonIso cut: {ee_tightnonIsoCut}")

        #eb_loosenonIsoCut = f"selPhoIdx_EBnonIsoTag != -1 && "+eb_promptcut+f" && selPho_nonIsoANNScore[selPhoIdx_EBnonIsoTag] >= {self._threshs['EB_nonIso']} && "+ebcut+" && selPhoIdx_EBIsoTag == -1" 
        eb_loosenonIsoCut = f"((selPhoIdx_EBnonIsoTag == 1 && selPhoIdx_EBIsoTag == -1) || (selPhoIdx_EBnonIsoTag == 0)) && "+eb_promptcut_noniso+f" && selPho_nonIsoANNScore[selPhoIdx_EBnonIsoTag] >= {self._threshs['EB_nonIso']} && "+ebcut_noniso+" && "+isopresel_cut_ebnoniso 
        print(f"loose EB_nonIso cut: {eb_loosenonIsoCut}")
        regions = {
            "earlyBHCR": df_regidxs.Filter(
                f"selPhoIdx_BHTag != -1 && selPhoWTimeSig[selPhoIdx_BHTag] < {self._threshs['early_timesig']}",
                f"earlyBHCR_{ch_name}"
            ),
            "lateBHCR": df_regidxs.Filter(
                f"selPhoIdx_BHTag != -1 && selPhoWTimeSig[selPhoIdx_BHTag] > {self._threshs['late_timesig']}",
                f"lateBHCR_{ch_name}"
            ),
            "loosenonIsoEECR": df_regidxs.Filter( #score passes EE_nonIso thresh
                ee_loosenonIsoCut,
                f"loosenonIsoEECR_{ch_name}"
            ),
            "looseNotTightnonIsoEECR": df_regidxs.Filter( #score passes EE_nonIso thresh, but less than EE_veryNonIso thresh
                ee_loosenonIsoCut,
                f"looseNotTightIsoEECR_{ch_name}"
            ),
            "tightnonIsoEECR": df_regidxs.Filter( #score needs to pass EE_veryNonIso thresh
                ee_tightnonIsoCut,
                f"tightnonIsoEECR_{ch_name}"
            ),
            "loosenonIsoEBCR": df_regidxs.Filter( #score passes EB_nonIso thresh
                eb_loosenonIsoCut,
                f"loosenonIsoEBCR_{ch_name}"
            ),
        }
        #define kinematic sideband regions
        regions[f"MsCR"] = df_regidxs.Filter(f"rjr_Ms0 < {self._msrs_bins['BHCR']['ms'][0]}")
        regions[f"RsCR"] = df_regidxs.Filter(f"rjr_Rs0 < {self._msrs_bins['BHCR']['rs'][0]}")
        regions[f"dxySigCR"] = df_regidxs.Filter("HadronicSV_dxySig[0] < 1000")       

        df_regidxs.Display(["selPhoIdx_EBnonIsoTag","selPhoIdx_EBIsoTag","selPho_isoANNScore","selPho_nonIsoANNScore"]).Print()
        #n0idx = df_regidxs.Filter(f"(selPhoIdx_EBnonIsoTag == 0) && "+eb_promptcut+f" && selPho_nonIsoANNScore[selPhoIdx_EBnonIsoTag] >= {self._threshs['EB_nonIso']} && "+ebcut).Count()
        #n1idx = df_regidxs.Filter(f"(selPhoIdx_EBnonIsoTag == 1 && selPhoIdx_EBIsoTag == -1) && "+eb_promptcut+f" && selPho_nonIsoANNScore[selPhoIdx_EBnonIsoTag] >= {self._threshs['EB_nonIso']} && "+ebcut).Count()
        #print("n0idx",n0idx.GetValue(),"n1idx",n1idx.GetValue())
 
        #if MC, de#fine PB early/late regions and iso SRs
        if mc:
            pbEarly = f"selPhoIdx_PBTag != -1 && selPhoWTimeSig[selPhoIdx_PBTag] < {self._threshs['early_timesig']} && selPhoIdx_BHTag == -1"
            pbLate = f"selPhoIdx_PBTag != -1 && selPhoWTimeSig[selPhoIdx_PBTag] > {self._threshs['late_timesig']} && selPhoIdx_BHTag == -1"
            regions["earlyPBCR"] = df_regidxs.Filter(pbEarly, f"earlyPBCR_{ch_name}")
            regions["latePBSR"] = df_regidxs.Filter(pbLate, f"latePBSR_{ch_name}")

            eeIso = f"selPhoIdx_EEIsoTag != -1 && "+ee_promptcut_iso+"  && "+eecut_iso
            regions["isoEESR"] = df_regidxs.Filter(eeIso,f"isoEESR_{ch_name}")
            ebIso = f"selPhoIdx_EBIsoTag != -1 && "+eb_promptcut_iso+" && "+ebcut_iso+" && selPhoIdx_EBnonIsoTag == -1 && "+isopresel_cut_ebiso
            print("ebIsoCut",ebIso)
            regions["isoEBSR"] = df_regidxs.Filter(ebIso,f"isoEBSR_{ch_name}")
            #do inclusive object-multiplicity-defined region
            regions[ch_name] = df_regidxs
        
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

                df00 = RDataFrame("kuSkimTree", file, self._branches)
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
    
            #print_photon_efficiencies(df_presel)
    
            channels = self.define_channels(df_presel)
   
            df_list = [df_presel]
            for ch_name, df_ch in channels.items():
                print(" doing channel",ch_name)
                regions = self.define_regions(df_ch, ch_name, mc)
                for reg_name, df_reg in regions.items():
                    #print("  doing region", reg_name)
                    #define columns with lowest Rs/Ms cuts
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
