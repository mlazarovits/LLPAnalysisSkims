from ROOT import RDataFrame, TChain, TFile, TH1, TH2, gInterpreter, std, gSystem
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
            "nPhotons",
            "SV_nLeptonic",
            "SV_nHadronic",
            "photon_isoANNScore",
            "photon_beamHaloCNNScore",
            "photon_E",
            "photon_Eta",
            "photon_Phi",
            "photon_Pt",
            "photon_WTimeSig"
        ]
        self._threshs = { #iso threshs gotten from SMS test sample
            "EE_veryLooseIso": "0.21944034", #large8 model, 40% FPR, 95% TPR, KUBaseLine
            #"EE_veryLooseIso": "0.05", #large8 model, 40% FPR, 95% TPR, KUBaseLine
            #"EE_medIso": "0.5220398", #large8 model, 10% FPR, 90% TPR 
            "EE_medIso": "0.05", #large8 model, test 
            "EE_medIso2": "0.45", #large8 model, test 
            #"EE_tightIso" : "0.744123", #large8 model, 5% FPR, 86% TPR 
            #"EE_tightIso" : "0.5220398", #large8 model, 10% FPR, 86% TPR 
            "EE_tightIso" : "0.9", #large8 model, 10% FPR, 86% TPR 
            "EE_veryTightIso" : "0.95", #large8 model, 10% FPR, 86% TPR 
            "EB_veryLooseIso": "0.2919292",  #16-8-4 model, 20% FPR, 82% TPR, KUBaseLine
            #"EB_veryLooseIso": "0.05",  #16-8-4 model, 20% FPR, 82% TPR, KUBaseLine
            "EB_medIso": "0.05", #test 16-8-4 model 
            "EB_medIso2": "0.45", #test 16-8-4 model 
            #"EB_tightIso" : "0.7167336", #16-8-4 model, 5% FPR, ~70% TPR 
            #"EB_tightIso" : "0.5209812", #16-8-4 model, 10% FPR, ~75% TPR 
            "EB_tightIso" : "0.9", #16-8-4 model, 10% FPR, ~75% TPR 
            "EB_veryTightIso" : "0.95", #16-8-4 model, 10% FPR, ~75% TPR 
            "prompt_time": "1",
            "prompt_timesig": "2.5",
            "early_timesig": "-2.5",
            "late_timesig": "2.5",
            "bhTag": " 0.917252",
            "notbhTag" : "0.18523645" #1 - 0.81476355
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
        self._kin_bins = {}
        self._kin_bins["*"] =   {"ms" : array("d", [0, 10000]), "rs": array("d", [0, 1.0])}
        self._kin_bins["BHCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["PBCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["PBSR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._kin_bins["geq1KUBaseLinePhotonMedIso0TightIsoCR"] = {"ms": array("d", [0, 1500, 2200, 10000]),"rs" :array("d", [0, 0.15, 0.4, 1.0]) , "ptisr":array("d",[0,3000])}
        self._kin_bins["geq1KUBaseLinePhotonTightIsoSR"] = {"ms": array("d", [0, 1500, 2200, 10000]),"rs" :array("d", [0, 0.15, 0.4, 1.0]), "ptisr":array('d',[0,3000]) }
        self._kin_bins["MsCR"] = {"ms" : array("d", [0, self._kin_bins["BHCR"]['ms'][0]]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["RsCR"] = {"ms" : array("d", [0, 10000]), "rs": array("d", [self._kin_bins["BHCR"]['rs'][0], 1.0])}
        self._kin_bins["dxySigCR"] = {"ms" : array("d", [1000, 10000]), "rs": array("d", [0.15, 1.0])}

    #includes barrel photon iso presel 
    def apply_preselection(self, df):
        presel = f"{self._metcut} && {self._ptscut} && {self._triggers} && {self._met_filters}"
        return df.Filter(presel, "presel")
       
    def define_channels(self, df):
        return {
            "presel": df,
            "ge1pho": df.Filter("nBaseLinePhotons > 0 && SV_nHadronic == 0", "ge1Pho"),
            #"1pho1HadSV" : df.Filter("nSelPhotons == 1 && SV_nHadronic == 1 && SV_nLeptonic == 0","1nSelPho1HadSV"),
            #"1HadSV" : df.Filter("nSelPhotons == 0 && SV_nHadronic == 1 && SV_nLeptonic == 0","1HadSV")
        }


    def slice_vec(self, vec):
        return vec[:2]

    def MakeNewBranches(self, df, sel, phoname):
        df_ret = df
        for branch in self._branches:
            if "photon_" not in branch:
                continue
            newBranchName = branch.replace("photon",phoname)
            #print("new branch name",newBranchName)
            #define selected lead + sublead (if exists) branches
            #df_ret = df_ret.Define(newBranchName,f"{branch}[{sel}]")
            #take two leading photons that pass selection
            brsel = "ROOT::RVec<float> result; for (size_t i = 0; i < "+branch+".size(); ++i) { if ("+sel+"&& result.size() < 3) result.push_back("+branch+"[i]); if(result.size() == 2) break; } return result;"
            df_ret = df_ret.Define(newBranchName,brsel)
        df_ret = df_ret.Define(f"n{phoname}",f"{newBranchName}.size()") #take size of last made new branch
        return df_ret



    def define_regions(self, df, ch_name, mc):
        noNonPromptPhotons = "( (nBaseLinePhotons == 1 && (baseLinePhoton_WTimeSig[0] < 2.5 && baseLinePhoton_WTimeSig[0] > -2.5)) || (nBaseLinePhotons == 2 && ( baseLinePhoton_WTimeSig[0] < 2.5 && baseLinePhoton_WTimeSig[0] > -2.5) && ( baseLinePhoton_WTimeSig[1] < 2.5 && baseLinePhoton_WTimeSig[1] > -2.5)) )"
        ##SV regions
        SVgeLep1_CR = "((SV_nLeptonic >= 1) && (LeptonicSV_dxySig[0] < 200))"
        SVgeLep1_SR = "((SV_nLeptonic >= 1) && (LeptonicSV_dxySig[0] >= 200))"
        
        Ch1CRGeLep = f"{SVgeLep1_CR}"
        Ch2SRGeLep = f"{SVgeLep1_SR}"
   
        #leptonic veto in hadronic regions
        noSVLep = "(SV_nLeptonic == 0)"
        SVgeHad1_CR = "((SV_nHadronic >= 1) && (HadronicSV_dxySig[0] < 800))"
        SVgeHad1_SR = "((SV_nHadronic >= 1) && (HadronicSV_dxySig[0] > 800))" 

        Ch3CRGeHad = f"{SVgeHad1_CR} && {noNonPromptPhotons}"
        Ch4SRGeHad = f"{SVgeHad1_SR} && {noNonPromptPhotons}"

	#TODO make mixed SVhad+delayed photons regions

        ##photon preselection
        noSV = "((SV_nHadronic==0) && (SV_nLeptonic==0))"
        ge1Pho = "(nBaseLinePhotons > 0)"
        noBHPhotons = "(( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.917252) || (nBaseLinePhotons == 2 && baseLinePhoton_beamHaloCNNScore[0] < 0.917252 && baseLinePhoton_beamHaloCNNScore[1] < 0.917252) ) )"

        ##delayed photon regions
        ge1BHPhoEarly = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] < -2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] < -2.5) || (baseLinePhoton_beamHaloCNNScore[1] >= 0.917252 && baseLinePhoton_WTimeSig[1] < -2.5)) ) )"
        ge1BHPhoLate = "(nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] >= 2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] >= 2.5) || (baseLinePhoton_beamHaloCNNScore[1] >= 0.917252 && baseLinePhoton_WTimeSig[1] >= 2.5)) )"
        ge1NotBHPhoEarly = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] < -2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] < -2.5) || (baseLinePhoton_beamHaloCNNScore[1] < 0.185 && baseLinePhoton_WTimeSig[1] < -2.5)) ) )"
        ge1NotBHPhoLate = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5) || (baseLinePhoton_beamHaloCNNScore[1] < 0.185 && baseLinePhoton_WTimeSig[1] >= 2.5)) ) )"
        Ch5CRgeq1PhoBHEarly = f"{noSV} && {ge1Pho} && {ge1BHPhoEarly}"        
        Ch6CRgeq1PhoBHLate = f"{noSV} && {ge1Pho} && {ge1BHPhoLate}"        
        Ch7CRgeq1PhoNotBHEarly = f"{noSV} && {noBHPhotons} && {ge1Pho} && {ge1NotBHPhoEarly}"        
        Ch8SRgeq1PhoNotBHLate = f"{noSV} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoLate}"        
      

        ##prompt photon regions
        eq1PhoMedIsoPrompt = "(nBaseLinePhotons == 1  && ((baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] < -0.000198*baseLinePhoton_Pt[0] + 1.0188) && ( baseLinePhoton_isoANNScore[0] >= -0.000198*baseLinePhoton_Pt[0] + 0.7698)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] < 0.999) && (baseLinePhoton_isoANNScore[0] >= 0.75)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479)))"
        eq1PhoTightIsoPrompt = "(nBaseLinePhotons == 1 && (baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.000198*baseLinePhoton_Pt[0] + 1.0188)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.999)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479))"
        eq2PhoMedIsoPrompt = "(nBaseLinePhotons == 2 && ((baseLinePhoton_Pt[0] > 100 && ( baseLinePhoton_isoANNScore[0] < -0.0001*baseLinePhoton_Pt[0] + 0.96 ) && (baseLinePhoton_isoANNScore[0] >= -0.0001*baseLinePhoton_Pt[0] + 0.76 )) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] < 0.95 && baseLinePhoton_isoANNScore[0] >= 0.75)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479 )) || ((baseLinePhoton_Pt[1] > 100 && (baseLinePhoton_isoANNScore[1] < -0.0001*baseLinePhoton_Pt[1] + 0.96) && ( (baseLinePhoton_isoANNScore[1] >= -0.0001*baseLinePhoton_Pt[1] + 0.76))) || (baseLinePhoton_Pt[1] <= 100 && baseLinePhoton_isoANNScore[1] < 0.95 && baseLinePhoton_isoANNScore[1] >= 0.75) && (baseLinePhoton_Eta[1] < 1.479 && baseLinePhoton_Eta[1] > -1.479 )) )"
        eq2PhoTightIsoPrompt = "(nBaseLinePhotons == 2 && (((baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.0001*baseLinePhoton_Pt[0] + 0.96)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.95))) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479)) && ((baseLinePhoton_Pt[1] > 100 && (baseLinePhoton_isoANNScore[1] >= -0.0001*baseLinePhoton_Pt[1] + 0.96)) || (baseLinePhoton_Pt[1] <= 100 && baseLinePhoton_isoANNScore[1] >= 0.95) && (baseLinePhoton_Eta[1] < 1.479 && baseLinePhoton_Eta[1] > -1.479)))" 

        Ch9CReq1PhoMedIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq1PhoMedIsoPrompt}"
        Ch10SReq1PhoTightIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq1PhoTightIsoPrompt}"
        Ch11CReq2PhoMedIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq2PhoMedIsoPrompt}"
        Ch12SReq2PhoTightIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq2PhoTightIsoPrompt}"

        regions = {
            "Ch1CRGeLep1" : df.Filter(Ch1CRGeLep),
            "Ch3CRGeHad1" : df.Filter(Ch3CRGeHad),
            "Ch5CRgeq1PhoBHEarly" : df.Filter(Ch5CRgeq1PhoBHEarly),
            "Ch6CRgeq1PhoBHLate" : df.Filter(Ch6CRgeq1PhoBHLate),
            "Ch7CRgeq1PhoNotBHEarly" : df.Filter(Ch7CRgeq1PhoNotBHEarly),
            "Ch9CReq1PhoMedIsoPrompt" : df.Filter(Ch9CReq1PhoMedIsoPrompt),
            "Ch11CReq2PhoMedIsoPrompt" : df.Filter(Ch11CReq2PhoMedIsoPrompt),
        }
        if mc:
            regions["Ch2SRGeLep1"] = df.Filter(Ch2SRGeLep);
            regions["Ch4SRGeHad1"] = df.Filter(Ch4SRGeHad);
            regions["Ch8SRgeq1PhoNotBHLate"] = df.Filter(Ch8SRgeq1PhoNotBHLate)
            regions["Ch10SReq1PhoTightIsoPrompt"] = df.Filter(Ch10SReq1PhoTightIsoPrompt)
            regions["Ch12SReq2PhoTightIsoPrompt"] = df.Filter(Ch12SReq2PhoTightIsoPrompt)
        return regions

    def do_ms_rs_cuts(self, df, reg_name):
        #in kinematic sidebands only do cuts on other obs
        if 'Ms' in reg_name:
            return df.Filter(f"rjr_Rs0 > {self._kin_bins['BHCR']['rs'][0]}")
        elif 'Rs' in reg_name:
            return df.Filter(f"rjr_Ms0 > {self._kin_bins['BHCR']['ms'][0]}")
        else:
            if any([reg_name not in binname for binname in self._kin_bins]):
                reg_name_key = "*"
            else:
                reg_name_key = next(k for k in self._kin_bins if k in reg_name)
            return df.Filter(f"rjr_Rs0 > {self._kin_bins[reg_name_key]['rs'][0]}").Filter(f"rjr_Ms0 > {self._kin_bins[reg_name_key]['ms'][0]}")


    def fill_region_hists(self, df, proc_name, reg_name, ch_name, h1d, h2d):
        if reg_name not in self._kin_bins:
            msbins = array('d',[0, 10000])
            rsbins = array('d',[0, 1])
            ptisrbins = array('d',[0,3000])
        else:
            reg_name_key = next(k for k in self._kin_bins if k in reg_name)
            msbins = self._kin_bins[reg_name_key]['ms']
            rsbins = self._kin_bins[reg_name_key]['rs']
            if 'ptisr' not in self._kin_bins[reg_name_key].keys():
                ptisrbins = array('d',[0,3000])
            else:
                ptisrbins = self._kin_bins[reg_name_key]['ptisr']
      
        msmax = max(msbins) 
        n_msbins = len(msbins) - 1
        n_rsbins = len(rsbins) - 1
        if ch_name == reg_name or reg_name == "MsCR":
            h1d.append(
                df.Histo1D(
                    (f"photonIsoScore_{proc_name}_{ch_name}_{reg_name}", f"photonIsoScore_{proc_name}", 80, 0., 1.07),
                    "photon_isoANNScore"
                )
            )

            highptcut = "200"
            vhighptcut = "500"
            df = df.Define("isoScoreHighPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt >= {highptcut} && baseLinePhoton_Pt < {vhighptcut}]").Define("isoScoreLowPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt < {highptcut}]").Define("isoScoreVeryHighPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt >= {vhighptcut}]")
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScore_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScore_{proc_name}", 80, 0.15, 1.07),
                    "baseLinePhoton_isoANNScore"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScoreHightPt_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScore_{proc_name}", 80, 0.15, 1.07),
                    "isoScoreHighPt"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScoreVeryHightPt_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScore_{proc_name}", 80, 0.15, 1.07),
                    "isoScoreVeryHighPt"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScoreLowtPt_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScore_{proc_name}", 80, 0.15, 1.07),
                    "isoScoreLowPt"
                )
            )

            df = df.Define("baseLinePhoton_isoANNScoreGen22","baseLinePhoton_isoANNScore[baseLinePhoton_SusyId == 22]")
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScoreGen22_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScoreGen22_{proc_name}", 80, 0.15, 1.07),
                    "baseLinePhoton_isoANNScoreGen22"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonTimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -5,5),
                    "baseLinePhoton_WTimeSig"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonPt_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonPt_{proc_name}", 50, 0, 1000),
                    "baseLinePhoton_Pt"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"nbaseLinePhoton_{proc_name}_{ch_name}_{reg_name}", "", 5,0,5),
                    f"nbaseLinePhoton"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"nbaseLinePhotonTightIso_{proc_name}_{ch_name}_{reg_name}", "", 3,0,3),
                    f"nbaseLinePhotonTightIso"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"nbaseLinePhotonBeamHalo_{proc_name}_{ch_name}_{reg_name}", "", 3,0,3),
                    f"nbaseLinePhotonBeamHalo"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonTimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -5,5),
                    f"baseLinePhoton_WTimeSig"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonPt_{proc_name}_{ch_name}_{reg_name}", "", 50, 0,1000),
                    f"baseLinePhoton_Pt"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonGenSusyId_{proc_name}_{ch_name}_{reg_name}", "", 60, -30,30),
                    "baseLinePhoton_SusyId"
                )
            )
         
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
                 50, 0, 10000, 50, 0, 1),
                "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
            )
        )
        #compressed kinematics
        h2d.append(
            df.Filter("rjrIsr_PtIsr != -1").Histo2D( #remove events with ill-defined ISR trees
             (f"RISRPtISR_{proc_name}_{ch_name}_{reg_name}", ";RISR;ptISR",
                 50, ptisrbins[0], ptisrbins[-1], 50, rsbins[0],rsbins[-1]+0.1 ),
                "rjrIsr_PtIsr", "rjrIsr_RIsr", "evtFillWgt"
            )
        )   
        
        h2d.append(
            df.Histo2D(
                 (f"yieldsComp_{proc_name}_{ch_name}_{reg_name}", ";Ms;Rs",
                 n_msbins, msbins, n_rsbins, rsbins),
                "rjrIsr_RIsr", "rjrIsr_PtIsr", "evtFillWgt"
            )
        )

        #do Ms cut on Rs    
        h1d.append(
            df.Filter("rjr_Ms0 > 2000").Histo1D(
               (f"Rs_{proc_name}_{ch_name}_{reg_name}", "", 25, rsbins[0], rsbins[-1]+0.01),
                "rjr_Rs0", "evtFillWgt"
            )
        )
   
        #do Rs cut on Ms    
        h1d.append(
            df.Filter("rjr_Rs0 > 0.15").Histo1D(
               (f"Ms_{proc_name}_{ch_name}_{reg_name}", "", 25, msbins[0], msbins[-1]),
                "rjr_Ms0", "evtFillWgt"
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
                if file == "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v48_testIsoID/SMS_Sig_SVHPM100_v34_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1500_ct1_rjrskim.root":
                    print("file",file,"does not exist. Skipping.")
                    continue
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
                if(mc):
                    if "SMS" in proc:
                        lumi_factor = args.lumi 
                    else: #bkg MC for 2018 only rn
                        lumi_factor = args.lumi / 67.9
                    df1 = df1.Redefine("evtFillWgt",f"evtFillWgt*{lumi_factor}")
    
                #do individual presel cuts here so they are printed out
                df_metcut = df1.Filter(self._metcut,self._metcut)
                df_pts = df1.Filter(self._ptscut, self._ptscut)
                df_triggers = df1.Filter(self._triggers,"triggers")
                df_filters = df1.Filter(self._met_filters,"met_filters")
    
                #do all presel cuts
                df_presel = self.apply_preselection(df1)
                df_list = [df_presel]
                
                #channels = self.define_channels(df_presel)
                #ch_name = "ge1blpho"
                #df_ge1pho = df_presel.Filter("nBaseLinePhotons > 0 && SV_nLeptonic == 0 && SV_nHadronic == 0")
                ch_name = "presel"
                regions = self.define_regions(df_presel, ch_name, mc)

                #get weighted evt counts
                wt_counts = {}
                for reg in regions:
                    wt_counts[reg] = regions[reg].Sum("evtFillWgt")

                report = df00.Report()
                # select cuts
                lines = self._eff_parser.report2str(report, wt_counts)
                denom_info = self._eff_parser.get_denom_line(lines, 'presel')
                total_eff = 0
                for line in lines:
                    parsed_eff = self._eff_parser.parse_eff_line(line, denom_info)
                    if parsed_eff is None:
                        continue
                    else:
                        if(show_output):
                            print(line)
                            #print("parsed_eff",parsed_eff)
                    if parsed_eff[0] != "presel" and parsed_eff[0] != "ge1KUBaseLinePhoton" and "Early" not in parsed_eff[0] and "Late" not in parsed_eff[0] and "Prompt" not in parsed_eff[0]:
                        #print("parsed_eff",parsed_eff)
                        total_eff += parsed_eff[2]
                    selected_data[infilename].append(parsed_eff)
                # write LaTeX table
                outfile = f"{procstr}_eff_table"
                if args.ofilename_extra is not None:
                    outfile += f"_{args.ofilename_extra}"
                if outfile[-1] == "_":
                    outfile = outfile[:-1]
                outfile += ".tex" 
                self._eff_parser.write_latex_table(outfile, selected_data, args.lumi, total_eff)
                print("Wrote efficiencies for process",proc,"to",outfile)

    def GetProcessName(self, proc, mGl = None, mN2 = None, mN1 = None, ctau = None):
        procstr = proc
        if proc == "METPD":
            procstr = "METFullRunII_RunIII22_23"   
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
                if "." in ctau:
                    ctau = ctau.replace(".","p")
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
            print("There are",df00.Filter("rjr_Ms.size() == 0 && nPhotons > 1").Count().GetValue(),"events without RJR info")
            df = df00.Filter("rjr_Rs.size() > 0 && rjr_Ms.size() > 0") #in case these have size 0, can lead to undefined behavior
            df1 = (
                df.Define("rjr_Rs0", "rjr_Rs[0]")
                  .Define("rjr_Ms0", "rjr_Ms[0]")
                  .Redefine("evtFillWgt",f"evtFillWgt*{args.lumi}")
            )
            #do individual presel cuts here so they are printed out
            df_metcut = df1.Filter(self._metcut,self._metcut)
            df_pts = df1.Filter(self._ptscut, self._ptscut)
            df_triggers = df1.Filter(self._triggers,"triggers")
            df_filters = df1.Filter(self._met_filters,"met_filters")
    
            #do all presel cuts
            df_presel = self.apply_preselection(df1)
            df_list = [df_presel]
            
            #channels = self.define_channels(df_presel)
            #ch_name = "ge1blpho"
            #df_ge1pho = df_presel.Filter("nBaseLinePhotons > 0 && SV_nLeptonic == 0 && SV_nHadronic == 0")
            ch_name = "presel"
            regions = self.define_regions(df_presel, ch_name, mc)
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
    #import kerberos credentials to conda env if not already there
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

    parser.add_argument(
        "--lumi",
        default=1.,
        help="set luminosity"
    )

    # optional argument with default
    parser.add_argument(
        "--ofilename-extra","-o",
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
