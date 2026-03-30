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
        #self._photypes = ["VLooseIsoPhoton","MyBaseLinePhoton", "KUBaseLinePhoton"]
        self._photypes = ["KUBaseLinePhoton"]
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
            #"EE_medIso": "0.5220398", #large8 model, 10% FPR, 90% TPR 
            "EE_medIso": "0.4", #large8 model, test 
            "EE_medIso2": "0.45", #large8 model, test 
            #"EE_tightIso" : "0.744123", #large8 model, 5% FPR, 86% TPR 
            "EE_tightIso" : "0.5220398", #large8 model, 10% FPR, 86% TPR 
            "EB_veryLooseIso": "0.2919292",  #16-8-4 model, 20% FPR, 82% TPR, KUBaseLine
            "EB_medIso": "0.4", #test 16-8-4 model 
            "EB_medIso2": "0.45", #test 16-8-4 model 
            #"EB_tightIso" : "0.7167336", #16-8-4 model, 5% FPR, ~70% TPR 
            "EB_tightIso" : "0.5209812", #16-8-4 model, 10% FPR, ~75% TPR 
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
        self._msrs_bins = {}
        self._msrs_bins["*"] =   {"ms" : array("d", [0, 10000]), "rs": array("d", [0, 1.0])}
        self._msrs_bins["BHCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._msrs_bins["PBCR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._msrs_bins["PBSR"] = {"ms" : array("d", [700, 1000, 10000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._msrs_bins["nonIsoEECR"] = {"ms": array("d", [700, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["isoEESR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["nonIsoEBCR"] = {"ms": array("d", [700, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["isoEBSR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["isoSR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["promptNonIsoCR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["promptVeryNonIsoCR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["promptNotVeryNonIsoCR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
        self._msrs_bins["promptIsoSR"] = {"ms": array("d", [500, 1000, 10000]),"rs" :array("d", [0.15, 0.4, 1.0]) }
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
                float lead_isoScore = (lead_nonIsoScore == -999 ? -999 : (1 - noniso_scores[0]));
                float lead_BHscore = bh_scores[0];
                float lead_PBscore = (1 - bh_scores[0]);
                //loosest isolation cut - photons not passing this do not make it into the analysis
                float EE_veryveryNonIso = 0.9939665;
                int lead_fail = 0;
                if(lead_prompt){
                    if(lead_endcap){
                        if(lead_nonIsoScore >= EE_veryveryNonIso){
                            lead_fail = -1;
                        }
                        else if(lead_nonIsoScore >= nonisoEE_scorethresh && lead_nonIsoScore < EE_veryveryNonIso){
                            return 1;
                        }
                        else if(lead_isoScore >= isoEE_scorethresh){
                            return 2;
                        }
                        else{
                            lead_fail = -2;
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
                            lead_fail = -3;
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
cout << "lead time sig " << timesigs[0] << endl;
                            lead_fail = -4;
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
cout << "lead time sig " << timesigs[0] << endl;
                            lead_fail = -5;
                        }
                    }
                    else{
                        lead_fail = -6;
                    }
                }
                //do only one photon case too!
                if(timesigs.size() < 2 && lead_fail < 0){
                    return lead_fail;
                }

                //check sublead
                bool sublead_prompt = (timesigs[1] < prompt_timesig) && (timesigs[1] > -prompt_timesig);
                bool sublead_endcap = (etas[1] > 1.479) || (etas[1] < -1.479);
                bool sublead_early = timesigs[1] < early_timesig;
                bool sublead_late = timesigs[1] > late_timesig;
                
                float sublead_nonIsoScore = noniso_scores[1];
                float sublead_isoScore = (sublead_nonIsoScore == -999 ? -999 : (1 - noniso_scores[1]));
                float sublead_BHscore = bh_scores[1];
                float sublead_PBscore = (1 - bh_scores[1]);
                if(sublead_prompt){
                    if(sublead_endcap){
                        if(sublead_nonIsoScore >= EE_veryveryNonIso){
                            return -100;
                        }
                        else if(sublead_nonIsoScore >= nonisoEE_scorethresh && sublead_nonIsoScore < EE_veryveryNonIso){
                            return 100;
                        }
                        else if(sublead_isoScore >= isoEE_scorethresh){
                            return 200;
                        }
                        else{
                            return -200;
                        }
                    }
                    else{ //(barrel)
                        if(sublead_nonIsoScore >= nonisoEB_scorethresh){
                            return 300;
                        }
                        else if(sublead_isoScore >= isoEB_scorethresh){
                            return 400;
                        }
                        else{
                            return -300;
                        }
                    }
                }
                else{ //(sublead nonprompt)
                    if(sublead_BHscore >= bh_scorethresh){
                        if(sublead_early){
                            return 500;
                        }
                        else if(sublead_late){
                            return 600;
                        }
                        else{
cout << "sublead time sig " << timesigs[1] << endl;
                            return -400;
                        }
                    }
                    else if(sublead_PBscore >= pb_scorethresh){
                        if(sublead_early){
                            return 700;
                        }
                        else if(sublead_late){
                            return 800;
                        }
                        else{
cout << "sublead time sig " << timesigs[1] << endl;
                            return -500;
                        }
                    }
                    else{
                        return -600;
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
        presel = f"{self._metcut} && {self._ptscut} && {self._triggers} && {self._met_filters}"
        #presel = f"{self._metcut} && {self._triggers} && {self._met_filters}"
        return df.Filter(presel, "presel")
       
    def define_channels(self, df):
        return {
            "presel": df,
            "ge1pho": df.Filter("nPhotons > 0 && SV_nLeptonic == 0 && SV_nHadronic == 0", "ge1Pho"),
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
        #df_regs = df.Define("regionIdx",f"getRegionIdx(selPhoWTimeSig, selPhoEta, selPho_beamHaloCNNScore, selPho_nonIsoANNScore, {self._threshs['prompt_timesig']}, {self._threshs['early_timesig']}, {self._threshs['late_timesig']}, {self._threshs['EE_nonIso']}, {self._threshs['EE_iso']}, {self._threshs['EB_nonIso']}, {self._threshs['EB_iso']}, {self._threshs['bh']}, {self._threshs['pb']})")
        #eepho = "(photon_Eta < -1.479 || photon_Eta > 1.479)"
        #ebpho = "(photon_Eta > -1.479 && photon_Eta < 1.479)"
        #eevliso = f"({eepho} && photon_isoANNScore > {self._threshs['EE_veryLooseIso']})"
        #ebvliso = f"({ebpho} && photon_isoANNScore > {self._threshs['EB_veryLooseIso']})"
        #ptmin = "(photon_Pt > 30)"
        eepho = "(photon_Eta[i] < -1.479 || photon_Eta[i] > 1.479)"
        ebpho = "(photon_Eta[i] > -1.479 && photon_Eta[i] < 1.479)"
        eevliso = f"({eepho} && photon_isoANNScore[i] > {self._threshs['EE_veryLooseIso']})"
        ebvliso = f"({ebpho} && photon_isoANNScore[i] > {self._threshs['EB_veryLooseIso']})"
        vliso = f"({eevliso} || {ebvliso})"       
        ptmin = "(photon_Pt[i] > 30)"
        pixveto = "(photon_PixSeed[i] == 0)"
        baseline_sel = f"{vliso} && {ptmin} && {pixveto}"
 
        df_new = df.Define("nPho","photon_Pt.size()")
        df_new = self.MakeNewBranches(df_new, "photon_baseline[i] == 1", "KUBaseLinePhoton")
        df_new = self.MakeNewBranches(df_new, f"{vliso}", "VLooseIsoPhoton")
        df_new = self.MakeNewBranches(df_new, f"{baseline_sel}", "MyBaseLinePhoton")


        regions = {} 
        for photype in self._photypes:
            #define medIso photons and tightIso photons from lead + sublead kubaseline photons
            eepho = f"({photype}_Eta < -1.479 || {photype}_Eta > 1.479)"
            ebpho = f"({photype}_Eta > -1.479 && {photype}_Eta < 1.479)"
            eemediso = f"({eepho} && {photype}_isoANNScore > {self._threshs['EE_medIso']} && {photype}_isoANNScore <= {self._threshs['EE_tightIso']})"
            ebmediso = f"({ebpho} && {photype}_isoANNScore > {self._threshs['EB_medIso']} && {photype}_isoANNScore <= {self._threshs['EB_tightIso']})"
            mediso = f"({eemediso} || {ebmediso})"
            #validation regions
            eemediso1 = f"({eepho} && {photype}_isoANNScore > {self._threshs['EE_medIso']} && {photype}_isoANNScore <= {self._threshs['EE_medIso2']})"
            ebmediso1 = f"({ebpho} && {photype}_isoANNScore > {self._threshs['EB_medIso']} && {photype}_isoANNScore <= {self._threshs['EB_medIso2']})"
            mediso1 = f"({eemediso1} || {ebmediso1})"
            eemediso2 = f"({eepho} && {photype}_isoANNScore > {self._threshs['EE_medIso2']} && {photype}_isoANNScore <= {self._threshs['EE_tightIso']})"
            ebmediso2 = f"({ebpho} && {photype}_isoANNScore > {self._threshs['EB_medIso2']} && {photype}_isoANNScore <= {self._threshs['EB_tightIso']})"
            mediso2 = f"({eemediso2} || {ebmediso2})"
            eetightiso = f"({eepho} && {photype}_isoANNScore > {self._threshs['EE_tightIso']})"
            ebtightiso = f"({ebpho} && {photype}_isoANNScore > {self._threshs['EB_tightIso']})"
            tightiso = f"({eetightiso} || {ebtightiso})"
            ##update to == 0 bh && == 0 nonprompt
            bh = f"{photype}_beamHaloCNNScore > {self._threshs['bhTag']}" 
            pb = f"{photype}_beamHaloCNNScore <= {self._threshs['notbhTag']}"
            nonpromptsel = f"{photype}_WTimeSig <= {self._threshs['early_timesig']} || {photype}_WTimeSig > {self._threshs['late_timesig']}" 

            df_phoTag = (df_new
                        .Define(f"{photype}BeamHalo_WTimeSig",f"{photype}_WTimeSig[{bh}]").Define(f"n{photype}BeamHalo",f"{photype}BeamHalo_WTimeSig.size()")
                        .Define(f"{photype}NotBeamHalo_WTimeSig",f"{photype}_WTimeSig[{pb}]").Define(f"n{photype}NotBeamHalo",f"{photype}NotBeamHalo_WTimeSig.size()")
                        .Define(f"{photype}TightIso_WTimeSig",f"{photype}_WTimeSig[{tightiso}]").Define(f"n{photype}TightIso",f"{photype}TightIso_WTimeSig.size()")
                        .Define(f"{photype}MedIso_WTimeSig",f"{photype}_WTimeSig[{mediso}]").Define(f"n{photype}MedIso",f"{photype}MedIso_WTimeSig.size()")
                        .Define(f"{photype}MedIso1_WTimeSig",f"{photype}_WTimeSig[{mediso1}]").Define(f"n{photype}MedIso1",f"{photype}MedIso1_WTimeSig.size()")
                        .Define(f"{photype}MedIso2_WTimeSig",f"{photype}_WTimeSig[{mediso2}]").Define(f"n{photype}MedIso2",f"{photype}MedIso2_WTimeSig.size()")
                        .Define(f"nonPromptTag",f"{photype}_WTimeSig[{nonpromptsel}]").Define(f"n{photype}Nonprompt","nonPromptTag.size()")
                )
            regions[f"ge1{photype}"] = df_phoTag.Filter(f"n{photype} > 0",f"ge1{photype}")
            ####can pass ge1pho df as input (st above cut has 100% efficiency) to see efficiency of 0tagged photon events (not just from Zs)

            #df_phoTag = df_phoTag.Define("leadNotBHTimeSig",f"{photype}NotBeamHalo_WTimeSig[0]").Define("BHCR",f"n{photype}BeamHalo > 0").Define(f"ge1{photype}NotBeamHalo0BHEarlyCR",f"n{photype}Nonprompt > 0 && {photype}NotBeamHalo_WTimeSig[0] <= {self._threshs['early_timesig']}").Define(f"ge1{photype}NotBeamHalo0BHLateSR",f"n{photype}Nonprompt > 0 && {photype}NotBeamHalo_WTimeSig[0] > {self._threshs['late_timesig']}").Define(f"geq1{photype}TightIsoSR",f"n{photype}Nonprompt == 0 && n{photype}TightIso > 0 && n{photype}BeamHalo == 0")
            #df_phoTag.Display([f"n{photype}TightIso",f"geq1{photype}TightIsoSR"]).Print()

            #bh regions - all CRs
            regions[f"ge1{photype}BeamHaloCR"] = df_phoTag.Filter(f"n{photype}BeamHalo > 0",f"ge1{photype}BeamHaloCR")
            bhearly = f"{photype}BeamHalo_WTimeSig[0] <= {self._threshs['early_timesig']}" #based on lead bh-tagged photon
            regions[f"ge1{photype}BeamHaloEarlyCR"] = df_phoTag.Filter(f"n{photype}BeamHalo > 0 && {bhearly}",f"  ge1{photype}BeamHaloEarlyCR")
            bhlate = f"{photype}BeamHalo_WTimeSig[0] > {self._threshs['late_timesig']}" #based on lead bh-tagged photon
            regions[f"ge1{photype}BeamHaloLateCR"] = df_phoTag.Filter(f"n{photype}BeamHalo > 0 && {bhlate}",f"  ge1{photype}BeamHaloLateCR")
            bhprompt = f"{photype}BeamHalo_WTimeSig[0] > -{self._threshs['prompt_timesig']} && {photype}BeamHalo_WTimeSig[0] < {self._threshs['prompt_timesig']}" #based on lead bh-tagged photon
            regions[f"ge1{photype}BeamHaloPrompt"] = df_phoTag.Filter(f"n{photype}BeamHalo > 0 && {bhprompt}",f"  ge1{photype}BeamHaloPrompt")

            ###############################################################
            ########## ALL REGIONS BELOW REQUIRE == 0 BH PHOTONS ##########
            ###############################################################
            nobh = f"n{photype}BeamHalo == 0"

            ##### >= 1 nonprompt pho #####
            ge1nonprompt = f"n{photype}Nonprompt > 0"
            #looks for at least 1 nonprompt baseline photon 

            #not bh regions
            ge1notbh0bh = f"n{photype}NotBeamHalo > 0 && {nobh} && {ge1nonprompt}"
            regions[f"ge1{photype}NotBeamHalo0BH"] = df_phoTag.Filter(f"{ge1notbh0bh}",f"ge1{photype}NotBeamHalo0BH")

            pbearly = f"{photype}NotBeamHalo_WTimeSig[0] <= {self._threshs['early_timesig']}"
            ge1notbh0bhearly = f"{ge1notbh0bh} && {pbearly}"
            regions[f"ge1{photype}NotBeamHalo0BHEarlyCR"] = df_phoTag.Filter(f"{ge1notbh0bhearly}",f"  ge1{photype}NotBeamHalo0BHEarlyCR")

            pblate = f"{photype}NotBeamHalo_WTimeSig[0] > {self._threshs['late_timesig']}"
            ge1notbh0bhlate = f"{ge1notbh0bh} && {pblate}"
            regions[f"ge1{photype}NotBeamHalo0BHLateSR"] = df_phoTag.Filter(f"{ge1notbh0bhlate}",f"  ge1{photype}NotBeamHalo0BHLateSR")
            
            pbprompt = f"{photype}NotBeamHalo_WTimeSig[0] > -{self._threshs['prompt_timesig']} && {photype}NotBeamHalo_WTimeSig[0] < {self._threshs['prompt_timesig']}" #based on lead pb-tagged photon
            #would like to put this region with iso prompt regions
            #ge1notbh0bh = f"n{photype}NotBeamHalo > 0 && {nobh} && {ge1nonprompt}"
            ge1notbh0bhprompt = f"({ge1notbh0bh} && {pbprompt})"
            regions[f"ge1{photype}NotBeamHalo0BHPrompt"] = df_phoTag.Filter(f"{ge1notbh0bhprompt}",f"  ge1{photype}NotBeamHalo0BHPrompt")

            #for completeness of region definitions
            regions[f"0{photype}NotBeamHalo0BH"] = df_phoTag.Filter(f"n{photype}NotBeamHalo == 0 && {nobh} && {ge1nonprompt}",f"0{photype}NotBeamHalo0BH")


            ##### == 0 nonprompt pho #####            
            #all iso regions are orthogonal to bh regions (req ==0 bh photons)
            #all iso regions are orthogonal to pb regions (req ==0 nonprompt photons)
            #eq0nonprompt = f"n{photype}Nonprompt == 0"
            #need to include ge1{photype}NotBeamHalo0BHPrompt in iso region
            #either 0 nonprompt photons (orthogonal to notbeamhalo early/late) OR (at least 1 nonprompt photon AND at least 1 notbeamhalo photon AND lead is prompt)
            eq0nonprompt = f"(n{photype}Nonprompt == 0 || {ge1notbh0bhprompt})"
    

            tightisoSR = f"n{photype}TightIso > 0 && {nobh} && {eq0nonprompt}"
            regions[f"geq1{photype}TightIsoSR"] = df_phoTag.Filter(f"{tightisoSR}", f"geq1{photype}TightIsoSR")
            #for genid validation
            #ge1tightiso0mediso = f"n{photype}TightIso > 0 && n{photype}MedIso == 0 && {nobh} && {eq0nonprompt}"
            #regions[f"geq1{photype}TightIso0Med"] = df_phoTag.Filter(f"{ge1tightiso0mediso}", f"  geq1{photype}TightIso0MedIso")
            
            medisoCR = f"n{photype}MedIso > 0 && n{photype}TightIso == 0 && {nobh} && {eq0nonprompt}"
            regions[f"geq1{photype}MedIso0TightIsoCR"] = df_phoTag.Filter(f"{medisoCR}", f"geq1{photype}MedIso0TightIsoCR")
            #validation regions
            mediso1CR = f"n{photype}MedIso1 > 0 && n{photype}TightIso == 0 && n{photype}MedIso2 == 0 && {nobh} && {eq0nonprompt} && {nobh}"
            regions[f"geq1{photype}MedIso1VR"] = df_phoTag.Filter(f"{mediso1CR}", f"geq1{photype}MedIso1VR")
            mediso2CR = f"n{photype}MedIso2 > 0 && n{photype}TightIso == 0 && {nobh} && {eq0nonprompt} && {nobh}"
            regions[f"geq1{photype}MedIso2VR"] = df_phoTag.Filter(f"{mediso2CR}", f"geq1{photype}MedIso2VR")

            eq0mediso0tightiso = f"n{photype}MedIso == 0 && n{photype}TightIso == 0 && {nobh} && {eq0nonprompt}"
            regions[f"0{photype}MedIso0{photype}TightIso"] = df_phoTag.Filter(f"{eq0mediso0tightiso}", f"0{photype}MedIso0{photype}TightIso")

            #1/2 photon channels
            #nphos = [1, 2]
            #for npho in nphos:
            #    regions[f"eq{npho}{photype}TightIso"] = df_phoTag.Filter(f"n{photype}TightIso == {npho} && n{photype}BeamHalo == 0", f"eq{npho}{photype}TightIso")
            #    #regions[f"eq{npho}{photype}TightIso0BH"] = df_phoTag.Filter(f"n{photype}TightIso == {npho} && n{photype}BeamHalo", f"eq0{photype}BeamHalo")
            #    regions[f"eq{npho}{photype}MedIso"] = df_phoTag.Filter(f"n{photype}MedIso == {npho} && n{photype}TightIso == 0 && n{photype}BeamHalo == 0", f"eq{npho}{photype}MedIso0TightIso")

        #if MC, de#fine PB early/late regions and iso SRs
        if mc:
            #regions["tighIsoEESR"] = df_regs.Filter("regionIdx == 2 || regionIdx == 200",f"isoEESR_{ch_name}")
            #regions["tightIsoEBSR"] = df_regs.Filter("regionIdx == 4 || regionIdx == 400",f"isoEBSR_{ch_name}")
            #do inclusive object-multiplicity-defined region
            regions[ch_name] = df_new

        #regions["fail"] = df_regs.Filter("regionIdx < 0",f"fail_{ch_name}")
        #if regions["fail"].Count().GetValue() > 10:
        #    for i in range(1,7):
        #        regions[f"failMode_{i}"] = regions["fail"].Filter(f"regionIdx == -{i}",f"fail_mode{i}_{ch_name}")
        #        regions[f"failMode_{i*100}"] = regions["fail"].Filter(f"regionIdx == -{i*100}",f"fail_mode{i*100}_{ch_name}")
        return regions
   
    def do_ms_rs_cuts(self, df, reg_name):
        #in kinematic sidebands only do cuts on other obs
        if 'Ms' in reg_name:
            return df.Filter(f"rjr_Rs0 > {self._msrs_bins['BHCR']['rs'][0]}")
        elif 'Rs' in reg_name:
            return df.Filter(f"rjr_Ms0 > {self._msrs_bins['BHCR']['ms'][0]}")
        else:
            if any([reg_name not in binname for binname in self._msrs_bins]):
                reg_name_key = "*"
            else:
                reg_name_key = next(k for k in self._msrs_bins if k in reg_name)
            return df.Filter(f"rjr_Rs0 > {self._msrs_bins[reg_name_key]['rs'][0]}").Filter(f"rjr_Ms0 > {self._msrs_bins[reg_name_key]['ms'][0]}")


    def fill_region_hists(self, df, proc_name, reg_name, ch_name, h1d, h2d):
        if reg_name not in self._msrs_bins:
            msbins = array('d',[0, 10000])
            rsbins = array('d',[0, 1])
        else:
            reg_name_key = next(k for k in self._msrs_bins if k in reg_name)
            msbins = self._msrs_bins[reg_name_key]['ms']
            rsbins = self._msrs_bins[reg_name_key]['rs']
       
        if "CR" in reg_name:
            msmax = 4000
        else:
            msmax = 5000 
        n_msbins = len(msbins) - 1
        n_rsbins = len(rsbins) - 1
        if ch_name == reg_name:
            h1d.append(
                df.Histo1D(
	                (f"isoScore_{proc_name}_{ch_name}_{reg_name}", "", 50, 0, 1.01),
                    "photon_isoANNScore"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"PhotonTimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -5,5),
                    "photon_WTimeSig"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"photonPt_{proc_name}_{ch_name}_{reg_name}", "", 50, 0, 1000),
                    "photon_Pt"
                )
            )
        for photype in self._photypes:
            if photype not in reg_name and ch_name != reg_name:
                continue
            h1d.append(
                df.Histo1D(
	                (f"n{photype}_{proc_name}_{ch_name}_{reg_name}", "", 5,0,5),
                    f"n{photype}"
                )
            )
            if photype not in reg_name:
                continue
            h1d.append(
                df.Histo1D(
	                (f"n{photype}TightIso_{proc_name}_{ch_name}_{reg_name}", "", 3,0,3),
                    f"n{photype}TightIso"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"n{photype}BeamHalo_{proc_name}_{ch_name}_{reg_name}", "", 3,0,3),
                    f"n{photype}BeamHalo"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"{photype}TimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -5,5),
                    f"{photype}_WTimeSig"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"{photype}Pt_{proc_name}_{ch_name}_{reg_name}", "", 50, 0,1000),
                    f"{photype}_Pt"
                )
            )
            h1d.append(
                df.Histo1D(
	                (f"{photype}GenSusyId_{proc_name}_{ch_name}_{reg_name}", "", 60, -30,30),
                    "baseLinePhoton_SusyId"
                )
            )
         
   
        #h2d.append(
        #    df.Histo2D(
	    #        (f"yields_{proc_name}_{ch_name}_{reg_name}", ";Ms;Rs",
        #         n_msbins, msbins, n_rsbins, rsbins),
        #        "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
        #    )
        #)
      
        #h2d.append(
        #    df.Histo2D(
	    #        (f"MsRs_{proc_name}_{ch_name}_{reg_name}", ";Ms;Rs",
        #         25, 1000, msmax, 25, rsbins[0], 1.0),
        #        "rjr_Ms0", "rjr_Rs0", "evtFillWgt"
        #    )
        #)
        ##compressed kinematics
        #h2d.append(
        #    df.Filter("rjrIsr_PtIsr != -1").Histo2D( #remove events with ill-defined ISR trees
	    #        (f"RISRPtISR_{proc_name}_{ch_name}_{reg_name}", ";RISR;ptISR",
        #         25, 0, rsbins[-1]+0.01, 25, 50, 1500),
        #        "rjrIsr_RIsr", "rjrIsr_PtIsr", "evtFillWgt"
        #    )
        #)
   
        #h1d.append(
        #    df.Histo1D(
	    #        (f"Rs_{proc_name}_{ch_name}_{reg_name}", "", 25, rsbins[0], rsbins[-1]+0.01),
        #        "rjr_Rs0", "evtFillWgt"
        #    )
        #)
   	
        #h1d.append(
        #    df.Histo1D(
	    #        (f"Ms_{proc_name}_{ch_name}_{reg_name}", "", 25, msbins[0], msbins[-1]),
        #        "rjr_Ms0", "evtFillWgt"
        #    )
        #)
   

    



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
