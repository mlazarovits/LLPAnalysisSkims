from ROOT import RDataFrame, TChain, TFile, TH1, TH2, gInterpreter, std, gSystem
from tools import EfficiencyParser, WriteBFIJSON, MakeBFIRegions
import sys
import yaml
import awkward as ak
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
        #self._met_filters = (
        #    "(Flag_BadPFMuonFilter == 1 && "
        #    "Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && "
        #    "Flag_HBHENoiseFilter == 1 && "
        #    "Flag_HBHENoiseIsoFilter == 1 && "
        #    "Flag_ecalBadCalibFilter == 1 && "
        #    "Flag_eeBadScFilter == 1 && "
        #    "Flag_goodVertices == 1)"
        #)
        #matches yaml presel
        self._met_filters = "(Flag_MetFilters == 1)"
        self._basekin = "((rjr_Ms[0] >= 2000) && (rjr_Rs[0] >= 0.15))"
        #10000 is upperlimit (ie inclusive)
        self._kin_bins = {}
        self._kin_bins["*"] =   {"ms" : array("d", [0, 10000]), "rs": array("d", [0, 1.0])}
        self._kin_bins["BHCR"] = {"ms" : array("d", [2000, 7000]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["PBCR"] = {"ms" : array("d", [2000, 7000]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["PBSR"] = {"ms" : array("d", [2000, 7000]), "rs": array("d", [0.15, 0.3, 1.0])}
        self._kin_bins["MsCR"] = {"ms" : array("d", [0, self._kin_bins["BHCR"]['ms'][0]]), "rs": array("d", [0.15, 1.0])}
        self._kin_bins["RsCR"] = {"ms" : array("d", [0, 10000]), "rs": array("d", [self._kin_bins["BHCR"]['rs'][0], 1.0])}
        self._kin_bins["dxySigCR"] = {"ms" : array("d", [1000, 10000]), "rs": array("d", [0.15, 1.0])}
        self._yaml_path = ''


        self._dfkinbins = {}
        gInterpreter.Declare(
            """
                float LifetimeReweight(float ctau, float tau_old, float tau_new){
                    double inv_old = 1.0 / tau_old;
                    double inv_new = 1.0 / tau_new;
                    return float( (tau_old / tau_new) * std::exp(ctau * (inv_old - inv_new)) );
                }
            """
        )

        gInterpreter.Declare(
            """
            using ROOT::RVecF;
            using ROOT::RVecI;
            using ROOT::VecOps::RVec;
            using ROOT::VecOps::Any;
            using ROOT::VecOps::Sum;
            using ROOT::VecOps::Nonzero;
            using std::endl;
            using std::cout;
            int getRegionIdxValidation(const int nSVHad, const RVecF& SVHadDxySig, const RVecF& SVHadMass, const int nSVLep, const RVecF& SVLepDxySig, const int nPhotons, const RVecF& timesigs, const RVecF& bh_scores, const RVecF& iso_scores, const RVecF& photon_eta, const RVecF& photon_pt, const RVecI& nRJRJetsA, const RVecI& nRJRJetsB){
                double st_bhearly = -3;
                double st_isoearly = -2.5;
                double noniso_cutoff = 0.9; //>= is iso, < is noniso

                if(nSVLep > 0){
                    if(SVLepDxySig[0] < 50)
                        return 1;
                    else if(SVLepDxySig[0] >= 50 && SVLepDxySig[0] < 200)
                        return 2;
                    else{
                        //cout << "lep dxy sig " << SVLepDxySig[0] << endl;
                        return -1;
                    }
                }
                else{ //nSVLep == 0
                    if(nSVHad > 0){
                        if(nPhotons > 0){
                            auto mask = bh_scores < 0.185;
                            auto notbh_timesigs = timesigs[mask];
                            
                            auto lead_timesig = notbh_timesigs[notbh_timesigs >= 2.5][0];
                            if(lead_timesig >= 2.5){ //mixed regions - had SV + late signal photon
                                if(SVHadMass[0] > 15){
                                    if(SVHadDxySig[0] < 200)
                                       return 13;
                                    else if(SVHadDxySig[0] >= 200 && SVHadDxySig[0] < 800)
                                       return 14;
                                    else{
                                        //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                        return -2;
                                    }
                                }
                                else{
                                    return -3;
                                }
                            }
                            else{ //sv only region
                                if(SVHadMass[0] > 15){
                                    if(SVHadDxySig[0] < 200)
                                       return 3;
                                    else if(SVHadDxySig[0] >= 200 && SVHadDxySig[0] < 800)
                                        return 4;
                                    else{
                                        //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                        return -4;
                                    }
                                }
                                else{
                                    return -5;
                                } 
                            }
                        }
                        else{
                            if(SVHadMass[0] > 15){
                                if(SVHadDxySig[0] < 200)
                                   return 3;
                                else if(SVHadDxySig[0] >= 200 && SVHadDxySig[0] < 800)
                                   return 4;
                                else{
                                    //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                    return -6;
                                }
                            }
                            else{
                                return -7;
                            }
                        }
                    }
                    else{ //nSVLep == 0 && nSVHad == 0 //TODO - move any failed mass cut SV events here
                        if(nPhotons < 1)
                            return -8;
                        if(Any(timesigs < -2.5 || timesigs >= 2.5)){ //nonprompt
                            if(Any(bh_scores >= 0.95 && (bh_scores < 0.99))){ //bh regions
                                auto mask = ((bh_scores >= 0.95) && (bh_scores < 0.99));
                                auto bh_timesigs = timesigs[mask];
                                
                                auto lead_timesig = bh_timesigs[bh_timesigs < st_bhearly || bh_timesigs >= 2.5][0];
                                    //cout << "bh scores " << endl;
                                    //for(auto b : bh_scores) cout << b << endl;
                                    //cout << "timesigs " << endl;
                                    //for(auto t : timesigs) cout << t << endl;
                                    //cout << "bh timesigs " << endl;
                                    //for(auto t : bh_timesigs) cout << t << endl;
                                    //cout << "lead_timesig " << lead_timesig << endl; cout << endl; 
                                if(lead_timesig < st_bhearly){ //early bh cr
                                    return 5;
                                }
                                else if(lead_timesig >= 2.5) //late bh cr
                                    return 6;
                                else
                                    return -9;
                            }
                            else{ //no bh regions
                                if(Any(bh_scores < 0.95 && bh_scores >= 0.917252)){ //!bh regions
                                    auto mask = ((bh_scores < 0.95) && (bh_scores >= 0.917252));
                                    auto notbh_timesigs = timesigs[mask];
                                    
                                    auto lead_timesig = notbh_timesigs[(notbh_timesigs >= -10 && notbh_timesigs < st_bhearly) || notbh_timesigs >= 2.5][0];

                                    if(lead_timesig < st_bhearly && lead_timesig >= -10){ //early !bh cr
                                        return 7;
                                    }
                                    else if(lead_timesig >= 2.5){ //late !bhiso sr
                                        return 8;
                                    }
                                    else{ //recover these in np iso regions
                                        if(Any(iso_scores >= 0.5 && iso_scores < 0.6)){ //noniso np region
                                            auto mask = ((iso_scores < 0.6) && (iso_scores >= 0.5));
                                            auto mediso_timesigs = timesigs[mask];
                                            
                                            auto lead_timesig = mediso_timesigs[(mediso_timesigs >= st_bhearly && mediso_timesigs < st_isoearly) || mediso_timesigs >= 2.5][0];
                                            if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early noniso cr
                                                return 15;
                                            else if(lead_timesig >= 2.5) //late noniso cr
                                                return 16;
                                            else
                                                return -10; //TODO - move to prompt
                                        }
                                        else{ //iso np region
                                            auto mask = ((iso_scores >= 0.6) && (iso_scores < 0.8));
                                            auto tightiso_timesigs = timesigs[mask];
                                            
                                            auto lead_timesig = tightiso_timesigs[(tightiso_timesigs >= st_bhearly && tightiso_timesigs < st_isoearly) || tightiso_timesigs >= 2.5][0];
                                            if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early iso cr
                                                return 17;
                                            else if(lead_timesig >= 2.5) //late !bhiso sr
                                                return 8;
                                            else
                                                return -11; //TODO - move to prompt
                                        }
                                    }
                                }
                                else{ //np iso regions
                                    if(Any(iso_scores < 0.6 && iso_scores >= 0.5)){ //noniso np region
                                        auto mask = iso_scores < 0.6 && iso_scores >= 0.5;
                                        auto lead_idx = Nonzero(mask)[0];
                                        
                                        auto lead_timesig = timesigs[lead_idx];
                                        if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early noniso cr
                                            return 15;
                                        else if(lead_timesig >= 2.5) //late noniso cr
                                            return 16;
                                        else
                                            return -12;
                                    }
                                    else{ //iso np region
                                        auto mask = ((iso_scores >= 0.6) && (iso_scores < 0.8));
                                        auto iso_timesigs = timesigs[mask];
                                        //CAUTION - if no elements of masked vector satisfy giving masking condition, the mask will return exactly 0 
                                        auto lead_timesig = iso_timesigs[(iso_timesigs < st_isoearly && iso_timesigs >= st_bhearly) || iso_timesigs >= 2.5][0];
                                        if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early iso cr
                                            return 17;
                                        else if(lead_timesig >= 2.5) //late !bhiso sr
                                            return 8;
                                        else
                                            return -13;
                                    }

                                }

                            }
                        }
                        else{ //prompt
                            if(Any(bh_scores >= 0.95 && (bh_scores < 0.99))) //no bh photons!!
                                return -20;
                            if(nPhotons == 1){
                                if(!(( nRJRJetsA[0] >= 3 && nRJRJetsB[0] >= 2 ) || ( nRJRJetsA[0] >= 2 && nRJRJetsB[0] >= 3))){
                                    return -14;
                                }
                                else{
                                    //require barrel only
                                    if((photon_eta[0] < 1.479 && photon_eta[0] > -1.479)){
                                        if(photon_pt[0] > 100){
                                            if((iso_scores[0] < -0.000198*photon_pt[0] + 1.0098) && ( iso_scores[0] >= -0.000198*photon_pt[0] + 0.7698)){
                                                   return 9;
                                            }
                                            else if((iso_scores[0] < -0.000198*photon_pt[0] + 1.0188) && ( iso_scores[0] >= -0.000198*photon_pt[0] + 1.0098))
                                                   return 10;
                                            else
                                                return -15;
                                        }
                                        else{
                                            if((iso_scores[0] < 0.99) && (iso_scores[0] >= 0.75)){
                                                return 9;
                                            }
                                            else if((iso_scores[0] >= 0.99) && (iso_scores[0] < 0.999))
                                                return 10;
                                            else
                                                return -16;
                                        }
                                    }
                                    else{ //endcap
                                            if((iso_scores[0] < 0.8) && (iso_scores[0] >= 0.75)){
                                                return 9;
                                            }
                                            else if((iso_scores[0] < 0.9) && (iso_scores[0] >= 0.8))
                                                return 10;
                                            else
                                                return -17;
                                    }
                                }
                            }
                            else if(nPhotons == 2){
                                //count number of mediso + tightiso photons
                                int nMedIso = 0;
                                int nTightIso = 0;
                                int nNoIso = 0;
                                if((photon_eta[0] < 1.479 && photon_eta[0] > -1.479)){
                                    if(photon_pt[0] > 100){
                                        if(( iso_scores[0] < -0.0001*photon_pt[0] + 0.86 ) && (iso_scores[0] >= -0.0001*photon_pt[0] + 0.76 ))
                                            nMedIso++;
                                        else if((iso_scores[0] < -0.0001*photon_pt[0] + 0.96) && ( iso_scores[0] >= -0.0001*photon_pt[0] + 0.86  ))
                                            nTightIso++;
                                        else
                                            nNoIso++;

                                    }
                                    else{
                                        if((iso_scores[0] < 0.85 && iso_scores[0] >= 0.75))
                                            nMedIso++;
                                        else if( (iso_scores[0] >= 0.85) && (iso_scores[0] < 0.95))
                                            nTightIso++;
                                        else
                                            nNoIso++;
                                    }
                                }
                                else{ //recover endcap photons
                                     if((iso_scores[0] < 0.8 && iso_scores[0] >= 0.75))
                                         nMedIso++;
                                     else if((iso_scores[0] < 0.99) && (iso_scores[0] >= 0.8 ) )
                                         nTightIso++;
                                     else
                                         nNoIso++;
                                }
                                if(photon_eta[1] < 1.479 && photon_eta[1] > -1.479){
                                    if(photon_pt[1] > 100){
                                        if(( iso_scores[1] < -0.0001*photon_pt[1] + 0.86 ) && (iso_scores[1] >= -0.0001*photon_pt[1] + 0.76 ))
                                            nMedIso++;
                                        else if((iso_scores[1] < -0.0001*photon_pt[1] + 0.96) && ( iso_scores[1] >= -0.0001*photon_pt[1] + 0.86  ))
                                            nTightIso++;
                                        else
                                            nNoIso++;

                                    }
                                    else{
                                        if((iso_scores[1] < 0.85 && iso_scores[1] >= 0.75))
                                            nMedIso++;
                                        else if((iso_scores[1] >= 0.85) && (iso_scores[1] < 0.95))
                                            nTightIso++;
                                        else
                                            nNoIso++;
                                    }

                                }
                                else{ //recover endcap photons
                                     if((iso_scores[1] < 0.8 && iso_scores[1] >= 0.75))
                                         nMedIso++;
                                     else if( (iso_scores[1] < 0.95) && (iso_scores[1] >= 0.8 ))
                                         nTightIso++;
                                     else
                                         nNoIso++;

                                }
                                if(nMedIso > 0)
                                    return 11;
                                else if(nMedIso == 0 && nTightIso > 0)
                                    return 12;
                                else
                                    return -18; 
                            }
                            else return -19;
                        }
                    }
                }
            }
            """
        )

        gInterpreter.Declare(
            """
            using ROOT::RVecF;
            using ROOT::RVecI;
            using ROOT::VecOps::RVec;
            using ROOT::VecOps::Any;
            using ROOT::VecOps::Sum;
            using ROOT::VecOps::Nonzero;
            using std::endl;
            using std::cout;
            int getRegionIdx(const int nSVHad, const RVecF& SVHadDxySig, const RVecF& SVHadMass, const int nSVLep, const RVecF& SVLepDxySig, const int nPhotons, const RVecF& timesigs, const RVecF& bh_scores, const RVecF& iso_scores, const RVecF& photon_eta, const RVecF& photon_pt, const RVecI& nRJRJetsA, const RVecI& nRJRJetsB){
                double st_bhearly = -3;
                double st_isoearly = -2.5;
                double noniso_cutoff = 0.9; //>= is iso, < is noniso

                if(nSVLep > 0){
                    if(SVLepDxySig[0] < 800)
                        return 1;
                    else if(SVLepDxySig[0] >= 800)
                        return 2;
                    else{
                        //cout << "lep dxy sig " << SVLepDxySig[0] << endl;
                        return -1;
                    }
                }
                else{ //nSVLep == 0
                    if(nSVHad > 0){
                        if(nPhotons > 0){
                            auto mask = bh_scores < 0.185;
                            auto notbh_timesigs = timesigs[mask];
                            
                            auto lead_timesig = notbh_timesigs[notbh_timesigs >= 2.5][0];
                            if(lead_timesig >= 2.5){ //mixed regions - had SV + late signal photon
                                if(SVHadMass[0] > 15){
                                    if(SVHadDxySig[0] < 800)
                                       return 13;
                                    else if(SVHadDxySig[0] > 800)
                                       return 14;
                                    else{
                                        //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                        return -2;
                                    }
                                }
                                else{
                                    return -3;
                                }
                            }
                            else{ //sv only region
                                if(SVHadMass[0] > 15){
                                    if(SVHadDxySig[0] < 800)
                                       return 3;
                                    else if(SVHadDxySig[0] > 800)
                                        return 4;
                                    else{
                                        //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                        return -4;
                                    }
                                }
                                else{
                                    return -5;
                                } 
                            }
                        }
                        else{
                            if(SVHadMass[0] > 15){
                                if(SVHadDxySig[0] < 800)
                                   return 3;
                                else if(SVHadDxySig[0] > 800)
                                   return 4;
                                else{
                                    //cout << "had dxy sig " << SVHadDxySig[0] << endl;
                                    return -6;
                                }
                            }
                            else{
                                return -7;
                            }
                        }
                    }
                    else{ //nSVLep == 0 && nSVHad == 0 //TODO - move any failed mass cut SV events here
                        if(nPhotons < 1)
                            return -8;
                        if(Any(timesigs < -2.5 || timesigs >= 2.5)){ //nonprompt
                            if(Any(bh_scores >= 0.917252)){ //bh regions
                                auto mask = bh_scores >= 0.917252;
                                auto bh_timesigs = timesigs[mask];
                                
                                auto lead_timesig = bh_timesigs[bh_timesigs < st_bhearly || bh_timesigs >= 2.5][0];
                                    //cout << "bh scores " << endl;
                                    //for(auto b : bh_scores) cout << b << endl;
                                    //cout << "timesigs " << endl;
                                    //for(auto t : timesigs) cout << t << endl;
                                    //cout << "bh timesigs " << endl;
                                    //for(auto t : bh_timesigs) cout << t << endl;
                                    //cout << "lead_timesig " << lead_timesig << endl; cout << endl; 
                                if(lead_timesig < st_bhearly){ //early bh cr
                                    return 5;
                                }
                                else if(lead_timesig >= 2.5) //late bh cr
                                    return 6;
                                else
                                    return -9;
                            }
                            else{ //no bh regions
                                if(Any(bh_scores < 0.185)){ //!bh regions
                                    auto mask = bh_scores < 0.185;
                                    auto notbh_timesigs = timesigs[mask];
                                    
                                    auto lead_timesig = notbh_timesigs[(notbh_timesigs >= -10 && notbh_timesigs < st_bhearly) || notbh_timesigs >= 2.5][0];

                                    if(lead_timesig < st_bhearly && lead_timesig >= -10){ //early !bh cr
                                        return 7;
                                    }
                                    else if(lead_timesig >= 2.5){ //late !bhiso sr
                                        return 8;
                                    }
                                    else{ //recover these in np iso regions
                                        if(Any(iso_scores < noniso_cutoff && iso_scores >= 0.4)){ //noniso np region
                                            auto mask = (iso_scores < noniso_cutoff && iso_scores >= 0.4);
                                            auto mediso_timesigs = timesigs[mask];
                                            
                                            auto lead_timesig = mediso_timesigs[(mediso_timesigs >= st_bhearly && mediso_timesigs < st_isoearly) || mediso_timesigs >= 2.5][0];
                                            if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early noniso cr
                                                return 15;
                                            else if(lead_timesig >= 2.5) //late noniso cr
                                                return 16;
                                            else
                                                return -10; //TODO - move to prompt
                                        }
                                        else{ //iso np region
                                            auto mask = iso_scores >= noniso_cutoff;
                                            auto tightiso_timesigs = timesigs[mask];
                                            
                                            auto lead_timesig = tightiso_timesigs[(tightiso_timesigs >= st_bhearly && tightiso_timesigs < st_isoearly) || tightiso_timesigs >= 2.5][0];
                                            if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early iso cr
                                                return 17;
                                            else if(lead_timesig >= 2.5) //late !bhiso sr
                                                return 8;
                                            else
                                                return -11; //TODO - move to prompt
                                        }
                                    }
                                }
                                else{ //np iso regions
                                    if(Any(iso_scores < noniso_cutoff && iso_scores >= 0.4)){ //noniso np region
                                        auto mask = iso_scores < noniso_cutoff && iso_scores >= 0.4;
                                        auto lead_idx = Nonzero(mask)[0];
                                        
                                        auto lead_timesig = timesigs[lead_idx];
                                        if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early noniso cr
                                            return 15;
                                        else if(lead_timesig >= 2.5) //late noniso cr
                                            return 16;
                                        else
                                            return -12;
                                    }
                                    else{ //iso np region
                                        auto mask = iso_scores >= noniso_cutoff;
                                        auto iso_timesigs = timesigs[mask];
                                        //CAUTION - if no elements of masked vector satisfy giving masking condition, the mask will return exactly 0 
                                        auto lead_timesig = iso_timesigs[(iso_timesigs < st_isoearly && iso_timesigs >= st_bhearly) || iso_timesigs >= 2.5][0];
                                        if(lead_timesig >= st_bhearly && lead_timesig < st_isoearly) //early iso cr
                                            return 17;
                                        else if(lead_timesig >= 2.5) //late !bhiso sr
                                            return 8;
                                        else
                                            return -13;
                                    }

                                }

                            }
                        }
                        else{ //prompt
                            if(Any(bh_scores >= 0.917252)) //no bh photons!!
                                return -20;
                            if(nPhotons == 1){
                                if(!(( nRJRJetsA[0] >= 3 && nRJRJetsB[0] >= 2 ) || ( nRJRJetsA[0] >= 2 && nRJRJetsB[0] >= 3))){
                                    return -14;
                                }
                                else{
                                    //require barrel only
                                    if((photon_eta[0] < 1.479 && photon_eta[0] > -1.479)){
                                        if(photon_pt[0] > 100){
                                            if((iso_scores[0] < -0.000198*photon_pt[0] + 1.0188) && ( iso_scores[0] >= -0.000198*photon_pt[0] + 0.7698)){
                                                   return 9;
                                            }
                                            else if(iso_scores[0] >= -0.000198*photon_pt[0] + 1.0188)
                                                   return 10;
                                            else
                                                return -15;
                                        }
                                        else{
                                            if((iso_scores[0] < 0.999) && (iso_scores[0] >= 0.75)){
                                                return 9;
                                            }
                                            else if((iso_scores[0] >= 0.999))
                                                return 10;
                                            else
                                                return -16;
                                        }
                                    }
                                    else{ //endcap
                                            if((iso_scores[0] < 0.9) && (iso_scores[0] >= 0.75)){
                                                return 9;
                                            }
                                            else if((iso_scores[0] >= 0.9))
                                                return 10;
                                            else
                                                return -17;
                                    }
                                }
                            }
                            else if(nPhotons == 2){
                                //count number of mediso + tightiso photons
                                int nMedIso = 0;
                                int nTightIso = 0;
                                int nNoIso = 0;
                                if((photon_eta[0] < 1.479 && photon_eta[0] > -1.479)){
                                    if(photon_pt[0] > 100){
                                        if(( iso_scores[0] < -0.0001*photon_pt[0] + 0.96 ) && (iso_scores[0] >= -0.0001*photon_pt[0] + 0.76 ))
                                            nMedIso++;
                                        else if((iso_scores[0] >= -0.0001*photon_pt[0] + 0.96))
                                            nTightIso++;
                                        else
                                            nNoIso++;

                                    }
                                    else{
                                        if((iso_scores[0] < 0.95 && iso_scores[0] >= 0.75))
                                            nMedIso++;
                                        else if( (iso_scores[0] >= 0.95) )
                                            nTightIso++;
                                        else
                                            nNoIso++;
                                    }
                                }
                                else{ //recover endcap photons
                                     if((iso_scores[0] < 0.99 && iso_scores[0] >= 0.75))
                                         nMedIso++;
                                     else if( (iso_scores[0] >= 0.99) )
                                         nTightIso++;
                                     else
                                         nNoIso++;
                                }
                                if(photon_eta[1] < 1.479 && photon_eta[1] > -1.479){
                                    if(photon_pt[1] > 100){
                                        if(( iso_scores[1] < -0.0001*photon_pt[1] + 0.96 ) && (iso_scores[1] >= -0.0001*photon_pt[1] + 0.76 ))
                                            nMedIso++;
                                        else if((iso_scores[1] >= -0.0001*photon_pt[1] + 0.96))
                                            nTightIso++;
                                        else
                                            nNoIso++;

                                    }
                                    else{
                                        if((iso_scores[1] < 0.95 && iso_scores[1] >= 0.75))
                                            nMedIso++;
                                        else if( (iso_scores[1] >= 0.95) )
                                            nTightIso++;
                                        else
                                            nNoIso++;
                                    }

                                }
                                else{ //recover endcap photons
                                     if((iso_scores[1] < 0.95 && iso_scores[1] >= 0.75))
                                         nMedIso++;
                                     else if( (iso_scores[1] >= 0.95) )
                                         nTightIso++;
                                     else
                                         nNoIso++;

                                }
                                if(nMedIso > 0)
                                    return 11;
                                else if(nMedIso == 0 && nTightIso > 0)
                                    return 12;
                                else
                                    return -18; 
                            }
                            else return -19;
                        }
                    }
                }
            }
            """
        )



    def set_yaml_path(self, val = False):
        if val:
            self._yaml_path = 'BigGuy_NonCompressed_Validation_AnalysisConfig_PromptBHDelayedContributions.yaml'
        else:
            self._yaml_path = 'BigGuy_NonCompressed_FullRegions_AnalysisConfig_PromptBHDelayedContributions.yaml'
            #'BigGuy_NonCompressed_FullRegions_AnalysisConfig_PhotonFlags.yaml'

    def get_kin_bins(self, val = False):
        with open(self._yaml_path,'r') as f:
            data = yaml.safe_load(f)
        SVlepbins = {"bLep00" : data['bLep00'][0], "bLep01" : data['bLep01'][0], "bLep10" : data['bLep10'][0], "bLep11" : data['bLep11'][0]}
        self._dfkinbins["Ch1CRGeLep1"] = SVlepbins
        self._dfkinbins["Ch2SRGeLep1"] = SVlepbins
        self._dfkinbins["Ch2CRGeLep1"] = SVlepbins
        
        SVhadbins = {"bHad00" : data['bHad00'][0], "bHad01" : data['bHad01'][0], "bHad10" : data['bHad10'][0], "bHad11" : data['bHad11'][0], "bHad21" : data['bHad21'][0]}
        self._dfkinbins["Ch3CRGeHad1"] = SVhadbins
        self._dfkinbins["Ch4SRGeHad1"] = SVhadbins
        self._dfkinbins["Ch4CRGeHad1"] = SVhadbins

        #if validation do this separately for Ms-Rs
        if args.val:
            split_yaml_ms = "BigGuy_NonCompressed_Validation_AnalysisConfig_DelayedPhotonMsShapes.yaml"
            split_yaml_rs = "BigGuy_NonCompressed_Validation_AnalysisConfig_DelayedPhotonRsShapes.yaml"
            with open(split_yaml_ms,'r') as f:
                data_split_bins_ms = yaml.safe_load(f)
            with open(split_yaml_rs,'r') as f:
                data_split_bins_rs = yaml.safe_load(f)
            PhoDelayedbins = {"DelayedBin00" : data_split_bins_ms['DelayedBin00'][0], "DelayedBin01" : data_split_bins_rs['DelayedBin01'][0], "DelayedBin10" : data_split_bins_ms['DelayedBin10'][0], "DelayedBin11" : data_split_bins_rs['DelayedBin11'][0]}
            self._dfkinbins["Ch5CRgeq1PhoBHEarly"] = PhoDelayedbins
            self._dfkinbins["Ch6CRgeq1PhoBHLate"] = PhoDelayedbins
            self._dfkinbins["Ch7CRgeq1PhoNotBHEarly"] = PhoDelayedbins
            self._dfkinbins["Ch8SRgeq1PhoNotBHLateTightIso"] = PhoDelayedbins
            self._dfkinbins["Ch8CRgeq1PhoNotBHLateTightIso"] = PhoDelayedbins
            self._dfkinbins["Ch15CRgeq1PhoMedIsoEarly"] = PhoDelayedbins 
            self._dfkinbins["Ch16CRgeq1PhoMedIsoLate"] = PhoDelayedbins
            self._dfkinbins["Ch17CRgeq1PhoTightIsoEarly"] = PhoDelayedbins
        else:
            PhoDelayedbins = {"DelayedBin00" : data['DelayedBin00'][0], "DelayedBin01" : data['DelayedBin01'][0], "DelayedBin10" : data['DelayedBin10'][0], "DelayedBin11" : data['DelayedBin11'][0]}
            self._dfkinbins["Ch5CRgeq1PhoBHEarly"] = PhoDelayedbins
            self._dfkinbins["Ch6CRgeq1PhoBHLate"] = PhoDelayedbins
            self._dfkinbins["Ch7CRgeq1PhoNotBHEarly"] = PhoDelayedbins
            self._dfkinbins["Ch8SRgeq1PhoNotBHLateTightIso"] = PhoDelayedbins
            self._dfkinbins["Ch8CRgeq1PhoNotBHLateTightIso"] = PhoDelayedbins
            self._dfkinbins["Ch15CRgeq1PhoMedIsoEarly"] = PhoDelayedbins 
            self._dfkinbins["Ch16CRgeq1PhoMedIsoLate"] = PhoDelayedbins
            self._dfkinbins["Ch17CRgeq1PhoTightIsoEarly"] = PhoDelayedbins

        eq1PhoPromptbins = {"eq1b00" : data['eq1b00'][0], "eq1b10" : data['eq1b10'][0], "eq1b01" : data['eq1b01'][0], "eq1b11" : data['eq1b11'][0], "eq1b02" : data['eq1b02'][0], "eq1b12" : data['eq1b12'][0]}
        self._dfkinbins["Ch9CReq1PhoMedIsoPrompt"] = eq1PhoPromptbins 
        self._dfkinbins["Ch10SReq1PhoTightIsoPrompt"] = eq1PhoPromptbins 
        self._dfkinbins["Ch10CReq1PhoTightIsoPrompt"] = eq1PhoPromptbins 

        eq2PhoPromptbins = {"eq2b00" : data['eq2b00'][0], "eq2b10" : data['eq2b10'][0], "eq2b01" : data['eq2b01'][0], "eq2b11" : data['eq2b11'][0]}
        self._dfkinbins["Ch11CReq2PhoMedIsoPrompt"] = eq2PhoPromptbins 
        self._dfkinbins["Ch12SReq2PhoTightIsoPrompt"] = eq2PhoPromptbins
        self._dfkinbins["Ch12CReq2PhoTightIsoPrompt"] = eq2PhoPromptbins
        
        SVPhoDelayedbins = {"SVDelayedBin00":data['SVDelayedBin00'][0], "SVDelayedBin01" :data['SVDelayedBin01'][0]} 
        self._dfkinbins["Ch13CRgeq1SVLowDxygeq1PhoNotBHLate"] = SVPhoDelayedbins 
        self._dfkinbins["Ch14SRgeq1SVHighDxygeq1PhoNotBHLate"] = SVPhoDelayedbins
        self._dfkinbins["Ch14CRgeq1SVHighDxygeq1PhoNotBHLate"] = SVPhoDelayedbins

        self._dfkinbins["Ch8SRgeq1PhoNotBHLate"] = PhoDelayedbins
        self._dfkinbins["Ch18CRGJetsLate"] = PhoDelayedbins
        self._dfkinbins["Ch19CRRsLate"] = PhoDelayedbins
        #self._dfkinbins["Ch18SRgeq1PhoTightIsoLate"] = PhoDelayedbins

    def do_kin_bins(self, region_dict):
        binned_regions = {}
        for regname, reg_df in region_dict.items():
            if "Ch" not in regname:
                continue
            kinbins = self._dfkinbins[regname]
            for binname, kbin in kinbins.items():
                binidx = binname[-2:]
                if binidx == "00" and "SR" in regname and "NotBHLateTightIso" not in regname: #dont do this for ABCD shape
                    regbin = regname.replace("SR","CR")+binidx 
                else:
                    regbin = regname+binidx
                binned_regions[regbin] = reg_df.Filter(kbin)
        return binned_regions
            

    def define_regions_ifelse(self, df, mc, val = False, unblind = False):
        if val:
            df = df.Define("regionIdx",f"getRegionIdxValidation(SV_nHadronic, HadronicSV_dxySig, HadronicSV_mass, SV_nLeptonic, LeptonicSV_dxySig, nBaseLinePhotons, baseLinePhoton_WTimeSig, baseLinePhoton_beamHaloCNNScore, baseLinePhoton_isoANNScore, baseLinePhoton_Eta, baseLinePhoton_Pt, rjrNJetsJa, rjrNJetsJb)")
        else:
            df = df.Define("regionIdx",f"getRegionIdx(SV_nHadronic, HadronicSV_dxySig, HadronicSV_mass, SV_nLeptonic, LeptonicSV_dxySig, nBaseLinePhotons, baseLinePhoton_WTimeSig, baseLinePhoton_beamHaloCNNScore, baseLinePhoton_isoANNScore, baseLinePhoton_Eta, baseLinePhoton_Pt, rjrNJetsJa, rjrNJetsJb)")
        MsCR = "(rjr_Ms0 <= 2000)"
        regions = {
            "Ch1CRGeLep1" : df.Filter("regionIdx == 1","Ch1CRGeLep1"),
            "Ch3CRGeHad1" : df.Filter("regionIdx == 3","Ch3CRGeHad1"),
            "Ch5CRgeq1PhoBHEarly" : df.Filter("regionIdx == 5","Ch5CRgeq1PhoBHEarly"),
            "Ch6CRgeq1PhoBHLate" : df.Filter("regionIdx == 6","Ch6CRgeq1PhoBHLate"),
            "Ch7CRgeq1PhoNotBHEarly" : df.Filter("regionIdx == 7","Ch7CRgeq1PhoNotBHEarly"),
            "Ch9CReq1PhoMedIsoPrompt" : df.Filter("regionIdx == 9","Ch9CReq1PhoMedIsoPrompt"),
            "Ch11CReq2PhoMedIsoPrompt" : df.Filter("regionIdx == 11","Ch11CReq2PhoMedIsoPrompt"),
            "Ch13CRgeq1SVLowDxygeq1PhoNotBHLate" : df.Filter("regionIdx == 13","Ch13CRgeq1SVLowDxygeq1PhoNotBHLate"),
            "Ch15CRgeq1PhoMedIsoEarly" : df.Filter("regionIdx == 15","Ch15CRgeq1PhoMedIsoEarly"),
            "Ch16CRgeq1PhoMedIsoLate" : df.Filter("regionIdx == 16","Ch16CRgeq1PhoMedIsoLate"),
            "Ch17CRgeq1PhoTightIsoEarly" : df.Filter("regionIdx == 17","Ch17CRgeq1PhoTightIsoEarly"),
            "Ch18CRGJetsLate" : df.Filter("((nBaseLinePhotons == 1 && (baseLinePhoton_WTimeSig[0] > 2.5 && baseLinePhoton_GJetsCR[0] == 1)) || (nBaseLinePhotons == 2 && ( baseLinePhoton_WTimeSig[0] > 2.5 && baseLinePhoton_GJetsCR[0] == 1) || ( baseLinePhoton_WTimeSig[1] > 2.5 && baseLinePhoton_GJetsCR[1] == 1)) )","Ch18CRGJetsLate"),
            "Ch19CRRsLate" : df.Filter("((rjr_Ms[0] >= 2000) && (rjr_Rs[0] < 0.15))","Ch19CRRsLate"),
        }
 
        if unblind or mc:
            regions["Ch2SRGeLep1"] = df.Filter("regionIdx == 2","Ch2SRGeLep1")
            regions["Ch4SRGeHad1"] = df.Filter("regionIdx == 4","Ch4SRGeHad1")
            regions["Ch8SRgeq1PhoNotBHLateTightIso"] = df.Filter("regionIdx == 8","Ch8SRgeq1PhoNotBHLateTightIso")
            #regions["Ch8SRgeq1PhoNotBHLate"] = df.Filter("regionIdx == 8","Ch8SRgeq1PhoNotBHLateTightIso")
            #regions["Ch18SRgeq1PhoTightIsoLate"] = df.Filter("regionIdx == 18","Ch18SRgeq1PhoNotBHTightIsoLate")
            regions["Ch10SReq1PhoTightIsoPrompt"] = df.Filter("regionIdx == 10","Ch10SReq1PhoTightIsoPrompt")
            regions["Ch12SReq2PhoTightIsoPrompt"] = df.Filter("regionIdx == 12","Ch12SReq2PhoTightIsoPrompt")
            regions["Ch14SRgeq1SVHighDxygeq1PhoNotBHLate"] = df.Filter("regionIdx == 14","Ch14SRgeq1SVHighDxygeq1PhoNotBHLate")
            regions["presel"] = df
        #create regions of 00 for data (!unblind and !mc)
        else:
            regions["Ch2SRGeLep1"] = df.Filter(f"(regionIdx == 2) && {self._dfkinbins['Ch2SRGeLep1']['bLep00']}","Ch2SRGeLep1")
            regions["Ch4SRGeHad1"] = df.Filter(f"(regionIdx == 4) && {self._dfkinbins['Ch4SRGeHad1']['bHad00']}","Ch4SRGeHad1")
            regions["Ch8SRgeq1PhoNotBHLateTightIso"] = df.Filter(f"(regionIdx == 8) && {self._dfkinbins['Ch8SRgeq1PhoNotBHLateTightIso']['DelayedBin00']}","Ch8SRgeq1PhoNotBHLateTightIso")
            regions["Ch10SReq1PhoTightIsoPrompt"] = df.Filter(f"(regionIdx == 10) && {self._dfkinbins['Ch10SReq1PhoTightIsoPrompt']['eq1b00']}","Ch10SReq1PhoTightIsoPrompt")
            regions["Ch12SReq2PhoTightIsoPrompt"] = df.Filter(f"(regionIdx == 12) && {self._dfkinbins['Ch12SReq2PhoTightIsoPrompt']['eq2b00']}","Ch12SReq2PhoTightIsoPrompt")
            regions["Ch14SRgeq1SVHighDxygeq1PhoNotBHLate"] = df.Filter(f"(regionIdx == 14) && {self._dfkinbins['Ch14SRgeq1SVHighDxygeq1PhoNotBHLate']['SVDelayedBin00']}","Ch14SRgeq1SVHighDxygeq1PhoNotBHLate")
        regions["failSel"] = df.Filter("regionIdx < 1","failSel")
        for i in range(-20,0):
            regions[f"failSel_{i}"] = df.Filter(f"regionIdx == {i}")#,f"failSel_{i}")       
        #regions["MsCR"] = df.Filter(MsCR,"MsCR")
        return regions

    def check_overlap(self, df):
            overlaps = {}
            #for i in range(-19, 15):
            for i in range(9, 10):
                for j in range(-20, i):
                    check_df = df.Filter(f"regionIdx == {i} && regionIdx == {j}")
                    overlaps[f"region{i}_region{j}"] = check_df.Count()
            return overlaps


    def define_regions_yaml(self, df, ch_name, mc, unblind):
        #assuming baseline and cleaning are the same across yamls
        baseline = None
        cleaning = None
        print("Parsing yaml",self._yaml_path)
        with open(self._yaml_path,'r') as f:
            data = yaml.safe_load(f)
        baseline = data['baseline_cuts'][0]
        cleaning = data['Cleaning'][0]
        kin = data['kin'][0]
        if 'regions' not in data.keys():
            print("Regions not defined in yaml",yaml_path,"Please define in yaml and rerun")
            exit()
        yaml_regions = data['regions']
        regions = {}
        regions["presel"] = df.Filter(f"({baseline}) && ({cleaning}) && ({kin})")
        for name, regdef in yaml_regions.items():
            cutstring = self.flatten_cutstring(regdef)
            if not mc and "SR" in name and not unblind:
                regions[name] = df.Filter("false",name)
            else:   
                regions[name] = df.Filter(cutstring, name)
        return regions

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

        Ch3CRGeHad = f"{noSVLep} && {SVgeHad1_CR} && {noNonPromptPhotons}"
        Ch4SRGeHad = f"{noSVLep} && {SVgeHad1_SR} && {noNonPromptPhotons}"


        ##photon preselection
        noSV = "((SV_nHadronic==0) && (SV_nLeptonic==0))"
        ge1Pho = "(nBaseLinePhotons > 0)"
        noBHPhotons = "(( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.917252) || (nBaseLinePhotons == 2 && baseLinePhoton_beamHaloCNNScore[0] < 0.917252 && baseLinePhoton_beamHaloCNNScore[1] < 0.917252) ) )"

        ##delayed photon regions
        ge1BHPhoEarly = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] < -2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] < -2.5) || (baseLinePhoton_beamHaloCNNScore[1] >= 0.917252 && baseLinePhoton_WTimeSig[1] < -2.5)) ) )"
        ge1BHPhoLate = "(nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] >= 2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] >= 0.917252 && baseLinePhoton_WTimeSig[0] >= 2.5) || (baseLinePhoton_beamHaloCNNScore[1] >= 0.917252 && baseLinePhoton_WTimeSig[1] >= 2.5)) )"
        ge1NotBHPhoEarly = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] < -2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] < -2.5) || (baseLinePhoton_beamHaloCNNScore[1] < 0.185 && baseLinePhoton_WTimeSig[1] < -2.5)) ) )"
        ge1NotBHPhoLate = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5) || (baseLinePhoton_beamHaloCNNScore[1] < 0.185 && baseLinePhoton_WTimeSig[1] >= 2.5)) ) )"
        ge1NotBHPhoLateTightIso = "( (nBaseLinePhotons == 1 && baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5 && (baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.000198*baseLinePhoton_Pt[0] + 1.0188)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.999))) || ( nBaseLinePhotons == 2 && ((baseLinePhoton_beamHaloCNNScore[0] < 0.185 && baseLinePhoton_WTimeSig[0] >= 2.5 && (((baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.0001*baseLinePhoton_Pt[0] + 0.96)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.95)))) || (baseLinePhoton_beamHaloCNNScore[1] < 0.185 && baseLinePhoton_WTimeSig[1] >= 2.5 && (baseLinePhoton_Pt[1] > 100 && (baseLinePhoton_isoANNScore[1] >= -0.0001*baseLinePhoton_Pt[1] + 0.96)) || (baseLinePhoton_Pt[1] <= 100 && baseLinePhoton_isoANNScore[1] >= 0.95))))))"
        Ch5CRgeq1PhoBHEarly = f"{noSV} && {ge1Pho} && {ge1BHPhoEarly}"        
        Ch6CRgeq1PhoBHLate = f"{noSV} && {ge1Pho} && {ge1BHPhoLate}"        
        Ch7CRgeq1PhoNotBHEarly = f"{noSV} && {noBHPhotons} && {ge1Pho} && {ge1NotBHPhoEarly}"        
        Ch8SRgeq1PhoNotBHLate = f"{noSV} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoLate}"  
        Ch8SRgeq1PhoNotBHLateTightIso = f"{noSV} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoLateTightIso}"      
      

        ##prompt photon regions
        eq1PhoMedIsoPrompt = "(nBaseLinePhotons == 1  && ((baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] < -0.000198*baseLinePhoton_Pt[0] + 1.0188) && ( baseLinePhoton_isoANNScore[0] >= -0.000198*baseLinePhoton_Pt[0] + 0.7698)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] < 0.999) && (baseLinePhoton_isoANNScore[0] >= 0.75)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479)))"
        eq1PhoTightIsoPrompt = "(nBaseLinePhotons == 1 && (baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.000198*baseLinePhoton_Pt[0] + 1.0188)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.999)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479))"
        eq2PhoMedIsoPrompt = "(nBaseLinePhotons == 2 && ((baseLinePhoton_Pt[0] > 100 && ( baseLinePhoton_isoANNScore[0] < -0.0001*baseLinePhoton_Pt[0] + 0.96 ) && (baseLinePhoton_isoANNScore[0] >= -0.0001*baseLinePhoton_Pt[0] + 0.76 )) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] < 0.95 && baseLinePhoton_isoANNScore[0] >= 0.75)) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479 )) || ((baseLinePhoton_Pt[1] > 100 && (baseLinePhoton_isoANNScore[1] < -0.0001*baseLinePhoton_Pt[1] + 0.96) && ( (baseLinePhoton_isoANNScore[1] >= -0.0001*baseLinePhoton_Pt[1] + 0.76))) || (baseLinePhoton_Pt[1] <= 100 && baseLinePhoton_isoANNScore[1] < 0.95 && baseLinePhoton_isoANNScore[1] >= 0.75) && (baseLinePhoton_Eta[1] < 1.479 && baseLinePhoton_Eta[1] > -1.479 )) )"
        eq2PhoTightIsoPrompt = "(nBaseLinePhotons == 2 && (((baseLinePhoton_Pt[0] > 100 && (baseLinePhoton_isoANNScore[0] >= -0.0001*baseLinePhoton_Pt[0] + 0.96)) || (baseLinePhoton_Pt[0] <= 100 && (baseLinePhoton_isoANNScore[0] >= 0.95))) && (baseLinePhoton_Eta[0] < 1.479 && baseLinePhoton_Eta[0] > -1.479)) && ((baseLinePhoton_Pt[1] > 100 && (baseLinePhoton_isoANNScore[1] >= -0.0001*baseLinePhoton_Pt[1] + 0.96)) || (baseLinePhoton_Pt[1] <= 100 && baseLinePhoton_isoANNScore[1] >= 0.95) && (baseLinePhoton_Eta[1] < 1.479 && baseLinePhoton_Eta[1] > -1.479)))" 

        Ch9CReq1PhoMedIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq1PhoMedIsoPrompt}"
        Ch10SReq1PhoTightIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq1PhoTightIsoPrompt}"
        Ch11CReq2PhoMedIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq2PhoMedIsoPrompt}"
        Ch12SReq2PhoTightIsoPrompt = f"{noSV} && {ge1Pho} && {noBHPhotons} && {noNonPromptPhotons} && {eq2PhoTightIsoPrompt}"
    
        ##mixed SV delayed photon regions
        Ch13CRgeq1SVgeq1PhoBHEarly =  f"{noSVLep} && {SVgeHad1_CR} && {ge1Pho} && {ge1BHPhoEarly}"
        Ch14CRgeq1SVgeq1PhoBHLate = f"{noSVLep} && {SVgeHad1_CR} && {ge1Pho} && {ge1BHPhoLate}"
        Ch15CRgeq1SVgeq1PhoNotBHEarly = f"{noSVLep} && {SVgeHad1_CR} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoEarly}"
        Ch16SRgeq1SVgeq1PhoNotBHLate = f"{noSVLep} && {SVgeHad1_CR} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoLate}"
        Ch16SRgeq1SVgeq1PhoNotBHLateMsCR = f"(rjr_Ms0 <= 2000) && {noSVLep} && {SVgeHad1_CR} && {ge1Pho} && {noBHPhotons} && {ge1NotBHPhoLate}"

        MsCR = "(rjr_Ms0 <= 2000)"
        regions = {
            "Ch1CRGeLep1" : df.Filter(Ch1CRGeLep,"Ch1CRGeLep"),
            "Ch3CRGeHad1" : df.Filter(Ch3CRGeHad,"Ch3CRGeHad1"),
            "Ch5CRgeq1PhoBHEarly" : df.Filter(Ch5CRgeq1PhoBHEarly,"Ch5CRgeq1PhoBHEarly"),
            "Ch6CRgeq1PhoBHLate" : df.Filter(Ch6CRgeq1PhoBHLate,"Ch6CRgeq1PhoBHLate"),
            "Ch7CRgeq1PhoNotBHEarly" : df.Filter(Ch7CRgeq1PhoNotBHEarly,"Ch7CRgeq1PhoNotBHEarly"),
            "Ch9CReq1PhoMedIsoPrompt" : df.Filter(Ch9CReq1PhoMedIsoPrompt,"Ch9CReq1PhoMedIsoPrompt"),
            "Ch11CReq2PhoMedIsoPrompt" : df.Filter(Ch11CReq2PhoMedIsoPrompt,"Ch11CReq2PhoMedIsoPrompt"),
            "Ch13CRgeq1SVgeq1PhoBHEarly" : df.Filter(Ch13CRgeq1SVgeq1PhoBHEarly,"Ch13CRgeq1SVgeq1PhoBHEarly"),
            "Ch14CRgeq1SVgeq1PhoBHLate" : df.Filter(Ch14CRgeq1SVgeq1PhoBHLate,"Ch14CRgeq1SVgeq1PhoBHLate"),
            "Ch15CRgeq1SVgeq1PhoNotBHEarly" : df.Filter(Ch15CRgeq1SVgeq1PhoNotBHEarly,"Ch15CRgeq1SVgeq1PhoNotBHEarly"),
            "Ch16SRgeq1SVgeq1PhoNotBHLate" : df.Filter(Ch16SRgeq1SVgeq1PhoNotBHLateMsCR,"Ch16SRgeq1SVgeq1PhoNotBHLate"),
            "MsCR" : df.Filter(MsCR,"MsCR")
        }
        if mc:
            regions["Ch2SRGeLep1"] = df.Filter(Ch2SRGeLep,"Ch2SRGeLep1");
            regions["Ch4SRGeHad1"] = df.Filter(Ch4SRGeHad,"Ch4SRGeHad1");
            regions["Ch8SRgeq1PhoNotBHLate"] = df.Filter(Ch8SRgeq1PhoNotBHLate,"Ch8SRgeq1PhoNotBHLate")
            regions["Ch8SRgeq1PhoNotBHLateTightIso"] = df.Filter(Ch8SRgeq1PhoNotBHLateTightIso,"Ch8SRgeq1PhoNotBHLateTightIso")
            regions["Ch10SReq1PhoTightIsoPrompt"] = df.Filter(Ch10SReq1PhoTightIsoPrompt,"Ch10SReq1PhoTightIsoPrompt")
            regions["Ch12SReq2PhoTightIsoPrompt"] = df.Filter(Ch12SReq2PhoTightIsoPrompt,"Ch12SReq2PhoTightIsoPrompt")
            regions["Ch16SRgeq1SVgeq1PhoNotBHLate"] = df.Filter(Ch16SRgeq1SVgeq1PhoNotBHLate,"Ch16SRgeq1SVgeq1PhoNotBHLate")
            regions["presel"] = df
        return regions

    #includes barrel photon iso presel 
    def apply_preselection(self, df, compressed=False):
        presel = f"{self._metcut} && {self._triggers} && {self._met_filters}"
        if not compressed:
            presel += f"&& {self._ptscut} && {self._basekin}"
        return df.Filter(presel, "presel")
       
    def flatten_cutstring(self, cutlist):
        flatlist = ak.flatten(cutlist,axis = None) 
        return ' && '.join(flatlist)


    def fill_region_hists(self, df, proc_name, reg_name, ch_name, h1d, h2d, val = False, unblind = False):
        if "failSel" in reg_name:
            return
        if not unblind and "MET" in proc_name and "SR" in ch_name:
            return
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
        if ch_name == reg_name or reg_name == "MsCR" or val:
            highptcut = "200"
            vhighptcut = "500"
            df = df.Define("isoScoreHighPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt >= {highptcut} && baseLinePhoton_Pt < {vhighptcut}]").Define("isoScoreLowPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt < {highptcut}]").Define("isoScoreVeryHighPt",f"baseLinePhoton_isoANNScore[baseLinePhoton_Pt >= {vhighptcut}]")
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonBeamHaloScore_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonBeamHaloScore_{proc_name}", 80, -0.07, 1.07),
                    "baseLinePhoton_beamHaloCNNScore"
                )
            )
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScore_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScore_{proc_name}", 80, 0.15, 1.07),
                    "baseLinePhoton_isoANNScore"
                )
            )
            df = df.Define("baseLinePhoton_isoANNScore_NotBHTagged","baseLinePhoton_isoANNScore[baseLinePhoton_beamHaloCNNScore < 0.185]")
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonIsoScoreNotBHTagged_{proc_name}_{ch_name}_{reg_name}", f"baseLinePhotonIsoScoreNotBHTagged_{proc_name}", 80, 0.15, 1.07),
                    "baseLinePhoton_isoANNScore_NotBHTagged"
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
            timesig_upperend = 20
            if ch_name == "presel":
                timesig_upperend = 2.
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhotonTimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -20,timesig_upperend),
                    "baseLinePhoton_WTimeSig"
                )
            )
            df = df.Define("baseLinePhoton_notBHTimeSig","baseLinePhoton_WTimeSig[baseLinePhoton_beamHaloCNNScore < 0.185]")
            h1d.append(
                df.Histo1D(
                    (f"baseLinePhoton_notBHTimeSig_{proc_name}_{ch_name}_{reg_name}", "", 50, -10,2.5),
                    "baseLinePhoton_notBHTimeSig"
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
                    (f"nBaseLinePhotons_{proc_name}_{ch_name}_{reg_name}", "", 5,0,5),
                    f"nBaseLinePhotons"
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
                 (f"isoScoreBHScore_{proc_name}_{ch_name}_{reg_name}", ";isoScore;bhScore",
                     50, 0, 1.01, 50, 0, 1.01),
                    "baseLinePhoton_isoANNScore", "baseLinePhoton_beamHaloCNNScore"
                )
            )
            df = df.Define("delayedIsoScore","baseLinePhoton_isoANNScore[baseLinePhoton_WTimeSig > 2.5]").Define("delayedBHScore","baseLinePhoton_beamHaloCNNScore[baseLinePhoton_WTimeSig > 2.5]")
            h2d.append(
                df.Histo2D(
                 (f"delayedIsoScoreBHScore_{proc_name}_{ch_name}_{reg_name}", ";isoScore;bhScore",
                     50, 0, 1.01, 50, 0, 1.01),
                    "delayedIsoScore", "delayedBHScore"
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
            if "MET" in proc:
                mc = False
            files = fileprocessor.GetFiles(proc, mGl, mN2, mN1, ctau)
            compressed = False
            if int(mGl) - int(mN1) <= 200:
                compressed = True
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
                    lumi_factor = args.lumi 
                    #if "SMS" in proc:
                    #    lumi_factor = args.lumi 
                    #else: #bkg MC for 2018 only rn
                    #    lumi_factor = args.lumi / 67.9
                else:
                    lumi_factor = 1   
                print(proc,"lumi",lumi_factor) 
                df1 = df1.Redefine("evtFillWgt",f"evtFillWgt*{lumi_factor}")
    
                #do individual presel cuts here so they are printed out
                df_metcut = df1.Filter(self._metcut,self._metcut)
                if not compressed:
                    df_pts = df1.Filter(self._ptscut, self._ptscut)
                    df_basekin = df1.Filter(self._basekin,"basekin") 
                df_triggers = df1.Filter(self._triggers,"triggers")
                df_filters = df1.Filter(self._met_filters,"met_filters")
                #do all presel cuts
                
                #channels = self.define_channels(df_presel)
                #ch_name = "ge1blpho"
                #df_ge1pho = df_presel.Filter("nBaseLinePhotons > 0 && SV_nLeptonic == 0 && SV_nHadronic == 0")
                ch_name = "presel"
                df_list = []
                if(args.yaml):
                    df_list.append(df1)
                    regions = self.define_regions_yaml(df1, ch_name, mc, args.unblind)
                else:
                    df_presel = self.apply_preselection(df1,compressed)
                    df_list.append(df_presel)
                    regions = self.define_regions_ifelse(df_presel, mc, args.val)

                #overlaps = self.check_overlap(df_presel)

                #get weighted evt counts
                wt_counts = {}
                for reg in regions:
                    wt_counts[reg] = regions[reg].Sum("evtFillWgt")

                report = df00.Report()
                #for key, count in overlaps.items():
                #    if(count.GetValue() > 0):
                #        print("overlap",key,count.GetValue())
                # select cuts
                lines = self._eff_parser.report2str(report, wt_counts)
                denom_info = self._eff_parser.get_denom_line(lines, 'presel')
                total_eff = 0
                total_ana_eff = 0
                total_sr_eff = 0
                for line in lines:
                    parsed_eff = self._eff_parser.parse_eff_line(line, denom_info)
                    if parsed_eff is None:
                        continue
                    else:
                        if(show_output):
                            print(line)
                            #print("parsed_eff",parsed_eff)
                    if parsed_eff[0] != "presel" and parsed_eff[0] != "failSel":
                        #print(parsed_eff[0],f"{parsed_eff[2]:.2f}%")
                        total_eff += parsed_eff[2]
                        if "failSel" not in parsed_eff[0]:
                            total_ana_eff += parsed_eff[2]
                            if "SR" in parsed_eff[0]:
                                total_sr_eff += parsed_eff[2]
                    selected_data[infilename].append(parsed_eff)
                # write LaTeX table
                outfile = f"{procstr}_eff_table"
                if args.ofilename_extra is not None:
                    outfile += f"_{args.ofilename_extra}"
                if outfile[-1] == "_":
                    outfile = outfile[:-1]
                outfile += ".tex"
                print(f"total_eff: {total_eff:.2f}%, total analysis eff: {total_ana_eff:.2f}%, total SR eff: {total_sr_eff:.2f}%") 
                self._eff_parser.write_latex_table(outfile, selected_data, args.lumi, total_eff)
                print("Wrote efficiencies for process",proc,"to",outfile)

    def ProcessCtauStr(self, ctau):
        if "p" in ctau:
            ctau = ctau.replace("p",".")
        return float(ctau) * 100 #go from m to cm

    def GetProcessName(self, proc, mGl = None, mN2 = None, mN1 = None, ctau = None):
        procstr = proc
        if proc == "METPD":
            procstr = "METFullRunII_RunIII"   
        if "rwtct" in proc:
            newct = proc[proc.find("rwtct")+5:]
            procstr = proc[:proc.find("_rwtct")]
            oldct = procstr[procstr.find("ct")+2:]
            procstr = procstr.replace(oldct, newct)
        return procstr

# -------------------------
# Main analysis
# -------------------------

    def runRJRAnalysis(self, args, ofilename_extra: str = ""):
        if args.unblind:
            print("UNBLINDING ANALYSIS - YIELDS IN SRs FOR DATA WILL BE SAVED!!!!")
        procs = args.proc
        mGl = args.mGl
        mN2 = args.mN2
        mN1 = args.mN1
        ctaus = args.ctau
        if ctaus is None:
            ctaus = [None] 
        new_taus = args.new_taus
        if new_taus is None:
            new_taus = ctaus
        if "SMS_gogoGZ" in procs:
            for ctaupair in zip(ctaus, new_taus):
                procs.append(f"SMS_gogoGZmGl{mGl}mN2{mN2}mN1{mN1}ct{ctaupair[0]}_rwtct{ctaupair[1]}")
            procs.remove("SMS_gogoGZ")
 
        fileprocessor = FileProcessor()
   
        self.set_yaml_path(args.val)
        if(args.yaml):
            print("Using yaml",self._yaml_path)
        else:
            print("Using IfElse logic")

 
        ofilename = ""
        hists1d, hists2d = [], []
        jsondict = {}
        jsonname = "BFI"
        if ofilename_extra:
            jsonname += f"_{ofilename_extra}"
        if args.val:
            jsonname += "_Validation"
        jsonname += ".json"
        self.get_kin_bins()
        for proc in procs:
            procstr = self.GetProcessName(proc)
            files = []
            compressed = False
            ctaupair = (None, None)
            if "SMS" in proc: #assume only one mass point given
                ctaupair = (proc[proc.find("ct")+2:proc.find("_rwtct")], proc[proc.find("rwtct")+5:])
                files += fileprocessor.GetFiles("SMS_gogoGZ", mGl, mN2, mN1, ctaupair[0])
                if int(mGl) - int(mN1) <= 200:
                    compressed = True
                if(len(files) < 1):
                    print("No files found for ",procstr)
                    continue
            else: 
                files += fileprocessor.GetFiles(proc)
            if len(files) == 0:
                print(f"No files found for proc {proc} with ctau={args.ctau}, mGl={args.mGl}, mN2={args.mN2}, mN1={args.mN1}. Going to next process")
                continue
            #remove underscores from procstr
            procstr = procstr.replace("_","")
            if "SMS" not in procstr:
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
            #print("There are",df00.Filter("rjr_Ms.size() == 0 && nPhotons > 1").Count().GetValue(),"events without RJR info")
            df = df00.Filter("rjr_Rs.size() > 0 && rjr_Ms.size() > 0") #in case these have size 0, can lead to undefined behavior
            df1 = (
                df.Define("rjr_Rs0", "rjr_Rs[0]")
                  .Define("rjr_Ms0", "rjr_Ms[0]")
            )
            if(mc):
                lumi_factor = args.lumi 
                df1 = df1.Redefine("evtFillWgt",f"evtFillWgt*{lumi_factor}")
            if "SMS" in proc and ctaupair[1] != ctaupair[0]:
                tau_old = self.ProcessCtauStr(ctaupair[0])
                tau_new = self.ProcessCtauStr(ctaupair[1])
                if tau_old != tau_new:
                    print(f"Reweighting ctau {tau_old} to {tau_new}")
                    #define xa xb ct weights
                    df1 = df1.Define("xawt",f"LifetimeReweight(Xa_ctau, {tau_old}, {tau_new})").Define("xbwt",f"LifetimeReweight(Xb_ctau, {tau_old}, {tau_new})")
                    df1 = df1.Redefine("evtFillWgt",f"evtFillWgt*xawt*xbwt")
            #do individual presel cuts here so they are printed out
            df_metcut = df1.Filter(self._metcut,self._metcut)
            df_pts = df1.Filter(self._ptscut, self._ptscut)
            df_triggers = df1.Filter(self._triggers,"triggers")
            df_filters = df1.Filter(self._met_filters,"met_filters")
    
            #do all presel cuts
            df_presel = self.apply_preselection(df1,compressed)
            df_list = [df_presel]
            
            #channels = self.define_channels(df_presel)
            #ch_name = "ge1blpho"
            #df_ge1pho = df_presel.Filter("nBaseLinePhotons > 0 && SV_nLeptonic == 0 && SV_nHadronic == 0")
            ch_name = "presel"
            if(args.yaml):
                regions = self.define_regions_yaml(df_presel, ch_name, mc, args.unblind)
            else:
                regions = self.define_regions_ifelse(df_presel, mc, args.val, args.unblind)


            binned_regions = self.do_kin_bins(regions)
            procjsondict = MakeBFIRegions(jsondict, procstr, binned_regions)
            jsondict = jsondict | procjsondict

            for reg_name, df_reg in regions.items():
                #print("  doing region", reg_name)
                #define columns with lowest Rs/Ms cuts
                self.fill_region_hists(
                    df_reg, procstr, reg_name, ch_name,
                    hists1d, hists2d, args.val, args.unblind
                )
                df_list.append(df_reg)
            report = df00.Report() 
            region = "presel"
            #region = "Ch13CRgeq1SVLowDxygeq1PhoNotBHLate"
            #region = "Ch15CRgeq1PhoMedIsoEarly"
            if region in regions.keys() and args.testLogic:
                print(region)
                regions[region] = (regions[region]
                    .Define("eta","baseLinePhoton_Eta").Define("iso_score","baseLinePhoton_isoANNScore").Define("bh_score","baseLinePhoton_beamHaloCNNScore")
                    .Define("pt","baseLinePhoton_Pt").Define("tsig","baseLinePhoton_WTimeSig")
                    .Define("LowPtEB","(pt < 100) && (eta < 1.479 && eta > -1.479) && (iso_score < 0.95 && iso_score >= 0.75)") 
                    .Define("Eq2HighPtEB","(pt >= 100) && (eta < 1.479 && eta > -1.479) && ( iso_score < -0.0001*pt + 0.96 ) && (iso_score >= -0.0001*pt + 0.76 )") 
                    .Define("Eq2Tight","(pt >= 100) && (eta < 1.479 && eta > -1.479) && ( iso_score >= -0.0001*pt + 0.96 )") 
                    .Define("Eq1HighPtEB","(pt >= 100) && (eta < 1.479 && eta > -1.479) && (iso_score < -0.000198*pt + 1.0188) && ( iso_score >= -0.000198*pt + 0.7698)") 
                    .Define("LeadEE","!(eta[0] < 1.479 && eta[0] > -1.479) && (iso_score[0] < 0.99)") 
                    .Define("SubleadEE","!(eta[1] < 1.479 && eta[1] > -1.479) && (iso_score[1] < 0.95)")
                    .Define("leadSVDxySig","HadronicSV_dxySig[0]").Define("leadSVMass","HadronicSV_mass[0]")
                    .Define("pholatesig","passNPhoGe1SelectionLateSignal")
                    .Define("earlytight","passNPhoGe1SelectionEarlyTightIsoCR")
                    .Define("earlytight0","passNPhoGe1SelectionEarlyTightIso0NotBHCR")
                    .Define("latemed","passNPhoGe1SelectionLateMedIsoCR")
                    .Define("latemed0","passNPhoGe1SelectionLateMedIso0NotBHCR")
                    .Define("npmediso_tsig","npmediso_tagged_lead_timesig")
                    .Define("yamlch9","((SV_nHadronic==0) && (SV_nLeptonic==0)) && (passNPhoEq1SelectionPromptMedIsoValCR == 1) && ( ( rjrNJetsJa[0] >= 3 && rjrNJetsJb[0] >= 2 ) || ( rjrNJetsJa[0] >= 2 && rjrNJetsJb[0] >= 3 ) )")
                    .Define("yamlch10","((SV_nHadronic==0) && (SV_nLeptonic==0)) && (passNPhoEq1SelectionPromptTightIsoValSR == 1) && ( ( rjrNJetsJa[0] >= 3 && rjrNJetsJb[0] >= 2 ) || ( rjrNJetsJa[0] >= 2 && rjrNJetsJb[0] >= 3 ) )")
                    .Define("yamlch11","((SV_nHadronic==0) && (SV_nLeptonic==0)) && (passNPhoEq2SelectionPromptMedIsoValCR == 1)")
                    .Define("yamlch13","(SV_nLeptonic == 0) && ((SV_nHadronic >= 1) && (HadronicSV_dxySig[0] < 200) && (HadronicSV_mass[0] > 15)) && (passNPhoGe1SelectionLateSignal == 1)")
                    .Define("yamlch15","((SV_nHadronic==0) && (SV_nLeptonic==0)) && ((passNPhoGe1SelectionEarlyMedIsoValCR == 1) || (passNPhoGe1SelectionEarlyMedIso0NotBHValCR == 1))")
                    .Define("yamlch16","((SV_nHadronic==0) && (SV_nLeptonic==0)) && ((passNPhoGe1SelectionLateMedIsoValCR == 1) || (passNPhoGe1SelectionLateMedIso0NotBHValCR == 1))")
                    .Define("yamlch17","((SV_nHadronic==0) && (SV_nLeptonic==0)) && ((passNPhoGe1SelectionEarlyTightIsoValCR == 1) || (passNPhoGe1SelectionEarlyTightIso0NotBHValCR == 1))")
                )
                chnum = "11"
                dispcols = [f"yamlch{chnum}","tsig","iso_score","val_npmediso_tagged_lead_timesig"]
                if not args.yaml:
                    dispcols.append("regionIdx")
                    regions[region].Filter(f"(yamlch{chnum}) && (regionIdx != {chnum})").Display(dispcols,20).Print()
                    print(f"events not in ifelse in yamlch{chnum}",regions[region].Filter(f"(yamlch{chnum}) && (regionIdx != {chnum})").Count().GetValue()) 
                    regions[region].Filter(f"!yamlch{chnum} && (regionIdx == {chnum})").Display(dispcols,20).Print()
                    print(f"events in ifelse not in yamlch{chnum}",regions[region].Filter(f"!yamlch{chnum} && (regionIdx == {chnum})").Count().GetValue()) 
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
        if any("SMS_gogoGZ" in proc for proc in args.proc):
            ofilename += "SMSgogoGZ" 
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

        WriteBFIJSON(jsonname, jsondict, args.unblind)
        print("Wrote BFI",jsonname)

if __name__ == "__main__":
    #import kerberos credentials to conda env if not already there
    kerb = os.getenv("KRB5CCNAME")
    if(kerb is None):
        print("Setting kerebos credentials")
        os.environ["KRB5CCNAME"] = "API:"
    parser = argparse.ArgumentParser(
        description="Run RJR analysis"
    )


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
        nargs='+',
        help="ctau(s) for SMS process"
    )
    
    parser.add_argument(
        "--new_taus",
        default=None,
        nargs='+',
        help="new taus for reweighting"
    )

    parser.add_argument(
        "--lumi",
        default=1,
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
    parser.add_argument(
        '--yaml',
        help='use definitions from yaml file',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--val',
        help='do validation regions',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--unblind',
        help='unblind regions',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--testLogic',
        help='test region logic',
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
