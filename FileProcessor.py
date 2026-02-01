import os
from ROOT import TFile
import re
import numpy as np

class FileProcessor():
    def __init__(self):
        self._procFiles = {}
        path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v46/"
        
        self._procFiles["SMS_gogoGZ"] = [
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2150_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2150_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-500_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1300_mN1-1000_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2350_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2400_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2300_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2350_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1000_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1500_ct0p5_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-500_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2150_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2150_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2200_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1500_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2300_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2350_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2350_ct0p1_rjrskim.root",
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2400_ct0p1_rjrskim.root",
        ]
        
        self._procFiles["SMS_test_prompt"] = [
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root"
        ]
        
        self._procFiles["SMS_test_delayed"] = [
            path+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-500_ct0p5_rjrskim.root"
        ]
        
	path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47/"
        self._procFiles["SMS_gogoG"] = [
            path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2500_mN2-1500_mN1-1000_ct0p1_rjrskim.root",
	]
        
        self._procFiles["WJets"] = [
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-100To200_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-1200To2500_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-200To400_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-2500ToInf_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-400To600_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-600To800_rjrskim.root",
        	path+"WJets_R18_SVIPM100_v31_WJetsToLNu_HT-800To1200_rjrskim.root"
        ]
        
        self._procFiles["ZJets"] = [
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-100To200_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-1200To2500_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-200To400_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-2500ToInf_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-400To600_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-600To800_rjrskim.root",
        	path+"ZJets_R18_SVIPM100_v31_ZJetsToNuNu_HT-800To1200_rjrskim.root"
        ]
        
        self._procFiles["TTXJets"] = [
        	path+"TTXJets_R18_SVIPM100_v31_TGJets_rjrskim.root",
        	path+"TTXJets_R18_SVIPM100_v31_TTGJets_rjrskim.root",
        	path+"TTXJets_R18_SVIPM100_v31_TTJets_rjrskim.root",
        	path+"TTXJets_R18_SVIPM100_v31_ttWJets_rjrskim.root",
        	path+"TTXJets_R18_SVIPM100_v31_ttZJets_rjrskim.root"
        ]
        
        path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v45/"
        self._procFiles["METPD16"] = [
            path+"MET_R16_SVIPM100_v31_MET_AOD_Run2016B_rjrskim_v45.root",
            path+"MET_R16_SVIPM100_v31_MET_AOD_Run2016C_rjrskim_v45.root",
            path+"MET_R16_SVIPM100_v31_MET_AOD_Run2016D_rjrskim_v45.root",
            path+"MET_R16_SVIPM100_v31_MET_AOD_Run2016G_rjrskim_v45.root",
            path+"MET_R16_SVIPM100_v31_MET_AOD_Run2016H_rjrskim_v45.root",
        ]
        
        self._procFiles["METPD17"] = [
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017A_rjrskim_v45.root",
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017B_rjrskim_v45.root",
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017C_rjrskim_v45.root",
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017D_rjrskim_v45.root",
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017E_rjrskim_v45.root",
            path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017F_rjrskim_v45.root",
        ]
        
        self._procFiles["METPD18"] = [
            path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018A_rjrskim_v45.root",
            path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018B_rjrskim_v45.root",
            path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018C_rjrskim_v45.root",
            path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018D_rjrskim_v45.root",
        ]

        self._procFiles["METPD_test"] = [
           path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018B_rjrskim_v45.root"
        ]
        
        path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v46/"
        self._procFiles["QCD"] = [
            path+"QCD_R18_SVIPM100_v31_QCD_HT1000to1500_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT100to200_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT1500to2000_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT2000toInf_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT200to300_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT300to500_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT500to700_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT50to100_rjrskim.root",
            path+"QCD_R18_SVIPM100_v31_QCD_HT700to1000_rjrskim.root",
        ]
        
        self._procFiles["GJets"] = [
            path+"GJets_R18_SVIPM100_v31_GJets_HT-40To100_rjrskim.root",
            path+"GJets_R18_SVIPM100_v31_GJets_HT-100To200_rjrskim.root",
            path+"GJets_R18_SVIPM100_v31_GJets_HT-200To400_rjrskim.root",
            path+"GJets_R18_SVIPM100_v31_GJets_HT-400To600_rjrskim.root",
            path+"GJets_R18_SVIPM100_v31_GJets_HT-600ToInf_rjrskim.root",
            path+"DiPJBox_R18_SVIPM100_v31_DiPhotonJetsBox_M40_80-sherpa_rjrskim.root",
            path+"DiPJBox_R18_SVIPM100_v31_DiPhotonJetsBox_MGG-0to40_rjrskim.root",
            path+"DiPJBox_R18_SVIPM100_v31_DiPhotonJetsBox_MGG-80toInf_rjrskim.root"
        ]

    def field(self, name, val):
        # match exact value if given, otherwise wildcard number
        return rf"{name}-{val}" if val is not None else rf"{name}-[^_]+"

    def GetFiles(self, proc, mGl = None, mN2 = None, mN1 = None, ctau = None):
        if proc not in self._procFiles.keys() and "PD" not in proc:
            print("Process",proc,"does not have files set. Currently set samples are",self._procFiles.keys())
            exit()

        if "SMS" not in proc and "PD" not in proc:
            return self._procFiles[proc]

        if "PD" in proc:
            if any(char.isdigit() for char in proc): #year specified
                return self._procFiles[proc]
            else: #all years
                #keys = [key in self._procFiles.keys() if proc in key]
                ret_files = []
                for key, val in self._procFiles.items():
                	if proc not in key:
                		continue
                	ret_files += val
                return ret_files

        pattern = (
            rf"mGl-{self.field('mGl', mGl).split('-',1)[1]}_"
            rf"mN2-{self.field('mN2', mN2).split('-',1)[1]}_"
            rf"mN1-{self.field('mN1', mN1).split('-',1)[1]}_"
            rf"ct{ctau if ctau is not None else '[^_]+'}"
        )
        regex = re.compile(pattern)
        if "SMS" in proc: #do regex matching
            return [f for f in self._procFiles[proc] if regex.search(f)]

    



