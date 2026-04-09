import subprocess

def double2str(dbl):
	ret = str(dbl)
	ret = ret.replace(".","p")
	return ret`

lead_iso_scores = [0.5, 0.7, 0.8, 0.9]
sublead_iso_scores = [0.5, 0.7, 0.8, 0.9]
lead_pts = [70,100]
sublead_pts = [70,100]


yaml_preamble = 
"""
analysis:
  name: "PhoOnly_FullRegions_NonCompressed_AnalysisConfig_CutMaster"
  luminosity: 165.27
  output_json: "PhoOnly_FullRegions_NonCompressed_AnalysisConfig_CutMaster.json"

samples:
  #backgrounds: []#"Wjets", "Zjets", "QCD", "Top", "Gjets", "Box"]
  signals: ["gogoGZ"] 
  data: ["MET16","MET17","MET18","MET23"]

baseline_cuts: &baseline_cuts
  - "(selCMet > 150)" 

Cleaning: &Cleaning
  - "(Flag_MetFilters == 1) && ( rjrPTS[0] < 150) && ( (Trigger_PFMET120_PFMHT120_IDTight == 1)||(Trigger_PFMET120_PFMHT120_IDTight_PFHT60 == 1)||(Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight)||(Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60))"
################################ SV VETO Strings ##################
#exactly 1 hadronic SR SV
noSV: &noSV
  - "SV_nHadronic==0 || (SV_nHadronic == 1 && HadronicSV_dxySig[0] < 1000) || (SV_nHadronic == 2 && HadronicSV_dxySig[0] < 1000 && HadronicSV_dxySig[1] < 1000) || SV_nLeptonic==0 || (SV_nLeptonic == 1 && LeptonicSV_dxySig[0] < 1000) || (SV_nLeptonic == 2 && LeptonicSV_dxySig[0] < 1000 && LeptonicSV_dxySig[1] < 1000)"
"""

#function1: write analysis config
#set new cutstring from arrays
#write and save analysis config with unique name
#return list of analysis config file names

def WriteAnalysisConfigs():
	lead_iso_score = lead_iso_scores[-1]
	sublead_iso_score = sublead_iso_scores[-1]
	lead_pt = lead_pts[0]
	sublead_pt = sublead_pts[0]

	leadiso = double2str(lead_iso_score) 
	subleadiso = double2str(sublead_iso_score) 
	leadiso = double2str(lead_pt) 
	subleadiso = double2str(sublead_pt) 

	filename = f"PhoOnly_PromptRegions_Uncompressed_AnalysisConfig_leadiso-{leadiso}_subleadiso-{subleadiso}_leadpt-{leadpt}_subleadpt-{subleadpt}.yaml"
	with open(filename,"w") as f:
		f.write(yaml_preamble)	
	print("wrote to",filename)	


#function2: run BFIs
#for config in analysis_config_file_names : run BFI with specified config
#return list of BFI json names

#function3: parse BFI jsons
#for BFI in BFIjsons :
#	get data in anchorch_anchorbin
#	get data in buoych_anchorbin
#	calculate transfer factor tf = buoych_anchorbin / anchorch_anchorbin
#	get data in anchorch_SRbin
#	calculate prediction buoych_SRbin_pred = tf * anchorch_SRbin
#	get signal in buoych_SRbin
#	calculate signf = buoych_SRbin / sqrt(buoych_SRbin_pred)
#	print out: sigpt signf sigyield bkgpred anchorbin leading_iso_score sublead_iso_score lead_pt sublead_pt
