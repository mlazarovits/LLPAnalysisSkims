import os
from ROOT import RDataFrame, TChain, TFile, gInterpreter
from ComparisonPlotting import create_comparison_canvas, dist, parse_data_name
from RJRAnalysis import RJRAnalysis
gInterpreter.Declare(
	"""
	using ROOT::RVecF;
	using ROOT::RVec;
	RVec<int> getIsoIdxs(const RVecF& scores, const RVecF& etas, const float ee_score_thresh, const float eb_score_thresh){
	RVec<int> retidxs;
		if(scores.size() < 1){
			return retidxs;
		}
	if(etas.size() != scores.size()){
		return retidxs;
	}
	for(int i = 0; i < scores.size(); i++){
		float score = scores[i];
		float eta = etas[i];
		if(eta < 1.479 && eta > -1.479){
			if(scores[i] > eb_score_thresh){
				retidxs.push_back(i);
			}
		}
		else{ //endcap
			if(scores[i] > ee_score_thresh){
				retidxs.push_back(i);
			}
		
		}
		
	}
	return retidxs;
		}
"""
)

def AddHists(df, hists, dfs, branches, histname = ""):
	#inclusive photons
	eta = df.Histo1D((branches["PhoEta"]+histname,branches["PhoEta"]+histname,50,-3,3),branches["PhoEta"])
	hists.append(eta)
	photime = df.Histo1D((branches["PhoWTime"]+histname,branches["PhoWTime"]+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),branches["PhoWTime"])
	hists.append(photime)
	timesig = df.Histo1D( (branches["PhoWTimeSig"]+histname,branches["PhoWTimeSig"]+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),branches["PhoWTimeSig"])
	hists.append(timesig)
	bh_score = df.Histo1D((branches["Pho_beamHaloCNNScore"]+histname,branches["Pho_beamHaloCNNScore"]+histname,50,-0.01,1.01),branches["Pho_beamHaloCNNScore"])
	hists.append(bh_score)
	iso_score = df.Histo1D((branches["Pho_isoANNScore"]+histname,branches["Pho_isoANNScore"]+histname,50,-0.01,1.01),branches["Pho_isoANNScore"])
	hists.append(iso_score)	
	pt = df.Histo1D((branches['PhoPt']+histname,branches['PhoPt']+histname,50,dist["PhoPt"]["xmin"],dist["PhoPt"]["xmax"]),branches["PhoPt"])
	hists.append(pt)


	#lead photon
	df_leadpho = df.Define("leadTime",f"{branches['PhoWTime']}[0]").Define("leadTimeSig",f"{branches['PhoWTimeSig']}[0]").Define("leadBHScore",f"{branches['Pho_beamHaloCNNScore']}[0]").Define("leadIsoScore",f"{branches['Pho_isoANNScore']}[0]").Define("leadPt",f"{branches['PhoPt']}[0]")
	dfs.append(df_leadpho)
	leadphotime = df_leadpho.Histo1D((f"{branches['PhoWTime']}Lead"+histname,"{branches['PhoWTime]'}Lead"+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),"leadTime")
	hists.append(leadphotime)
	leadtimesig = df_leadpho.Histo1D( (f"{branches['PhoWTimeSig']}Lead"+histname,f"{branches['PhoWTimeSig']}Lead"+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"leadTimeSig")
	hists.append(leadtimesig)
	bh_scorelead = df_leadpho.Histo1D((f"{branches['Pho_beamHaloCNNScore']}Lead"+histname,f"{branches['Pho_beamHaloCNNScore']}Lead"+histname,50,-0.01,1.01),"leadBHScore")
	hists.append(bh_scorelead)
	iso_scorelead = df_leadpho.Histo1D((f"{branches['Pho_isoANNScore']}Lead"+histname,f"{branches['Pho_isoANNScore']}Lead"+histname,50,-0.01,1.01),"leadIsoScore")
	hists.append(iso_scorelead)	
	ptlead = df_leadpho.Histo1D((f"{branches['PhoPt']}Lead"+histname,f"{branches['PhoPt']}Lead"+histname,50,dist["PhoPt"]["xmin"],dist["PhoPt"]["xmax"]),"leadPt")
	hists.append(ptlead)

	#sublead photon
	df_ge2pho = df.Filter(f"{branches['nPhotons']} > 1").Define("subleadTime",f"{branches['PhoWTime']}[1]").Define("subleadTimeSig",f"{branches['PhoWTimeSig']}[1]").Define("subleadBHScore",f"{branches['Pho_beamHaloCNNScore']}[1]").Define("subleadIsoScore",f"{branches['Pho_isoANNScore']}[1]").Define("subleadPt",f"{branches['PhoPt']}[1]")
	dfs.append(df_ge2pho)
	subleadphotime = df_ge2pho.Histo1D((f"{branches['PhoWTime']}Sublead"+histname,f"{branches['PhoWTime']}Sublead"+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),"subleadTime")
	hists.append(subleadphotime)
	subleadtimesig = df_ge2pho.Histo1D( (f"{branches['PhoWTimeSig']}Sublead"+histname,f"{branches['PhoWTimeSig']}Sublead"+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"subleadTimeSig")
	hists.append(subleadtimesig)
	bh_scoresublead = df_ge2pho.Histo1D((f"{branches['Pho_beamHaloCNNScore']}Sublead"+histname,f"{branches['Pho_beamHaloCNNScore']}Sublead"+histname,50,-0.01,1.01),"subleadBHScore")
	hists.append(bh_scoresublead)
	iso_scoresublead = df_ge2pho.Histo1D((f"{branches['Pho_isoANNScore']}Sublead"+histname,f"{branches['Pho_isoANNScore']}Sublead"+histname,50,-0.01,1.01),"subleadIsoScore")
	hists.append(iso_scoresublead)	
	ptsublead = df_ge2pho.Histo1D((f"{branches['PhoPt']}Sublead"+histname,f"{branches['PhoPt']}Sublead"+histname,50,dist["PhoPt"]["xmin"],dist["PhoPt"]["xmax"]),"subleadPt")
	hists.append(ptsublead)



	return hists, dfs
	

def MakeBranchMapping(branches, tag):
	branch_map = {}
	for branch in branches:
		if branch == "nPhotons":
			branch_map[branch] = f"n{tag.capitalize()}Photons"
		else:
			branch_map[branch] = tag+branch
	return branch_map


#import kerebos credentials to conda env if not already there
kerb = os.getenv("KRB5CCNAME")
if(kerb is None):
	print("Setting kerebos credentials")
	os.environ["KRB5CCNAME"] = "API:"

eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47/"
eosdir_mini = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47_mini/"
files = [eosdir_mini+"SMS_SVHPM100_v33_gogoGZ_MINIAODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root",eosdir+"SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root",eosdir+"JetMET_R23_SVHPMet100_v31_JetMET1_AOD_Run2023C-19Dec2023-v1_rjrskim.root",eosdir_mini+"JetMET_R23_SVHPMet100_v31_JetMET1_MINIAOD_Run2023C-19Dec2023-v1_rjrskim.root"] 
ofilename = "photonMiniCompPlots.root"
ofile = TFile(ofilename,"RECREATE")
hists_aod = {}
hists_mini = {}
isopho = "(selPhoEta < 1.479 && selPhoEta > -1.479 && selPhoEcalRHSumEtConeDR04 < 10 && selPhoHadTowOverEM < 0.02 && selPhoTrkSumPtSolidConeDR04 < 6.)"

rjrana = RJRAnalysis()
triggers = rjrana._triggers
filters = rjrana._met_filters
presel = f"{triggers} && {filters}"

branches = ["PhoEta","nPhotons","PhoWTime","PhoWTimeSig","Pho_beamHaloCNNScore","Pho_isoANNScore","PhoPt"] 
for file in files:
	print("Doing file",file)
	isdata = False
	if "SIM" in file: #mc
		if "SMS" in file:
			histname = file[file.find("AODSIM")+6:file.rfind(".")]
			histname = histname.replace("-","_")
		else:
			histname = file[:file.find("_rjrskim")]
		isdata = False
	else: #data
		histname = parse_data_name(file[file.rfind("/"):])
		isdata = True
	ch = TChain("kuSkimTree")
	ch.Add(file)
	df0 = RDataFrame(ch)
	df_trig = df0.Filter(triggers, "triggers")
	df_filter = df0.Filter(filters,"filters")
	df_presel = df0.Filter(presel,"triggers+filters").Filter("nSelPhotons > 0","nselpho > 0")
	if isdata:
		df_presel = df_presel.Filter("selCMet < 150","met < 150")
	if "SMS" in file:	
		#do at least 1 isolated barrel photon selection for SMS files to match selection
		df = df_presel.Define("nIsoPho",f"selPhoPt[{isopho}].size()").Filter("nIsoPho > 0",">0 iso barrel pho")
	else:
		df = df_presel
	
	hists = []
	dfs = []
	branch_map = MakeBranchMapping(branches,"sel")
	hists, dfs = AddHists(df, hists, dfs, branch_map, histname)
	#isolation cuts
	eepho = "(selPhoEta < -1.479 || selPhoEta > 1.479)"
	ebpho = "(selPhoEta > -1.479 && selPhoEta < 1.479)"
	#eeisopho = f"({eepho} && selPho_isoANNScore > 0.9994431)"
	ebisopho = f"({ebpho} && {isopho})"
	eeisopho = f"({eepho} && selPho_isoANNScore > 0.003383696)"
	#ebisopho = f"({ebpho} && selPho_isoANNScore > 0.)"
	isopho = f"{eeisopho} || {ebisopho}"
	df_vvlooseisophos = df.Define("isoPhoIdxs","getIsoIdxs(selPho_isoANNScore, selPhoEta, 0.9994431, 0.003383696)").Define("nBaselinePhotons",f"selPho_isoANNScore[{isopho}].size()").Define("baselinePhoEta",f"selPhoEta[{isopho}]").Define("baselinePhoWTime",f"selPhoWTime[{isopho}]").Define("baselinePhoWTimeSig",f"selPhoWTimeSig[{isopho}]").Define("baselinePho_beamHaloCNNScore",f"selPho_beamHaloCNNScore[{isopho}]").Define("baselinePho_isoANNScore",f"selPho_isoANNScore[{isopho}]").Define("baselinePhoPt",f"selPhoPt[{isopho}]")
	df_isophos = df_vvlooseisophos.Filter("nBaselinePhotons > 0")
	branch_map = MakeBranchMapping(branches,"baseline")
	hists, dfs = AddHists(df_isophos, hists, dfs, branch_map, histname)
	if "MINI" in file:
		hists_mini[histname] = hists
	else:
		hists_aod[histname] = hists
	df0.Report().Print()

#make hist pairs
key_pairs = []
for mini_key in hists_mini.keys():
	for aod_key in hists_aod.keys():
		if aod_key == mini_key:
			#make sure order is [aod, mini] for labels
			key_pairs.append((aod_key, mini_key))
if len(key_pairs) < 1:
	print("Did not find any mass points to compare. Returning")
	exit()
ofile.cd()
for key_pair in key_pairs:
	aod_key = key_pair[0]
	mini_key = key_pair[1]
	mini_hists = hists_mini[mini_key]
	aod_hists = hists_aod[aod_key]
	#make sure order is [aod, mini] for labels
	for hists in zip(aod_hists,mini_hists):
		aod_hist = hists[0]
		mini_hist = hists[1]
		name = aod_hist.GetName()
		if "Score" in name:
			name = name.replace("_","",1)
		name = name[:name.find("_")]
		#print("doing hists for obs",name)
		sig_name = aod_key[aod_key.find("_"):]
		if "rjrskim" in sig_name:
			sig_name = sig_name[:sig_name.find("_rjrskim")]
		logy = True
		if "Eta" in name:
			logy = False
		canvas, objs = create_comparison_canvas(name, aod_hist, mini_hist, sig_name, logy=logy, labels = ["AOD","MINIAOD"])
		canvas.Write()

print("wrote file",ofilename)
ofile.Close()





