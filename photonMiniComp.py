import os
from ROOT import RDataFrame, TChain, TFile, gInterpreter
from ComparisonPlotting import create_comparison_canvas, dist

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

#TODO - add discriminator distributions
def AddHists(df, hists, dfs, histname = ""):
	photime = df.Histo1D(("weightedTime"+histname,"weightedTime"+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),"selPhoWTime")
	hists.append(photime)

	timesig = df.Histo1D( ("timeSig"+histname,"timeSig"+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"selPhoWTimeSig")
	hists.append(timesig)
	
	df_leadpho = df.Define("leadTime","selPhoWTime[0]").Define("leadTimeSig","selPhoWTimeSig[0]")
	dfs.append(df_leadpho)
	leadphotime = df_leadpho.Histo1D(("weightedTimeLead"+histname,"weightedTimeLead"+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),"leadTime")
	hists.append(leadphotime)
	leadtimesig = df_leadpho.Histo1D( ("timeSigLead"+histname,"timeSigLead"+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"leadTimeSig")
	hists.append(leadtimesig)


	df_ge2pho = df.Filter("nSelPhotons > 1").Define("subleadTime","selPhoWTime[1]").Define("subleadTimeSig","selPhoWTimeSig[1]")
	dfs.append(df_ge2pho)
	subleadphotime = df_ge2pho.Histo1D(("weightedTimeSublead"+histname,"weightedTimeSublead"+histname,50,dist["weightedTime"]["xmin"],dist["weightedTime"]["xmax"]),"subleadTime")
	hists.append(subleadphotime)
	subleadtimesig = df_ge2pho.Histo1D( ("timeSigSublead"+histname,"timeSigSublead"+histname,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"subleadTimeSig")
	hists.append(subleadtimesig)

	eta = df.Histo1D(("selPhoEta"+histname,"selPhoEta"+histname,50,-3,3),"selPhoEta")
	hists.append(eta)
	return hists, dfs
	

#import kerebos credentials to conda env if not already there
kerb = os.getenv("KRB5CCNAME")
if(kerb is None):
	print("Setting kerebos credentials")
	os.environ["KRB5CCNAME"] = "API:"

eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47/"
eosdir_mini = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47_mini/"
files = [eosdir_mini+"SMS_SVHPM100_v33_gogoGZ_MINIAODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root",eosdir+"SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root"] 
ofilename = "photonMiniCompPlots.root"
ofile = TFile(ofilename,"RECREATE")
hists_aod = {}
hists_mini = {}
isopho = "(selPhoEta < 1.479 && selPhoEta > -1.479 && selPhoEcalRHSumEtConeDR04 < 10 && selPhoHadTowOverEM < 0.02 && selPhoTrkSumPtSolidConeDR04 < 6.)"
for file in files:
	print("Doing file",file)
	masspt = file[file.find("AODSIM")+6:file.rfind(".")]
	masspt = masspt.replace("-","_")
	ch = TChain("kuSkimTree")
	ch.Add(file)
	df = RDataFrame(ch)
	#do at least 1 isolated photon selection
	df = df.Define("nIsoPho",f"selPhoPt[{isopho}].size()").Filter("nIsoPho > 1")

	hists = []
	dfs = []
	hists, dfs = AddHists(df, hists, dfs, masspt)
	#isolation cuts
	eepho = "(selPhoEta < -1.479 || selPhoEta > 1.479)"
	ebpho = "(selPhoEta > -1.479 && selPhoEta < 1.479)"
	#eeisopho = f"({eepho} && selPho_isoANNScore > 0.9994431)"
	eeisopho = f"({eepho} && selPho_isoANNScore > 0.9)"
	ebisopho = f"({ebpho} && selPho_isoANNScore > 0.003383696)"
	#ebisopho = f"({ebpho} && selPho_isoANNScore > 0.)"
	isopho = f"{eeisopho} || {ebisopho}"
	df_isophos0 = df.Define("isoPhoIdxs","getIsoIdxs(selPho_isoANNScore, selPhoEta, 0.9994431, 0.003383696)").Define("nSigPho",f"selPho_isoANNScore[{isopho}].size()").Define("sigPhoEta",f"selPhoEta[{isopho}]").Define("sigPhoTime",f"selPhoWTime[{isopho}]").Define("sigPhoTimeSig",f"selPhoWTimeSig[{isopho}]")
	isophos = df_isophos0.Histo1D( ("nIsoPho"+masspt,"nIsoPho"+masspt,5,0,5),"nSigPho")
	hists.append(isophos)
	nselpho = df_isophos0.Histo1D( ("nSelPhotons"+masspt,"nSelPhotons"+masspt,5,0,5),"nSelPhotons")
	hists.append(nselpho)

	df_isophos = df_isophos0.Filter("nSigPho > 0").Define("leadSigTime","sigPhoTime[0]").Define("leadSigTimeSig","sigPhoTimeSig[0]")
	hists, dfs = AddHists(df_isophos, hists, dfs, "Iso"+masspt)
	if "MINI" in file:
		hists_mini[masspt] = hists
	else:
		hists_aod[masspt] = hists
	df.Report().Print()

#make hist pairs
key_pairs = []
for mini_key in hists_mini.keys():
	for aod_key in hists_aod.keys():
		if aod_key == mini_key:
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
	for hists in zip(aod_hists,mini_hists):
		aod_hist = hists[0]
		mini_hist = hists[1]
		name = aod_hist.GetName()
		name = name[:name.find("_")]
		print("doing hists for obs",name)
		sig_name = aod_key[:aod_key.find("rjrskim")]
		logy = True
		if "Eta" in name:
			logy = False
		canvas, objs = create_comparison_canvas(name, aod_hist, mini_hist, sig_name, logy=logy, labels = ["AOD","MINIAOD"])
		canvas.Write()

print("wrote file",ofilename)
ofile.Close()





