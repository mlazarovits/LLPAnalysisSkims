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

#import kerebos credentials to conda env if not already there
kerb = os.getenv("KRB5CCNAME")
if(kerb is None):
	print("Setting kerebos credentials")
	os.environ["KRB5CCNAME"] = "API:"

eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v47_noSMSGenFlags_noTimingAdjustments/"
files = ["SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p5_FASTSIMAOD_rjrskim.root","SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p5_rjrskim.root","SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2200_ct0p1_FASTSIMAOD_rjrskim.root","SMS_SVIPM100_v33_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2200_ct0p1_rjrskim.root"]
ofilename = "photimesFastSimComp.root"
ofile = TFile(ofilename,"RECREATE")
hists_full = {}
hists_fast = {}
for file in files:
	print("Doing file",file)
	ch = TChain("kuSkimTree")
	ch.Add(eosdir+file)
	df = RDataFrame(ch)
	masspt = file[file.find("AODSIM")+6:file.find(".")]
	masspt = masspt.replace("-","_")
	hists = []
	photime = df.Histo1D(("rawSeedTime"+masspt,"rawSeedTime"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"selPhoRawTime")
	hists.append(photime)

	timesig = df.Histo1D( ("timeSig"+masspt,"timeSig"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"selPhoWTimeSig")
	hists.append(timesig)
	
	df_leadpho = df.Define("leadRawPhoTime","selPhoRawTime[0]").Define("leadTimeSig","selPhoWTimeSig[0]")
	leadphotime = df_leadpho.Histo1D(("rawSeedTimeLead"+masspt,"rawSeedTimeLead"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"leadRawPhoTime")
	hists.append(leadphotime)
	leadtimesig = df_leadpho.Histo1D( ("timeSigLead"+masspt,"timeSigLead"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"leadTimeSig")
	hists.append(leadtimesig)


	df_ge2pho = df.Filter("nSelPhotons > 1").Define("subleadRawPhoTime","selPhoRawTime[1]").Define("subleadTimeSig","selPhoWTimeSig[1]")
	subleadphotime = df_ge2pho.Histo1D(("rawSeedTimeSublead"+masspt,"rawSeedTimeSublead"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"subleadRawPhoTime")
	hists.append(subleadphotime)
	subleadtimesig = df_ge2pho.Histo1D( ("timeSigSublead"+masspt,"timeSigSublead"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"subleadTimeSig")
	hists.append(subleadtimesig)

	eta = df.Histo1D(("selPhoEta"+masspt,"selPhoEta"+masspt,50,-3,3),"selPhoEta")
	hists.append(eta)

	#isolation cuts
	eepho = "(selPhoEta < -1.479 || selPhoEta > 1.479)"
	ebpho = "(selPhoEta > -1.479 && selPhoEta < 1.479)"
	#eeisopho = f"({eepho} && selPho_isoANNScore > 0.9994431)"
	eeisopho = f"({eepho} && selPho_isoANNScore > 0.9)"
	ebisopho = f"({ebpho} && selPho_isoANNScore > 0.003383696)"
	#ebisopho = f"({ebpho} && selPho_isoANNScore > 0.)"
	isopho = f"{eeisopho} || {ebisopho}"
	df_isophos0 = df.Define("isoPhoIdxs","getIsoIdxs(selPho_isoANNScore, selPhoEta, 0.9994431, 0.003383696)").Define("nSigPho",f"selPho_isoANNScore[{isopho}].size()").Define("sigPhoEta",f"selPhoEta[{isopho}]").Define("sigPhoRawTime",f"selPhoRawTime[{isopho}]").Define("sigPhoTimeSig",f"selPhoWTimeSig[{isopho}]")
	isophos = df_isophos0.Histo1D( ("nIsoPho"+masspt,"nIsoPho"+masspt,5,0,5),"nSigPho")
	hists.append(isophos)
	nselpho = df_isophos0.Histo1D( ("nSelPhotons"+masspt,"nSelPhotons"+masspt,5,0,5),"nSelPhotons")
	hists.append(nselpho)

	df_isophos = df_isophos0.Filter("nSigPho > 0").Define("leadSigTime","sigPhoRawTime[0]").Define("leadSigTimeSig","sigPhoTimeSig[0]")
	
	#df_isophos.Display(["selPhoEta","selPho_isoANNScore","sigPhoEta","isoPhoIdxs","nSigPho","nSelPhotons"],20).Print()
	isotime = df_isophos.Histo1D( ("rawSeedTimeIso"+masspt,"rawSeedTimeIso"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"sigPhoRawTime")
	hists.append(isotime)
	isotimesig = df_isophos.Histo1D( ("timeSigIso"+masspt,"timeSigIso"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"sigPhoTimeSig")
	hists.append(isotimesig)

	leadisotime = df_isophos.Histo1D( ("rawSeedTimeLeadIso"+masspt,"rawSeedTimeLeadIso"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"leadSigTime")
	hists.append(leadisotime)
	leadisotimesig = df_isophos.Histo1D( ("timeSigLeadIso"+masspt,"timeSigLeadIso"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"leadSigTimeSig")
	hists.append(leadisotimesig)

	df_ge2isophos = df_isophos.Filter("nSigPho > 1").Define("subleadSigTime","sigPhoRawTime[1]").Define("subleadSigTimeSig","sigPhoTimeSig[1]")
	subleadisotime = df_ge2isophos.Histo1D(("rawSeedTimeSubleadIso"+masspt,"rawSeedTimeSubleadIso"+masspt,50,dist["rawSeedTime"]["xmin"],dist["rawSeedTime"]["xmax"]),"subleadSigTime")
	hists.append(subleadisotime)
	subleadisotimesig = df_ge2isophos.Histo1D( ("timeSigSubleadIso"+masspt,"timeSigSubleadIso"+masspt,50,dist["timeSig"]["xmin"],dist["timeSig"]["xmax"]),"subleadSigTimeSig")
	hists.append(subleadisotimesig)


	isoeta = df_isophos.Histo1D( ("isoPhoEta"+masspt,"isoPhoEta"+masspt,50,-3,3),"sigPhoEta")
	hists.append(isoeta)
	if "FAST" in masspt:
		hists_fast[masspt] = hists
	else:
		hists_full[masspt] = hists
	df.Report().Print()
#make hist pairs
key_pairs = []
for fast_key in hists_fast.keys():
	for full_key in hists_full.keys():
		if full_key[:full_key.find("rjrskim")] in fast_key[:fast_key.find("FASTSIMAOD")]:
			key_pairs.append((full_key, fast_key))
ofile.cd()
for key_pair in key_pairs:
	full_key = key_pair[0]
	fast_key = key_pair[1]
	fast_hists = hists_fast[fast_key]
	full_hists = hists_full[full_key]
	for hists in zip(full_hists,fast_hists):
		full_hist = hists[0]
		fast_hist = hists[1]
		name = full_hist.GetName()
		name = name[:name.find("_")]
		#print("doing hists for obs",name)
		sig_name = full_key[:full_key.find("rjrskim")]
		logy = True
		if "Eta" in name:
			logy = False
		canvas, objs = create_comparison_canvas(name, full_hist, fast_hist, sig_name, logy=logy)
		canvas.Write()

print("wrote file",ofilename)
ofile.Close()





