#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

double extract_after(const std::string& s, const std::string& key) {
    auto pos = s.find(key);
    if (pos == std::string::npos) return 0.0;

    pos += key.size();
    return std::stod(s.substr(pos));
}


map<string, vector<string>> procFiles;

void MakeFileMap(){
	string path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v46/";

	procFiles["SMS"] = {};
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-1500_mN2-500_mN1-100_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1000_mN1-1_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1000_mN1-250_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1000_mN1-500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1500_mN1-1000_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1500_mN1-1_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1500_mN1-250_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1500_mN1-500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1900_mN1-1000_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1900_mN1-1500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1900_mN1-1_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1900_mN1-250_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1900_mN1-500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-1000_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-1500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-1900_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-1_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-250_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-1950_mN1-500_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-500_mN1-1_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2000_mN2-500_mN1-250_ct0p1_rjrskim_v46.root");
	procFiles.at("SMS").push_back(path+"SMS_SVIPM100_v31_gogoG_AODSIM_mGl-2500_mN2-1500_mN1-1000_ct0p1_rjrskim_v46.root");

	procFiles["MET"] = {};
	//procFiles.at("MET").push_back(path+"MET_R17_SVIPM100_v31_MET_AOD_Run2017A_rjrskim_v46.root");
	procFiles.at("MET").push_back(path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018B_rjrskim_v46.root");
	procFiles.at("MET").push_back(path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018C_rjrskim_v46.root");
	procFiles.at("MET").push_back(path+"MET_R18_SVIPM100_v31_MET_AOD_Run2018D_rjrskim_v46.root");
}




//pass skimmed skim file
void RJRAnalysis(string proc, string ofilename_extra = ""){
	MakeFileMap();

	TChain chain = TChain("kuSkimTree");
	int nfiles = procFiles.at(proc).size();
	for(int f = 0; f < nfiles; f++)
		chain.Add(procFiles.at(proc)[f].c_str());
	

	vector<string> branches = {"rjr_Rs","rjr_Ms","selCMet","rjrPTS","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight","Trigger_PFMET120_PFMHT120_IDTight_PFHT60","Trigger_PFMET120_PFMHT120_IDTight","Flag_BadChargedCandidateFilter","Flag_BadPFMuonDzFilter","Flag_BadPFMuonFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_ecalBadCalibFilter","Flag_eeBadScFilter","Flag_goodVertices","Flag_hfNoisyHitsFilter","nSelPhotons","selPhoWTime","selPho_beamHaloCNNScore","selPho_physBkgCNNScore","selPhoWTimeSig"};
	ROOT::RDataFrame df(chain, branches);
	auto df0 = df.Define("rjr_Rs0","rjr_Rs[0]").Define("rjr_Ms0","rjr_Ms[0]");


	//cleaning cuts and preselection
	//MET > 150
	string metcut = "(selCMet > 150)";
	//pTs < 150 - always use index 0 for rjr variables
	string ptscut = "(rjrPTS[0] < 150)";
	//triggers applied - 2018
	string triggers = "Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || Trigger_PFMET120_PFMHT120_IDTight_PFHT60 == 1 || Trigger_PFMET120_PFMHT120_IDTight == 1";
	//event cleaning filters applied
	string met_filters = "Flag_BadChargedCandidateFilter == 1 && Flag_BadPFMuonDzFilter == 1 && Flag_BadPFMuonFilter == 1 && Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_ecalBadCalibFilter == 1 && Flag_eeBadScFilter == 1 && Flag_goodVertices && Flag_hfNoisyHitsFilter == 1";
	//hem veto for 2018 (if hemVeto then remove)
	//string hemveto = "Flag_hemVeto == 0";
	//nselphotons selection
	string nselpho = "nSelPhotons >= 1";
	//regions need to be exclusive with each other
	//beam halo region: >= 1 beam halo photon
	string bh_thresh = "0.917252";
	string pb_thresh = "0.81476355";
	string early_sel = "selPhoWTime < 0";
	string late_sel =  "selPhoWTime > 0";
	string barrelSel = "selPhoEta < 1.5 && selPhoEta > -1.5";
	
	string bh_sel = "selPho_beamHaloCNNScore > "+bh_thresh;
	string bh_sel_early = bh_sel+"&&"+early_sel;
	string bh_sel_late = bh_sel+"&&"+late_sel;

	//phys bkg region: >= 1 phys bkg photon && == 0 beam halo photons
	string pb_sel = "selPho_physBkgCNNScore > "+pb_thresh;
	string pb_sel_early = pb_sel+ "&&"+early_sel;

	

	//# photons that pass BH discr cut and are early and in the barrel
	string nPho_bhCR_early = "selPhoWTime["+bh_sel_early+"].size()";
	//# photons that pass BH discr cut and are late and in the barrel
	string nPho_bhCR_late = "selPhoWTime["+bh_sel_late+"].size()";
	string nPho_bhCR = "selPhoWTime["+bh_sel+"].size()";


	string nPho_pbCR_early = "selPhoWTime["+pb_sel_early+"].size()";
	string nPho_pbCR = "selPhoWTime["+pb_sel+"].size()";
	auto df1 = df0.Define("nPho_pbCR",nPho_pbCR).Define("nPho_bhCR",nPho_bhCR).
				Define("nPho_pbCR_early",nPho_pbCR_early).Define("nPho_bhCR_early",nPho_bhCR_early).
				Define("nPho_bhCR_late",nPho_bhCR_late);

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

	//update presel to include hem veto with new fix
	auto df_presel = df1.Filter(metcut + " && " + ptscut + " && " + triggers + " && "+met_filters,"baseline");
	auto df_ge1pho = df_presel.Filter(nselpho,"ge1nSelPho");
	auto bh_discr_score = df_ge1pho.Histo1D({"bh_discr_score","bh_discr_score",50,-0.01,1.01},"selPho_beamHaloCNNScore");
	hists1d.push_back(bh_discr_score);
	auto pb_discr_score = df_ge1pho.Histo1D({"pb_discr_score","pb_discr_score",50,-0.01,1.01},"selPho_physBkgCNNScore");
	hists1d.push_back(pb_discr_score);

	//there are photons in the middle region that are not classified as beam halo OR phys bkg that can sneak in
	//at least 1 photon in BH CR and no photons in phys bkg CR
	auto df_bhCR = df_ge1pho.Filter("nPho_bhCR >= 1 && nPho_pbCR == 0","bhCR");
	auto df_bhCR_rsCut = df_bhCR.Filter("rjr_Rs0 > 0.15");
	auto df_bhCR_msCut = df_bhCR.Filter("rjr_Ms0 > 1000");
	auto df_bhCR_msCut_rsCut = df_bhCR.Filter("rjr_Ms0 > 1000 && rjr_Rs0 > 0.15");

	string bh_cr_early_sel = "nPho_bhCR_early >= 1 && nPho_pbCR == 0";
	string bh_cr_late_sel = "nPho_bhCR_late >= 1 && nPho_pbCR == 0";
	//at least 1 photon in BH CR and early and no photons in phys bkg CR
	auto df_bhCR_early = df_bhCR.Filter(bh_cr_early_sel,"beamHaloCR_early");
	//at least 1 photon in BH CR and late and no photons in phys bkg CR
	auto df_bhCR_late = df_bhCR.Filter(bh_cr_late_sel,"beamHaloCR_late");

	string pbCR_cutname;	
	if(proc == "SMS")
		pbCR_cutname = "pbCR";
	auto df_pbCR = df_ge1pho.Filter("nPho_pbCR >= 1 && nPho_bhCR == 0",pbCR_cutname);

	vector<double> msbins = {1000,1500,2000,3000};
	vector<double> rsbins = {0.15,0.3,0.5,1};
	int n_msbins = msbins.size()-1;
	int n_rsbins = rsbins.size()-1;
	auto yields_early = df_bhCR_early.
		Histo2D({"yields_bhEarly","yields_bhEarly;Ms;Rs",n_msbins,&msbins[0],n_rsbins,&rsbins[0]},"rjr_Ms0","rjr_Rs0");
	hists2d.push_back(yields_early);
	auto yields_late = df_bhCR_late.
		Histo2D({"yields_bhLate","yields_bhLate;Ms;Rs",n_msbins,&msbins[0],n_rsbins,&rsbins[0]},"rjr_Ms0","rjr_Rs0");
	hists2d.push_back(yields_late);
	
	auto df_bhCR_early_rsCut = df_bhCR_early.Filter("rjr_Rs0 > 0.15");
	auto df_bhCR_early_msCut = df_bhCR_early.Filter("rjr_Ms0 > 1000");
	auto df_bhCR_early_msCut_rsCut = df_bhCR_early.Filter("rjr_Ms0 > 1000 && rjr_Rs0 > 0.15");
	auto df_bhCR_late_rsCut = df_bhCR_late.Filter("rjr_Rs0 > 0.15");
	auto df_bhCR_late_msCut = df_bhCR_late.Filter("rjr_Ms0 > 1000");
	auto df_bhCR_late_msCut_rsCut = df_bhCR_late.Filter("rjr_Ms0 > 1000 && rjr_Rs0 > 0.15");
	
	auto h_rjrRs_bhCR_early = df_bhCR_early_msCut.
			Histo1D({"Rs_bhCR_early_msCut","Rs_bhCR_early_msCut",50,0,1.01},"rjr_Rs0");
	hists1d.push_back(h_rjrRs_bhCR_early);
	auto h_rjrRs_bhCR_late = df_bhCR_late_msCut.
			Histo1D({"Rs_bhCR_late_msCut","Rs_bhCR_late_msCut",50,0,1.01},"rjr_Rs0");
	hists1d.push_back(h_rjrRs_bhCR_late);
	//auto h_time_bhCR_early = df_bhCR_early.
	//		Histo1D({"time_bhCR_early","time_bhCR_early",100,-20,20},"selPhoWTime");
	//hists1d.push_back(h_time_bhCR_early);
	//auto h_score_bhCR_early = df_bhCR_early.
	//		Histo1D({"score_bhCR_early","score_bhCR_early",100,0,1.02},"selPho_beamHaloCNNScore");
	//hists1d.push_back(h_score_bhCR_early);
	//auto h_score_bhCR = df_bhCR.
	//		Histo1D({"score_bhCR","score_bhCR",100,0,1.02},"selPho_beamHaloCNNScore");
	//hists1d.push_back(h_score_bhCR);
	auto h_rjrMs_bhCR_early = df_bhCR_early_rsCut.
			Histo1D({"Ms_bhCR_early_rsCut","Ms_bhCR_early_rsCut",50,0,3000},"rjr_Ms0");
	hists1d.push_back(h_rjrMs_bhCR_early);
	auto h_rjrMs_bhCR_late = df_bhCR_late_rsCut.
			Histo1D({"Ms_bhCR_late_rsCut","Ms_bhCR_late_rsCut",50,0,3000},"rjr_Ms0");
	hists1d.push_back(h_rjrMs_bhCR_late);

	auto h_rjrMs_rjrRs_bhCR_early = df_bhCR_early.
			Histo2D({"MsRs_bhCR_early","MsRs_bhCR_early;Ms [GeV];Rs;a.u.",50,0,3000,50,0,1.01},"rjr_Ms0","rjr_Rs0");
	hists2d.push_back(h_rjrMs_rjrRs_bhCR_early);
	auto h_rjrMs_rjrRs_bhCR_late = df_bhCR_late.
			Histo2D({"MsRs_bhCR_late","MsRs_bhCR_late;Ms [GeV];Rs;a.u.",50,0,3000,50,0,1.01},"rjr_Ms0","rjr_Rs0");
	hists2d.push_back(h_rjrMs_rjrRs_bhCR_late);

	//time significance
	auto h_timesig_bhCR = df_bhCR. 
			Histo1D({"timeSig_bhCR","timeSig_bhCR",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR);
	auto h_timesig_bhCR_rsCut = df_bhCR_rsCut. 
			Histo1D({"timeSig_bhCR_rsCut","timeSig_bhCR_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_rsCut);
	auto h_timesig_bhCR_msCut = df_bhCR_msCut. 
			Histo1D({"timeSig_bhCR_msCut","timeSig_bhCR_msCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_msCut);
	auto h_timesig_bhCR_msCut_rsCut = df_bhCR_msCut_rsCut. 
			Histo1D({"timeSig_bhCR_msCut_rsCut","timeSig_bhCR_msCut_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_msCut_rsCut);
	
	//time significance - early
	auto h_timesig_bhCR_early = df_bhCR_early. 
			Histo1D({"timeSig_bhCR_early","timeSig_bhCR_early",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_early);
	auto h_timesig_bhCR_early_rsCut = df_bhCR_early_rsCut. 
			Histo1D({"timeSig_bhCR_early_rsCut","timeSig_bhCR_early_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_early_rsCut);
	auto h_timesig_bhCR_early_msCut = df_bhCR_early_msCut. 
			Histo1D({"timeSig_bhCR_early_msCut","timeSig_bhCR_early_msCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_early_msCut);
	auto h_timesig_bhCR_early_msCut_rsCut = df_bhCR_early_msCut_rsCut. 
			Histo1D({"timeSig_bhCR_early_msCut_rsCut","timeSig_bhCR_early_msCut_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_early_msCut_rsCut);
	
	//time significance - late
	auto h_timesig_bhCR_late = df_bhCR_late. 
			Histo1D({"timeSig_bhCR_late","timeSig_bhCR_late",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_late);
	auto h_timesig_bhCR_late_rsCut = df_bhCR_late_rsCut. 
			Histo1D({"timeSig_bhCR_late_rsCut","timeSig_bhCR_late_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_late_rsCut);
	auto h_timesig_bhCR_late_msCut = df_bhCR_late_msCut. 
			Histo1D({"timeSig_bhCR_late_msCut","timeSig_bhCR_late_msCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_late_msCut);
	auto h_timesig_bhCR_late_msCut_rsCut = df_bhCR_late_msCut_rsCut. 
			Histo1D({"timeSig_bhCR_late_msCut_rsCut","timeSig_bhCR_late_msCut_rsCut",50,-20,20},"selPhoWTimeSig");
	hists1d.push_back(h_timesig_bhCR_late_msCut_rsCut);

	string ofilename = procFiles.at(proc)[0];
	ofilename = ofilename.substr(0,ofilename.find("_v31_"));
	ofilename = ofilename.substr(ofilename.rfind("/")+1);
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += "_rjrObsTEST.root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
	ofile->Close();
	df.Report()->Print();
	cout << "Writing output to " << ofilename << endl;



}

