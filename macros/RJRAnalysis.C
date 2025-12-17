#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

//pass skimmed skim file
void RJRAnalysis(string file, string ofilename_extra = ""){

	string ofilename = file.substr(0,file.find(".root"));
	ROOT::RDataFrame df("kuSkimTree",file.c_str());
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
	
	string bh_sel = "selPho_beamHaloCNNScore > "+bh_thresh;
	string bh_sel_early = bh_sel;//+"&&"+early_sel;
	string bh_sel_late = bh_sel+"&&"+late_sel;

	//phys bkg region: >= 1 phys bkg photon && == 0 beam halo photons
	string pb_sel = "selPho_physBkgCNNScore > "+pb_thresh;
	string pb_sel_early = pb_sel+ "&&"+early_sel;

	//# photons that pass BH discr cut and are early
	string nPho_bhCR_early = "selPhoWTime["+bh_sel_early+"].size()";
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

	//there are photons in the middle region that are not classified as beam halo OR phys bkg that can sneak in
	//at least 1 photon in BH CR and no photons in phys bkg CR
	auto df_bhCR = df_ge1pho.Filter("nPho_bhCR >= 1 && nPho_pbCR == 0","bhCR");

	string bh_cr_early_sel = "nPho_bhCR_early >= 1 && nPho_pbCR == 0";
	string bh_cr_late_sel = "nPho_bhCR_late >= 1 && nPho_pbCR == 0";
	//at least 1 photon in BH CR and early and no photons in phys bkg CR
	auto df_bhCR_early = df_ge1pho.Filter(bh_cr_early_sel,"beamHaloCR_early");
	//at least 1 photon in BH CR and late and no photons in phys bkg CR
	auto df_bhCR_late = df_ge1pho.Filter(bh_cr_late_sel,"beamHaloCR_late");
	
	//auto bf_pbCR_early = df_ge1pho.Filter(pb_cr_sel,"physBkgCR_early");


	auto df_bhCR_early_rsCut = df_bhCR_early.Filter("rjr_Rs0 > 0.15");
	auto df_bhCR_early_msCut = df_bhCR_early.Filter("rjr_Ms0 > 1000");
	auto df_bhCR_late_rsCut = df_bhCR_late.Filter("rjr_Rs0 > 0.15");
	auto df_bhCR_late_msCut = df_bhCR_late.Filter("rjr_Ms0 > 1000");
	
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




	
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += "_rjrObs.root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
	ofile->Close();
	df.Report()->Print();
	cout << "Writing output to " << ofilename << endl;




}

