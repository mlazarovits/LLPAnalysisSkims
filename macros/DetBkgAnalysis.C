#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

void DetBkgAnalysis(){
	string file = "skims/condor_superclusters_defaultv8_beta0-1e-5_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NperGeV-0p0333333_emAlpha-1e-5_MET_R22_AL1NpSC_v31_MET.root";
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	//rh maps	
	TTreeInterface TI(file,"tree");
	vector<string> branchlist;
	branchlist.push_back("SC_trueLabel");
	branchlist.push_back("SC_predScore_physBkg");
	branchlist.push_back("SC_predScore_BH");
	branchlist.push_back("SC_predScore_spike");
	//branchlist.push_back("SC_EtaVar");
	//branchlist.push_back("SC_PhiVar");
	//branchlist.push_back("SC_TimeVar");
	//branchlist.push_back("SC_Energy");

	vector<string> subbranchlist;
	TI.SetNSubBranch("SC_nRHs_grid");
	subbranchlist.push_back("rh_iEta");
	subbranchlist.push_back("rh_iPhi");
	subbranchlist.push_back("rh_energy");

	string csvname = "SCsUnrolled.csv";
	if(gSystem->AccessPathName(file.c_str())){
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname);
	}
	
	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ');

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;
	auto truelab1_map = df.Filter("SC_trueLabel == 1").Histo2D({"physBkgMap_true","physBkgMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(truelab1_map);
	auto truelab2_map = df.Filter("SC_trueLabel == 2").Histo2D({"beamHaloMap_true","beamHaloMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(truelab2_map);
	auto truelab3_map = df.Filter("SC_trueLabel == 3").Histo2D({"spikeMap_true","spikeMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(truelab3_map);


	string physbkg_thresh = "0.67";	
	auto predlab1_map = df.Filter("SC_predScore_physBkg > "+physbkg_thresh).Histo2D({"physBkgMap_pred","physBkgMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(predlab1_map);
	
	string bh_thresh = "0.87";	
	auto predlab2_map = df.Filter("SC_predScore_BH > "+bh_thresh).Histo2D({"beamHaloMap_pred","beamHaloMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(predlab2_map);
	
	string spike_thresh = "0.65";	
	auto predlab3_map = df.Filter("SC_predScore_spike > "+spike_thresh).Histo2D({"spikeMap_pred","spikeMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
	hists2d.push_back(predlab3_map);



	//do PU subclusters
	//TI = TTreeInterface(file,"tree");
	//branchlist.clear();
	//branchlist.push_back("SC_EtaVar");
	//branchlist.push_back("SC_PhiVar");
	//branchlist.push_back("SC_TimeVar");
	//branchlist.push_back("SC_Energy");



	string ofilename = "test.root";
	cout << "Writing output to " << ofilename << endl;
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
	ofile->Close();

}
