#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

void DetBkgAnalysis(string file, bool recreatecsv = true, string ofilename_extra = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	bool rhmaps = true;
	string physbkgLabel = "SC_trueLabel == 1";
	//string physbkgLabel = "SC_seedTime > -0.5 && SC_seedTime < 0.5 && SC_isoPresel == 1 && SC_seedTimeSignificance < 1 && SC_seedTimeSignificance > -1 && SC_dR_track > 0.03";
	string bhLabel = "SC_trueLabel == 2";
	string phicut = "0.1";
	//string bhLabel = "SC_isoPresel == 1 && SC_dR_track > 0.03 && SC_seedTime > -7 && SC_seedTime <= -2 && SC_seedTimeSignificance < -3 && (SC_PhiCenter < "+phicut+" || (acos(-1) - "+phicut+" < SC_PhiCenter && SC_PhiCenter < acos(-1) + "+phicut+") || 2*acos(-1) - "+phicut+" < SC_PhiCenter )";
	//string bhLabel_old = "SC_isoPresel == 1 && SC_dR_track > 0.03 && SC_TimeCenter > -7 && SC_seedTime <= -2 && (SC_PhiCenter < "+phicut+" || (acos(-1) - "+phicut+" < SC_PhiCenter && SC_PhiCenter < acos(-1) + "+phicut+") || 2*acos(-1) - "+phicut+" < SC_PhiCenter )";
	string bhLabel_old = bhLabel;
	string spikeLabel = "SC_trueLabel == 3";
	//string spikeLabel = "SC_isoPresel == 1 && SC_dR_track <= 0.02 && SC_seedTime <= -10 && SC_seedTime != -999 && SC_seedTimeSignificance < -3 && !(SC_PhiCenter < 0.3 || (acos(-1) - 0.3 < SC_PhiCenter && SC_PhiCenter < acos(-1) + 0.3) || 2*acos(-1) - 0.3 < SC_PhiCenter )";

	//rh maps	
	TTreeInterface TI(file,"tree");
	vector<string> branchlist;
	branchlist.push_back("SC_trueLabel");
	branchlist.push_back("SC_predScore_physBkg");
	branchlist.push_back("SC_predScore_BH");
	branchlist.push_back("SC_predScore_spike");
	branchlist.push_back("SC_TimeCenter");
	branchlist.push_back("SC_EtaCenter");
	branchlist.push_back("SC_PhiCenter");
	branchlist.push_back("SC_EtaVar");
	branchlist.push_back("SC_PhiVar");
	branchlist.push_back("SC_dR_trackSubcl");
	branchlist.push_back("SC_EovP_trackSubcl");
	//branchlist.push_back("SC_dR_track");
	//branchlist.push_back("SC_EovP_track");
	//branchlist.push_back("SC_seedTime");
	//branchlist.push_back("SC_seedTimeSignificance");
	//branchlist.push_back("SC_isoPresel");
	//branchlist.push_back("SC_PassGJetsCR_Obj");

	vector<string> subbranchlist;
	if(rhmaps){
		TI.SetNSubBranch("SC_nRHs_grid");
		subbranchlist.push_back("rh_iEta");
		subbranchlist.push_back("rh_iPhi");
		subbranchlist.push_back("rh_energy");
		//subbranchlist.push_back("SC_rh_iEta");
		//subbranchlist.push_back("SC_rh_iPhi");
		//subbranchlist.push_back("SC_rh_energy");
	}
	
	map<string, double> skipbranch;
	//skipbranch["SC_isoPresel"] = 0;
	//skipbranch["SC_trueLabel"] = -1;
	TI.SetSkipBranches(skipbranch);

	string csvname_base = file;//"SCsUnrolled.csv";
	csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	string csvname = "csv/"+csvname_base+"_SCsUnrolledRhs.csv";
	vector<string> evtbranchlist;
	evtbranchlist.push_back("Flag_globalSuperTightHalo2016Filter");
	if(gSystem->AccessPathName(csvname.c_str())){
		cout << "Creating CSV" << endl;
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);
	}
	else{
		if(recreatecsv){
			cout << "Recreating CSV" << endl;
			TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);
		}
		else cout << "Using CSV " << csvname << endl;
	}
	

	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ');

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

	string predBH = "SC_predScore_BH > 0.67";
	string earlytimes = "SC_TimeCenter < -2";
	if(rhmaps){	
		auto truelab1_map = df.Filter(physbkgLabel).Histo2D({"physBkgMap_true","physBkgMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(truelab1_map);
		auto truelab2_map = df.Filter(bhLabel).Histo2D({"beamHaloMap_true","beamHaloMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(truelab2_map);
		auto truelab3_map = df.Filter(spikeLabel).Histo2D({"spikeMap_true","spikeMap_true;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(truelab3_map);

		string physbkg_thresh = "0.67";	
		auto predlab1_map = df.Filter("SC_predScore_physBkg > "+physbkg_thresh).Histo2D({"physBkgMap_pred","physBkgMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(predlab1_map);
		
		auto predlab2_map = df.Filter(predBH).Histo2D({"beamHaloMap_pred","beamHaloMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(predlab2_map);
		
		string spike_thresh = "0.65";	
		auto predlab3_map = df.Filter("SC_predScore_spike > "+spike_thresh).Histo2D({"spikeMap_pred","spikeMap_pred;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(predlab3_map);
		
		auto predlab2early_map = df.Filter(predBH+" && "+earlytimes).Histo2D({"beamHaloMap_predEarly","beamHaloMap_predEarly;rh_ieta;rh_iphi;energy",7,-3,4,7,-3,4}, "rh_iEta","rh_iPhi","rh_energy");
		hists2d.push_back(predlab2early_map);
	}

	auto dfphi = df.Define("SC_PhiSig","sqrt(SC_PhiVar)");
	auto df0 = dfphi.Define("SC_EtaSig","sqrt(SC_EtaVar)");
	//auto df0 = dfeta.Filter("SC_trueLabel != -1","inCR");
	auto df1 = df0.Filter(physbkgLabel,"physbkgCR");
	auto df2 = df0.Filter(bhLabel,"bhCR");
	auto df3 = df0.Filter(spikeLabel,"spikeCR");
	auto df_bhOld = df0.Filter(bhLabel_old,"bhLabel_old");


	//bh filter testing
	auto h_nom_tc_ec = df0. 
				Histo2D({"SC_TimeCenter_EtaCenter","SC_TimeCenter_EtaCenter;TimeCenter;EtaCenter;a.u.",50,-16,16,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
	hists2d.push_back(h_nom_tc_ec);
	auto h_bhfilter_tc_ec = df0.Filter("Flag_globalSuperTightHalo2016Filter == 1"). 
				Histo2D({"SC_TimeCenter_EtaCenter_BHFilter","SC_TimeCenter_EtaCenter_BHFilter;TimeCenter;EtaCenter;a.u.",50,-16,16,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
	hists2d.push_back(h_bhfilter_tc_ec);





	//spikes in phys bkg
	auto h_dr_EovP_label1 = df1.
				Histo2D({"SC_dRtrack_EovPtrack_physBkgCR","SC_dRtrack_EovPtrack_physBkgCR;dRtrack;EovPtrack",50,0,0.1,50,0,5},"SC_dR_trackSubcl","SC_EovP_trackSubcl");
				//Histo2D({"SC_dRtrack_EovPtrack_physBkgCR","SC_dRtrack_EovPtrack_physBkgCR;dRtrack;EovPtrack",50,0,0.1,50,0,5},"SC_dR_track","SC_EovP_track");
	hists2d.push_back(h_dr_EovP_label1);	
	auto h_dr_EovP_label2 = df2.
				Histo2D({"SC_dRtrack_EovPtrack_beamHaloCR","SC_dRtrack_EovPtrack_beamHaloCR;dRtrack;EovPtrack",50,0,0.1,50,0,250},"SC_dR_trackSubcl","SC_EovP_trackSubcl");
				//Histo2D({"SC_dRtrack_EovPtrack_beamHaloCR","SC_dRtrack_EovPtrack_beamHaloCR;dRtrack;EovPtrack",50,0,0.1,50,0,250},"SC_dR_track","SC_EovP_track");
	hists2d.push_back(h_dr_EovP_label2);	
	auto h_dr_EovP_label3 = df3.
				Histo2D({"SC_dRtrack_EovPtrack_spikeCR","SC_dRtrack_EovPtrack_spikeCR;dRtrack;EovPtrack",50,0,0.1,50,0,5},"SC_dR_trackSubcl","SC_EovP_trackSubcl");
				//Histo2D({"SC_dRtrack_EovPtrack_spikeCR","SC_dRtrack_EovPtrack_spikeCR;dRtrack;EovPtrack",50,0,0.1,50,0,5},"SC_dR_track","SC_EovP_track");
	hists2d.push_back(h_dr_EovP_label3);	
	auto h_ec_tc_label1 = df1.
				Histo2D({"SC_TimeCenter_EtaCenter_physBkgCR","SC_TimeCenter_EtaCenter_physBkgCR;TimeCenter;EtaCenter;a.u.",50,-0.6,0.6,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
				//Histo2D({"SC_seedTime_EtaCenter_physBkgCR","SC_seedTime_EtaCenter_physBkgCR;seedTime;EtaCenter",50,-0.6,0.6,50,-1.6,1.6},"SC_seedTime","SC_EtaCenter");
	hists2d.push_back(h_ec_tc_label1);	
	auto h_ec_tc_label2 = df2.
				Histo2D({"SC_TimeCenter_EtaCenter_beamHaloCR","SC_TimeCenter_EtaCenter_beamHaloCR;TimeCenter;EtaCenter;a.u.",50,-8,-1,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
				//Histo2D({"SC_seedTime_EtaCenter_beamHaloCR","SC_seedTime_EtaCenter_beamHaloCR;seedTime;EtaCenter",50,-8,-1,50,-1.6,1.6},"SC_seedTime","SC_EtaCenter");
	hists2d.push_back(h_ec_tc_label2);	
	auto h_ec_tc_label3 = df3.
				Histo2D({"SC_TimeCenter_EtaCenter_spikeCR","SC_TimeCenter_EtaCenter_spikeCR;TimeCenter;EtaCenter;a.u.",50,-20,-7,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
				//Histo2D({"SC_seedTime_EtaCenter_spikeCR","SC_seedTime_EtaCenter_spikeCR;seedTime;EtaCenter",50,-20,-7,50,-1.6,1.6},"SC_seedTime","SC_EtaCenter");
	hists2d.push_back(h_ec_tc_label3);	

	//phys bkg in spikes
	auto h_phisig_label1 = df1.
				Histo1D({"SC_PhiSig_physBkgCR", "SC_PhiSig_physBkgCR", 50, 0,0.1},"SC_PhiSig");
	hists1d.push_back(h_phisig_label1);	
	auto h_phisig_label2 = df2.
				Histo1D({"SC_PhiSig_beamHaloCR", "SC_PhiSig_beamHaloCR", 50, 0,0.1},"SC_PhiSig");
	hists1d.push_back(h_phisig_label2);	
	auto h_phisig_label3 = df3.
				Histo1D({"SC_PhiSig_spikeCR", "SC_PhiSig_spikeCR", 50, 0,0.1},"SC_PhiSig");
	hists1d.push_back(h_phisig_label3);	

	//phys bkg in beam halo
	auto h_etasig_label2 = df2.
				Histo1D({"SC_EtaSig_beamHaloCR","SC_EtaSig_beamHaloCR",50,0,0.03},"SC_EtaSig");
	hists1d.push_back(h_etasig_label2);	
	
	//spikes in bh
	auto h_etasig_label3 = df3.
				Histo1D({"SC_EtaSig_spikeCR","SC_EtaSig_spikeCR",50,0,0.03},"SC_EtaSig");
	hists1d.push_back(h_etasig_label3);	

	//bh in phys bkg
	auto h_etasig_label1 = df1.
				Histo1D({"SC_EtaSig_physBkgCR","SC_EtaSig_physBkgCR",50,0,0.03},"SC_EtaSig");
	hists1d.push_back(h_etasig_label1);	

	//closure test
	auto df_predBH = df0.Filter(predBH,"predictedBH");
	auto h_time_eta_predBH = df_predBH. 
				Histo2D({"SC_TimeCenter_EtaCenter_predBH","SC_TimeCenter_EtaCenter_predBH;TimeCenter;EtaCenter;a.u.",50,-16,16,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
	hists2d.push_back(h_time_eta_predBH);

	string ofilename = csvname_base+"_detbkg";
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += ".root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
	ofile->Close();
	df.Report()->Print();
	cout << "Writing output to " << ofilename << endl;

}
