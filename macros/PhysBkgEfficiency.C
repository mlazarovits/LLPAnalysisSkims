#include <iostream>
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"


void PhysBkgEfficiency(string file, bool recreatecsv = true){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}


	TTreeInterface TI(file,"tree");
	vector<string> evtbranchlist;
	evtbranchlist.push_back("PassGJetsCR");
	evtbranchlist.push_back("Flag_BadChargedCandidateFilter");
	evtbranchlist.push_back("Flag_BadPFMuonDzFilter");
	evtbranchlist.push_back("Flag_BadPFMuonFilter");
	evtbranchlist.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	evtbranchlist.push_back("Flag_HBHENoiseFilter");
	evtbranchlist.push_back("Flag_HBHENoiseIsoFilter");
	evtbranchlist.push_back("Flag_ecalBadCalibFilter");
	evtbranchlist.push_back("Flag_goodVertices");
	evtbranchlist.push_back("Flag_hfNoisyHitsFilter");
	evtbranchlist.push_back("Flag_globalSuperTightHalo2016Filter");
	vector<string> branchlist;
	branchlist.push_back("SC_seedTime");
	branchlist.push_back("SC_seedTimeSignificance");
	branchlist.push_back("SC_PhiCenter");
	branchlist.push_back("SC_dR_track");
	branchlist.push_back("SC_PassGJetsCR_Obj");
	branchlist.push_back("SC_dPhi_PhoJetSys");
	branchlist.push_back("SC_JetObjPtAsym");
	branchlist.push_back("SC_Photon_Pt");
	branchlist.push_back("SC_isoPresel");
	vector<string> subbranchlist;

	map<string, double> skipbranch;
	skipbranch["SC_isoPresel"] = 0;
	TI.SetSkipBranches(skipbranch);

	
	string csvname_base = file;//"SCsUnrolled.csv";
	csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	string csvname = "csv/"+csvname_base+"_SCsUnrolledRhs.csv";
	if(gSystem->AccessPathName(csvname.c_str())){
		cout << "Creating CSV" << endl;
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);
	}
	
	if(!gSystem->AccessPathName(csvname.c_str()) && recreatecsv){
		cout << "Recreating CSV" << endl;
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);
	}
	else cout << "Using CSV " << csvname << endl;


	
	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(), true,' ');
	//do event level GJets selection
	auto df_passAllEvtFilters = df.Filter("Flag_BadChargedCandidateFilter == 1 && Flag_BadPFMuonDzFilter == 1 && Flag_BadPFMuonFilter == 1 && Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_ecalBadCalibFilter == 1 && Flag_goodVertices == 1 && Flag_hfNoisyHitsFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1","passEvtFilters");
	auto df_passGJetsEvt = df_passAllEvtFilters.Filter("PassGJetsCR == 1","passGJetsEvtCR");

	//do object level GJets selection
	auto df_passGJetsObj = df_passGJetsEvt.Filter("SC_PassGJetsCR_Obj == 1","passGJetsObjCR"); 
	auto df_passDphiGJetsObj = df_passGJetsEvt.Filter("SC_dPhi_PhoJetSys > acos(-1) - 0.5","passDphiGJetsObj");
	auto df_passPtAsymGJetsObj = df_passDphiGJetsObj.Filter("SC_JetObjPtAsym > 0.6","passPtAsymGJetsObj");
	auto df_passPtMinGJetsObj = df_passPtAsymGJetsObj.Filter("SC_Photon_Pt > 50","passPtMinGJetsObj");


	auto df_passTime = df_passPtMinGJetsObj.Filter("SC_seedTime < 0.5 && SC_seedTime > -0.5","passTimeNew");
	auto df_passIso = df_passTime.Filter("SC_isoPresel == 1","passIsoNew");
	auto df_passTimeSig = df_passIso.Filter("SC_seedTimeSignificance < 3 && SC_seedTimeSignificance > -3","passTimeSigNew");
	auto df_passPhiCenter = df_passTimeSig.Filter("!(SC_PhiCenter < 0.1 || (acos(-1) - 0.1 < SC_PhiCenter && SC_PhiCenter < acos(-1) + 0.1) || 2*acos(-1) - 0.1 < SC_PhiCenter )","passPhiCenterNew");
	auto df_passTrackMatchVeto = df_passPhiCenter.Filter("SC_dR_track > 0.03","passTrackMatchVetoNew");
	
	//do label selection
	auto df_passGJetsObj_passTime = df_passGJetsObj.Filter("SC_seedTime < 0.5 && SC_seedTime > -0.5","passTime");
	auto df_passGJetsObj_passIso = df_passGJetsObj_passTime.Filter("SC_isoPresel == 1","passIso");
	auto df_passGJetsObj_passTimeSig = df_passGJetsObj_passIso.Filter("SC_seedTimeSignificance < 1 && SC_seedTimeSignificance > -1","passTimeSig");
	auto df_passGJetsObj_passPhiCenter = df_passGJetsObj_passTimeSig.Filter("!(SC_PhiCenter < 0.1 || (acos(-1) - 0.1 < SC_PhiCenter && SC_PhiCenter < acos(-1) + 0.1) || 2*acos(-1) - 0.1 < SC_PhiCenter )","passPhiCenter");
	auto df_passGJetsObj_passTrackMatchVeto = df_passGJetsObj_passPhiCenter.Filter("SC_dR_track > 0.03","passTrackMatchVeto");

cout << "With GJetsObj new selection." << endl;
	df_passTrackMatchVeto.Report()->Print();
cout << "With GJetsObj nominal selection." << endl;
	df_passGJetsObj_passTrackMatchVeto.Report()->Print();	
}	
