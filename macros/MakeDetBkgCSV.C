#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

void MakeDetBkgCSV(string file, string tag = "CMS", bool recreatecsv = false, string ofilename_extra = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	TTreeInterface TI(file,"tree");
	vector<string> branchlist;
	branchlist.push_back("SC_trueLabel_"+tag);
	branchlist.push_back("SC_predScore_physBkg_"+tag);
	branchlist.push_back("SC_predScore_BH_"+tag);
	branchlist.push_back("SC_predScore_spike_"+tag);
	branchlist.push_back("SC_TimeCenter_"+tag);
	branchlist.push_back("SC_EtaCenter_"+tag);
	branchlist.push_back("SC_PhiCenter_"+tag);
	branchlist.push_back("SC_EtaVar_"+tag);
	branchlist.push_back("SC_PhiVar_"+tag);
	branchlist.push_back("SC_dR_track_"+tag);
	branchlist.push_back("SC_EovP_track_"+tag);
	branchlist.push_back("SC_seedTime_"+tag);
	branchlist.push_back("SC_seedTimeSignificance_"+tag);
	if(file.find("SMS_Sig") == string::npos){
		branchlist.push_back("SC_isoPresel");
		branchlist.push_back("SC_PassGJetsCR_Obj");
	}

	vector<string> subbranchlist;
	map<string, double> skipbranch;
	//skipbranch["SC_isoPresel"] = 0;
	//skipbranch["SC_trueLabel_"+tag] = -1;
	//skipbranch["SC_seedTime_"+tag] = -999;
	skipbranch["SC_trueLabel_CMS"] = -1;
	TI.SetSkipBranches(skipbranch);

	if(file.find("root://") != string::npos){ //extra parsing for eos file
		string match = "root://cmseos.fnal.gov//store/user/malazaro/LLPMVA_TrainingSamples/";
		file = file.substr(match.size());
	}

	string csvname_base = file;//"SCsUnrolled.csv";
	csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	string csvname = "csv/"+csvname_base+"_"+tag+"_SCsUnrolled.csv";
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

}
