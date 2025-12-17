#include <iostream>
#include "TTreeInterface.h"

void MakeDetBkgCSVs(string file, string tag = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	if(tag != "")
		tag = "_"+tag;

	TTreeInterface TI(file,"tree");
	vector<string> branchlist;
	branchlist.push_back("SC_trueLabel"+tag);

	vector<string> subbranchlist;
	TI.SetNSubBranch("SC_nRHs_grid"+tag);
	subbranchlist.push_back("SC_rh_iEta"+tag);
	subbranchlist.push_back("SC_rh_iPhi"+tag);
	subbranchlist.push_back("SC_rh_Energy"+tag);
	subbranchlist.push_back("SC_rh_Weight"+tag);
	
	map<string, double> skipbranch;
	//skipbranch["SC_isoPresel"] = 0;
	skipbranch["SC_trueLabel"+tag] = -1;
	TI.SetSkipBranches(skipbranch);

	string csvname_base = file.substr(file.rfind("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));

	string csvname = "csv/"+csvname_base+tag+"_SCsRhCNNMaps.csv";
	vector<string> evtbranchlist;
	cout << "Creating CSV" << endl;
	TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);

}
