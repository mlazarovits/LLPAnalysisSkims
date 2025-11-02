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

	TTreeInterface TI(file,"tree");


	vector<string> branchlist;
	branchlist.push_back("SC_trueLabel");
	branchlist.push_back("SC_predLabel");

	vector<string> subbranchlist;
	TI.SetNSubBranch("SC_nRHs");
	subbranchlist.push_back("rh_iEta");
	subbranchlist.push_back("rh_iPhi");
	subbranchlist.push_back("rh_energy");

	string csvname = "SCsUnrolled.csv";
	TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname);
	
	//ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ',-1LL,std::move(coltypes));


}
