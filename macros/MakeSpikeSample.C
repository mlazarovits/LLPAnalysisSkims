#include <iostream>
#include <vector>
#include "TTree.h"


void MakeSpikeSample(){
	string file = "skims/condor_superclusters_defaultv8p2_beta0-1e-5_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NperGeV-0p0333333_emAlpha-1e-5_MET_R18_AL1NpSC_DEOnly_v31_MET.root";


	TFile* f = new TFile(file.c_str(),"READ");
	TTree* tree = (TTree*)f->Get("tree");

	double evt;
	tree->SetBranchAddress("evt",&evt);	

	vector<vector<double>>* rh_ieta = nullptr;
	vector<vector<double>>* rh_iphi = nullptr;
	vector<vector<double>>* rh_e = nullptr;
	vector<double>* nrhs_grid = nullptr;
	tree->SetBranchAddress("SC_rh_iEta",&rh_ieta);
	tree->SetBranchAddress("SC_rh_iPhi",&rh_iphi);
	tree->SetBranchAddress("SC_rh_energy",&rh_e);
	tree->SetBranchAddress("SC_nRHs_grid",&nrhs_grid);
	vector<double>* dr_track = nullptr;
	vector<double>* phicenter = nullptr;
	vector<double>* seedtime = nullptr;
	vector<double>* timesig = nullptr;
	vector<double>* isopresel = nullptr;
	tree->SetBranchAddress("SC_dR_track",&dr_track);
	tree->SetBranchAddress("SC_PhiCenter",&phicenter);
	tree->SetBranchAddress("SC_seedTime",&seedtime);
	tree->SetBranchAddress("SC_seedTimeSignificance",&timesig);
	tree->SetBranchAddress("SC_isoPresel", &isopresel);


	string csvname = file;
	csvname = csvname.substr(csvname.find("/")+1);
	csvname = csvname.substr(0,csvname.find(".root"));
	csvname = "csv/"+csvname+"_spikeOnly.csv";

	std::ofstream ocsv;
  	ocsv.open(csvname);
	
	//write header
	ocsv << "sample,event,event_weight,object,CNNgrid_cell-3_-3,CNNgrid_cell-3_-2,CNNgrid_cell-3_-1,CNNgrid_cell-3_0,CNNgrid_cell-3_1,CNNgrid_cell-3_2,CNNgrid_cell-3_3,CNNgrid_cell-2_-3,CNNgrid_cell-2_-2,CNNgrid_cell-2_-1,CNNgrid_cell-2_0,CNNgrid_cell-2_1,CNNgrid_cell-2_2,CNNgrid_cell-2_3,CNNgrid_cell-1_-3,CNNgrid_cell-1_-2,CNNgrid_cell-1_-1,CNNgrid_cell-1_0,CNNgrid_cell-1_1,CNNgrid_cell-1_2,CNNgrid_cell-1_3,CNNgrid_cell0_-3,CNNgrid_cell0_-2,CNNgrid_cell0_-1,CNNgrid_cell0_0,CNNgrid_cell0_1,CNNgrid_cell0_2,CNNgrid_cell0_3,CNNgrid_cell1_-3,CNNgrid_cell1_-2,CNNgrid_cell1_-1,CNNgrid_cell1_0,CNNgrid_cell1_1,CNNgrid_cell1_2,CNNgrid_cell1_3,CNNgrid_cell2_-3,CNNgrid_cell2_-2,CNNgrid_cell2_-1,CNNgrid_cell2_0,CNNgrid_cell2_1,CNNgrid_cell2_2,CNNgrid_cell2_3,CNNgrid_cell3_-3,CNNgrid_cell3_-2,CNNgrid_cell3_-1,CNNgrid_cell3_0,CNNgrid_cell3_1,CNNgrid_cell3_2,CNNgrid_cell3_3,label" << endl;

	string sample = "METPD";

	int nEntries = tree->GetEntries();
	int ngrid = 7;
	for(Long64_t i=0; i< nEntries; i++){
	cout << "Event " << i << " of " << nEntries << "\r";
	cout.flush();
		tree->GetEntry(i);

		int nSCs = nrhs_grid->size();
		//apply spike selection
		if(nSCs < 1) continue;

		vector<double> rh_es(ngrid*ngrid);
		for(int n = 0; n < nSCs; n++){
			//skip if not spike
			if(isopresel->at(n) != 1) continue;
			if(dr_track->at(n) > 0.02) continue;
			if(seedtime->at(n) > -10 || seedtime->at(n) == -999) continue;
			if(timesig->at(n) >= -3) continue;
			double pc = phicenter->at(n);
			if((pc < 0.3 || (acos(-1) - 0.3 < pc && pc < acos(-1) + 0.3) || 2*acos(-1) - 0.3 < pc )) continue;
			int nrhs = nrhs_grid->at(n);
			for(int rh = 0; rh < nrhs; rh++){
				double ieta = rh_ieta->at(n).at(rh);
				double iphi = rh_iphi->at(n).at(rh);
				double e = rh_e->at(n).at(rh);
				//fill rh_es vector
				//need to convert (ieta, iphi) to flat, column-major index
				int index = (ieta + 3) * 7 + (iphi + 3);
				rh_es[index] = e;
			}
			//TODO - check that these are written out correctly by cross-checking with root file	
			//write out all flat responses
			ocsv << sample << "," << evt << ",1,"; 
			ocsv << n << ",";
			for(int rh = 0; rh < rh_es.size(); rh++){
				ocsv << rh_es[rh] << ","; 
			}
			//write out label
			ocsv << "3" << endl;	

		}


	}
	ocsv.close();
	cout << "Writing spike CSV file to " << csvname << endl;
}
