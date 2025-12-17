#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

void ShapeRatioStudy(string ofilename_extra = ""){
	//bool rhmaps = true;
	string physbkgLabel = "SC_trueLabel_CMS == 1";
	string bhLabel = "SC_trueLabel_CMS == 2";


	//string csvname_base = file;//"SCsUnrolled.csv";
	//csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	//csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	string csvname_base = "condor_superclusters_defaultv9p2_MET_R18_AL1NpSC_v31_MET_AOD_Run2018C-15Feb2022_UL2018-v1";
	string csvname = "csv/"+csvname_base+"_CMS_SCsUnrolled.csv";
	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ');

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;


	auto df0 = df.Define("shapeRatio","sqrt(SC_EtaVar_CMS) / sqrt(SC_PhiVar_CMS)");
	auto df_trueBH = df0.Filter(bhLabel,"trueBH");
	auto df_truePhysBkg = df0.Filter(physbkgLabel,"truePhysBkg");
	

	auto h_truebh_shaperatio = df_trueBH.
				Histo1D({"trueBH_shapeRatio","trueBH_shapeRatio",50,0,5}, "shapeRatio");
	hists1d.push_back(h_truebh_shaperatio);	
	auto h_truephysbkg_shaperatio = df_truePhysBkg.
				Histo1D({"truePhysBkg_shapeRatio","truePhysBkg_shapeRatio",50,0,5}, "shapeRatio");
	hists1d.push_back(h_truephysbkg_shaperatio);	


	string ofilename = csvname_base+"_shapeRatio";
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
