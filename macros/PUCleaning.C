#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

void PUCleaning(string file, bool recreatecsv = true){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	string physbkgLabel = "SC_trueLabel == 1";
	string bhLabel = "SC_trueLabel == 2";
	//string spikeLabel = "SC_trueLabel == 3";
	string spikeLabel = "SC_isoPresel == 1 && SC_dR_track <= 0.02 && SC_seedTime <= -10 && SC_seedTime != -999 && SC_seedTimeSignificance < -3 && !(SC_PhiCenter < 0.3 || (acos(-1) - 0.3 < SC_PhiCenter && SC_PhiCenter < acos(-1) + 0.3) || 2*acos(-1) - 0.3 < SC_PhiCenter )";

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;
	//do PU subclusters
	vector<string> branchlist, subbranchlist;
	TTreeInterface TI = TTreeInterface(file,"tree");
	TI.SetNSubBranch("SC_nSubclusters_prePUcleaning");
	branchlist.push_back("SC_EtaVar_prePUcleaning");
	branchlist.push_back("SC_PhiVar_prePUcleaning");
	branchlist.push_back("SC_TimeVar_prePUcleaning");
	branchlist.push_back("SC_Energy_prePUcleaning");
	branchlist.push_back("SC_nSubclusters");
	branchlist.push_back("SC_seedTime");
	subbranchlist.push_back("SC_subcluster_EtaVar_prePUcleaning");
	subbranchlist.push_back("SC_subcluster_PhiVar_prePUcleaning");
	subbranchlist.push_back("SC_subcluster_TimeVar_prePUcleaning");
	subbranchlist.push_back("SC_subcluster_Energy_prePUcleaning");

	map<string, double> skipbranch;
	skipbranch["SC_nSubclusters_prePUcleaning"] = 1;
	skipbranch["SC_nSubclusters"] = 0;
	TI.SetSkipBranches(skipbranch);

	string csvname_base = file;//"SCsUnrolled.csv";
	csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	string csvname = "csv/"+csvname_base+"_SCsUnrolledSubclusters.csv";
	if(gSystem->AccessPathName(csvname.c_str())){
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname);
	}
	
	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ');


	string relEtaVar = "SC_subcluster_EtaVar_prePUcleaning / SC_EtaVar_prePUcleaning";
	string relPhiVar = "SC_subcluster_PhiVar_prePUcleaning / SC_PhiVar_prePUcleaning";
	string relTimeVar = "SC_subcluster_TimeVar_prePUcleaning / SC_TimeVar_prePUcleaning";
	string relGeoAvgFuncStr = "pow( "+relEtaVar+" * "+relPhiVar+" * "+relTimeVar+", 1./3.)";
	auto df4 = df.Define("subclRelEnergy","SC_subcluster_Energy_prePUcleaning / SC_Energy_prePUcleaning"); //relative energy
	auto df5 = df4.Define("subclRelGeoAvgVar",relGeoAvgFuncStr).Define("subclRelTimeVar",relTimeVar).Define("subclRelEtaVar",relEtaVar).Define("subclRelPhiVar",relPhiVar); //relative geo avg var
	//auto filtered_df5 = df5.Filter("subclRelEnergy != 1 && subclRelGeoAvgVar != 1","1subcl_cut");
	auto filtered_df5 = df5.Filter("SC_nSubclusters_prePUcleaning != 1","1subcl_cut");
	
	
	auto filtered_df5p5 = filtered_df5.Filter("subclRelEnergy < 1. && subclRelGeoAvgVar < 10.","histbounds_cut");
	

	//auto hist2d_pu_cleaning = filtered_df5.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy","SC_subclRelGeoAvgVar_subclRelEnergy;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
	auto hist2d_pu_cleaning = filtered_df5.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy","SC_subclRelGeoAvgVar_subclRelEnergy;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,1.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
		//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
	hists2d.push_back(hist2d_pu_cleaning);

	//filter filtered_df5 based on seed times corresponding to CRs
	auto df_spikeTimes = filtered_df5.Filter("SC_seedTime <= -10 && SC_seedTime != -999");
	//auto hist2d_pu_cleaning_spikeTimes = df_spikeTimes.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy_spikeTimes","SC_subclRelGeoAvgVar_subclRelEnergy_spikeTimes;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
	auto hist2d_pu_cleaning_spikeTimes = df_spikeTimes.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy_spikeTimes","SC_subclRelGeoAvgVar_subclRelEnergy_spikeTimes;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,1.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
	hists2d.push_back(hist2d_pu_cleaning_spikeTimes);

	auto df_bhTimes = filtered_df5.Filter("SC_seedTime <= -2 && SC_seedTime > -7");
	//auto hist2d_pu_cleaning_bhTimes = df_bhTimes.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy_bhTimes","SC_subclRelGeoAvgVar_subclRelEnergy_bhTimes;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
	auto hist2d_pu_cleaning_bhTimes = df_bhTimes.Histo2D({"SC_subclRelGeoAvgVar_subclRelEnergy_bhTimes","SC_subclRelGeoAvgVar_subclRelEnergy_bhTimes;subclRelGeoAvgVar;subclRelEnergy;a.u.",100,0,1.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
	hists2d.push_back(hist2d_pu_cleaning_bhTimes);


	auto timevar = filtered_df5.Histo1D({"subclRelTimeVar","subclRelTimeVar",50,0,10},"subclRelTimeVar");
	hists1d.push_back(timevar);
	auto etavar = filtered_df5.Histo1D({"subclRelEtaVar","subclRelEtaVar",50,0,10},"subclRelEtaVar");
	hists1d.push_back(etavar);
	auto phivar = filtered_df5.Histo1D({"subclRelPhiVar","subclRelPhiVar",50,0,10},"subclRelPhiVar");
	hists1d.push_back(phivar);

	auto hist2d_relTimeVar_relE = filtered_df5.Histo2D({"SC_subclRelTimeVar_subclRelEnergy","SC_subclRelTimeVar_subclRelEnergy;subclRelTimeVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelTimeVar","subclRelEnergy");
	//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
	hists2d.push_back(hist2d_relTimeVar_relE);
	
	auto hist2d_relEtaVar_relE = filtered_df5.Histo2D({"SC_subclRelEtaVar_subclRelEnergy","SC_subclRelEtaVar_subclRelEnergy;subclRelEtaVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelEtaVar","subclRelEnergy");
	//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
	hists2d.push_back(hist2d_relEtaVar_relE);
	
	auto hist2d_relPhiVar_relE = filtered_df5.Histo2D({"SC_subclRelPhiVar_subclRelEnergy","SC_subclRelPhiVar_subclRelEnergy;subclRelPhiVar;subclRelEnergy;a.u.",100,0,10.2,100,0,1.2},"subclRelPhiVar","subclRelEnergy");
	//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
	hists2d.push_back(hist2d_relPhiVar_relE);
	
	//filted filtered_df5p5 based on # of starting subclusters and ending subclusters
	auto df_2to1 = filtered_df5p5.Filter("SC_nSubclusters_prePUcleaning == 2 && SC_nSubclusters == 1","2to1");
	auto df_ge2to1 = filtered_df5p5.Filter("SC_nSubclusters_prePUcleaning > 2 && SC_nSubclusters == 1","ge2to1");
	auto df_2toge1 = filtered_df5p5.Filter("SC_nSubclusters_prePUcleaning == 2 && SC_nSubclusters > 1","2toge1");
	auto df_ge2toge1 = filtered_df5p5.Filter("SC_nSubclusters_prePUcleaning > 2 && SC_nSubclusters > 1","ge2toge1");

	string pucut = "subclRelGeoAvgVar - subclRelEnergy <= 0";
	auto filtered_df6 = filtered_df5p5.Filter(pucut,"pu_cleaning_cut");

	string ofilename = csvname_base+"_pucleaning.root";
	cout << "Writing output to " << ofilename << endl;
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
	ofile->Close();
	df.Report()->Print();

}
