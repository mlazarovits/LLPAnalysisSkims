#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

//TODO - decide whether this is running on KUCMSSkims or my skims...

void MakeCSVs(string& egamma_csv, string& met_csv){
	string file_met = "root://cmseos.fnal.gov//store/user/mlazarov/LLPMVA_TrainingSamples/LLPSkims_test/condor_superclusters_defaultv9p2_testDebugSelectClass_MET_R22_AL1NpSC_v31_MET_AOD_Run2022C-27Jun2023-v2.root";
	//EGamma test skim for this version does not exist....
	string file_egamma = "root://cmseos.fnal.gov//store/user/mlazarov/LLPMVA_TrainingSamples/LLPSkims_test/";

	//if using file: condor_superclusters_defaultv9p2_testDebugSelectClass_MET_R22_AL1NpSC_v31_MET_AOD_Run2022C-27Jun2023-v2.root
	//use pred score from: SC_KU_CNN_detector_withSelectClassDebugged_1000epochs_small3_CMS_predScore_physBkg/BH_CMS
	//^^this 2022 file is where the good eta-time plots are from (with that discr)!!

	vector<string> files = {file_met, file_egamma};


	string physbkgLabel = "SC_trueLabel_CMS == 1";
	string bhLabel = "SC_trueLabel_CMS == 2";
	
	vector<string> branchlist, subbranchlist, skipbranch, evtbranchlist;
	branchlist.push_back("SC_trueLabel_CMS");
	branchlist.push_back("SC_KU_CNN_detector_withSelectClassDebugged_1000epochs_small3_CMS_predScore_physBkg_CMS");
	branchlist.push_back("SC_KU_CNN_detector_withSelectClassDebugged_1000epochs_small3_CMS_predScore_BH_CMS");
	branchlist.push_back("SC_EtaVar_CMS");
	branchlist.push_back("SC_PhiVar_CMS");

	for(int f = 0; f < files.size(); f++){
		string file = files[f];
		TTreeInterface TI(file,"tree");
		
		if(file.find("root://") != string::npos){ //extra parsing for eos file
			string match = "root://cmseos.fnal.gov//store/user/malazaro/LLPMVA_TrainingSamples/";
			file = file.substr(match.size());
		}

		string csvname_base = file;//"SCsUnrolled.csv";
		csvname_base = csvname_base.substr(csvname_base.find("/")+1);
		csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
		string csvname = "csv/"+csvname_base+"_CMS_SCsUnrolled.csv";
		cout << "Creating CSV" << endl;
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname, evtbranchlist);
		if(f == 0)
			egamma_csv = csvname;
		if(f == 1)
			met_csv = csvname;
	}
	

}

//TODO - run over ku skims for latest network discriminator
void ShapeRatioStudy(string ofilename_extra = ""){

	//   vector<float>   *selPhoCovEtaEta;
   	//vector<float>   *selPhoCovEtaPhi;
   	//vector<float>   *selPhoCovPhiPhi;

	string csvname_egamma = "csv/condor_superclusters_defaultv9p1_EGamma_R18_InvMetPho30_NoSV_v31_EGamma_AOD_Run2018C_CMS_SCsUnrolled.csv";
	ROOT::RDataFrame df_egamma = ROOT::RDF::FromCSV(csvname_egamma.c_str(),true,' ');
	auto df_physbkg = df_egamma.Filter("SC_trueLabel_CMS == 1","physbkgSel");
	auto df_physbkg0 = df_physbkg.Define("shapeRatio","sqrt(SC_EtaVar_CMS) / sqrt(SC_PhiVar_CMS)");
	


	//do for MET
	string csvname_met = "csv/condor_superclusters_defaultv9p2_MET_R18_AL1NpSC_v31_MET_AOD_Run2018C-15Feb2022_UL2018-v1_CMS_SCsUnrolled.csv";
	ROOT::RDataFrame df_met = ROOT::RDF::FromCSV(csvname_met.c_str(),true,' ');
	auto df_bh = df_met.Filter("SC_trueLabel_CMS == 2","beamhaloSel");
	auto df_bh0 = df_bh.Define("shapeRatio","sqrt(SC_EtaVar_CMS) / sqrt(SC_PhiVar_CMS)");


	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;

	auto h_truebh_shaperatio = df_bh0.
				Histo1D({"trueBH_shapeRatio","trueBH_shapeRatio",50,0,8.}, "shapeRatio");
	hists1d.push_back(h_truebh_shaperatio);	
	auto h_truephysbkg_shaperatio = df_physbkg0.
				Histo1D({"truePhysBkg_shapeRatio","truePhysBkg_shapeRatio",50,0,8.}, "shapeRatio");
	hists1d.push_back(h_truephysbkg_shaperatio);	

	auto h_truebh_bhdiscrscore = df_bh0.
				Histo1D({"trueBH_predScore_beamHaloCNN","trueBH_predScore_beamHaloCNN",50,0,1}, "SC_predScore_BH_CMS");
	hists1d.push_back(h_truebh_bhdiscrscore);
	auto h_truebh_pbdiscrscore = df_bh0.
				Histo1D({"trueBH_predScore_physBkgCNN","trueBH_predScore_physBkgCNN",50,0,1}, "SC_predScore_physBkg_CMS");
	hists1d.push_back(h_truebh_pbdiscrscore);
	auto h_truepb_bhdiscrscore = df_physbkg0.
				Histo1D({"truePhysBkg_predScore_physBkgCNN","truePhysBkg_predScore_physBkgCNN",50,0,1}, "SC_predScore_physBkg_CMS");
	hists1d.push_back(h_truepb_bhdiscrscore);
	auto h_truepb_pbdiscrscore = df_physbkg0.
				Histo1D({"truePhysBkg_predScore_beamHaloCNN","trueBH_predScore_beamHaloCNN",50,0,1}, "SC_predScore_BH_CMS");
	hists1d.push_back(h_truepb_pbdiscrscore);
		

	//make shape ratio ROC curve
	//TODO - do same for score
	double ratio_max = 8;
	double step = 0.1;
	vector<std::unique_ptr<ROOT::RDF::RNode>> bh_cuts, pb_cuts, bh_discrcuts, pb_discrcuts;
	vector<double> bheffs, pbeffs, bheffs_discr, pbeffs_discr;
	for(double i = step; i < ratio_max; i += step){
		//shape ratio cuts
		auto df_bhEff = std::make_unique<ROOT::RDF::RNode>( df_bh0.Filter("shapeRatio > "+std::to_string(i),"shapeRatioCut"+to_string(i)) );
		bh_cuts.push_back(std::move(df_bhEff));
		auto df_pbEff = std::make_unique<ROOT::RDF::RNode>( df_physbkg0.Filter("shapeRatio > "+std::to_string(i),"shapeRatioCut"+to_string(i)) );
		pb_cuts.push_back(std::move(df_pbEff));
	}
	for(double i = step; i < 1; i += step){ 
		//discriminator score cuts
		auto df_bhEff_discr = std::make_unique<ROOT::RDF::RNode>( df_bh0.Filter("SC_predScore_BH_CMS > "+std::to_string(i),"predScoreCut"+to_string(i)) );
		bh_discrcuts.push_back(std::move(df_bhEff_discr));
		auto df_pbEff_discr = std::make_unique<ROOT::RDF::RNode>( df_physbkg0.Filter("SC_predScore_BH_CMS > "+std::to_string(i),"predScoreCut"+to_string(i)) );
		pb_discrcuts.push_back(std::move(df_pbEff_discr));
	}

//get cut efficiencies
	auto allCutsReport = df_met.Report();
cout << "Beam halo" << endl;
	for(auto &&cutInfo : allCutsReport){
		if(cutInfo.GetName().find("shapeRatioCut") != string::npos){
			double bhEff = (double)cutInfo.GetEff()/100.;
			bheffs.push_back(bhEff);
		}
		else if(cutInfo.GetName().find("predScoreCut") != string::npos){
			double bhEff = (double)cutInfo.GetEff()/100.;
			bheffs_discr.push_back(bhEff);
		}
		else continue;
cout << " cut name " << cutInfo.GetName() << " efficiency " << (double)cutInfo.GetEff()/100. << " all " << cutInfo.GetAll() << " pass " << cutInfo.GetPass() << endl;
	}
	allCutsReport = df_egamma.Report();
cout << "Physics Bkg" << endl;
	for(auto &&cutInfo : allCutsReport){
		if(cutInfo.GetName().find("shapeRatioCut") != string::npos){
			double pbEff = (double)cutInfo.GetEff()/100.;
			pbeffs.push_back(pbEff);
		}
		else if(cutInfo.GetName().find("predScoreCut") != string::npos){
			double pbEff = (double)cutInfo.GetEff()/100.;
			pbeffs_discr.push_back(pbEff);
		}
		else continue;
cout << " cut name " << cutInfo.GetName() << " efficiency " << (double)cutInfo.GetEff()/100. << " all " << cutInfo.GetAll() << " pass " << cutInfo.GetPass() << endl;
	}
	
	TGraph* gr = new TGraph(bheffs.size(),&pbeffs[0],&bheffs[0]);
	gr->SetTitle("Shape Ratio ROC curve;Physics Background Efficiency;Beam Halo Efficiency");
	gr->SetName("shaperatio_roc_curve");
	
	TGraph* gr_discr = new TGraph(bheffs_discr.size(),&pbeffs_discr[0],&bheffs_discr[0]);
	gr_discr->SetTitle("Discriminator ROC curve;Physics Background Efficiency;Beam Halo Efficiency");
	gr_discr->SetName("discr_roc_curve");

	string ofilename = "test_shapeRatio";
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += ".root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	cout << "Writing hists" << endl;
	for(auto hist : hists1d) hist->Write();
	gr->Write();
	gr_discr->Write();
	ofile->Close();
	df_bh0.Report()->Print();
	df_physbkg0.Report()->Print();
	cout << "Writing output to " << ofilename << endl;

}
