#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

string GetOutputName(string file){
	string csvname_base = file;
	csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	return csvname_base;
}


ROOT::RDataFrame MakeRDataFrame(string file, vector<string> branchlist, string nsubobjbranch, vector<string> subbranchlist, vector<string> evtbranchlist, bool recreatecsv){
	//rh maps	
	TTreeInterface TI(file,"tree");

	TI.SetNSubBranch(nsubobjbranch);
	string csvname_base = GetOutputName(file);
	string csvname = "csv/"+csvname_base+"_PhotonsUnrolledSubcls.csv";
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
	
	return df;
}


void FillPUplots(ROOT::RDataFrame& df, vector<ROOT::RDF::RResultPtr<TH1D>>& hists1d, vector<ROOT::RDF::RResultPtr<TH2D>>& hists2d, string hextra = ""){
	string relEtaVar = "Photon_subcluster_EtaVar_prePUcleaning / Photon_EtaVar_prePUcleaning";
	string relPhiVar = "Photon_subcluster_PhiVar_prePUcleaning / Photon_PhiVar_prePUcleaning";
	string relTimeVar = "Photon_subcluster_TimeVar_prePUcleaning / Photon_TimeVar_prePUcleaning";
	string relGeoAvgFuncStr = "pow( "+relEtaVar+" * "+relPhiVar+" * "+relTimeVar+", 1./3.)";
	auto df1 = df.Define("subclRelEnergy","Photon_subcluster_Energy_prePUcleaning / Photon_Energy_prePUcleaning"); //relative energy
	auto df2 = df1.Define("subclRelGeoAvgVar",relGeoAvgFuncStr); //relative geo avg var
	auto filtered_df2 = df2.Filter("subclRelEnergy != 1 && subclRelGeoAvgVar != 1","1subcl_cut");
	
	
	auto filtered_df2p5 = filtered_df2.Filter("subclRelEnergy < 1. && subclRelGeoAvgVar < 10.","histbounds_cut");
	
	string pucut = "subclRelGeoAvgVar - subclRelEnergy <= 0";
	auto filtered_df3 = filtered_df2p5.Filter(pucut,"pu_cleaning_cut");

	string histname = "Photon_subclRelGeoAvgVar_subclRelEnergy";
	if(hextra != "")
		histname += "_"+hextra;

	auto hist2d_pu_cleaning = filtered_df2.Histo2D({histname.c_str(),(histname+";subclRelGeoAvgVar;subclRelEnergy;a.u.").c_str(),100,0,10.2,100,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
		//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
	hists2d.push_back(hist2d_pu_cleaning);


}

//void PhotonIDAnalysis(string file, bool recreatecsv = true, string ofilename_extra = ""){
void PhotonIDAnalysis(bool recreatecsv = true, string ofilename_extra = ""){


	vector<string> branchlist;
	branchlist.push_back("Photon_trueLabel");
	branchlist.push_back("Photon_predScore_isoBkg");
	branchlist.push_back("Photon_predScore_nonIsoBkg");
	branchlist.push_back("Photon_TimeCenter");
	branchlist.push_back("Photon_EtaCenter");
	branchlist.push_back("Photon_PhiCenter");
	branchlist.push_back("Photon_EtaVar");
	branchlist.push_back("Photon_PhiVar");
	branchlist.push_back("Photon_TimeVar");
	branchlist.push_back("Photon_EtaPhiCov");
	branchlist.push_back("Photon_majorLength");
	branchlist.push_back("Photon_minorLength");
	branchlist.push_back("Photon_Energy");
	branchlist.push_back("Photon_Pt");
	branchlist.push_back("Photon_nSubclusters");
	branchlist.push_back("Photon_EtaVar_prePUcleaning");
	branchlist.push_back("Photon_PhiVar_prePUcleaning");
	branchlist.push_back("Photon_TimeVar_prePUcleaning");
	branchlist.push_back("Photon_Energy_prePUcleaning");


	string nsubobjbr = "Photon_nSubclusters_prePUcleaning";
	vector<string> subbranchlist;
	subbranchlist.push_back("Photon_subcluster_EtaVar_prePUcleaning");
	subbranchlist.push_back("Photon_subcluster_PhiVar_prePUcleaning");
	subbranchlist.push_back("Photon_subcluster_TimeVar_prePUcleaning");
	subbranchlist.push_back("Photon_subcluster_Energy_prePUcleaning");


	string file_gjets = "skims/condor_photons_defaultv4p1_noIso_GJetsCR_EGamma_R18_InvMetPho30_NoSV_v31_EGamma.root";
	if(gSystem->AccessPathName(file_gjets.c_str())){
		cout << "File " << file_gjets << " does not exist." << endl;
		return;
	}
	string file_jetht = "skims/condor_photons_defaultv4p1_noIso_DijetsCR_JetHT_R18_InvMET100_nolumimask_v31_JetHT.root";	
	if(gSystem->AccessPathName(file_jetht.c_str())){
		cout << "File " << file_jetht << " does not exist." << endl;
		return;
	}

	vector<string> evtbranchlist;
	ROOT::RDataFrame df_gjets = MakeRDataFrame(file_gjets, branchlist, nsubobjbr, subbranchlist, evtbranchlist, recreatecsv);
	ROOT::RDataFrame df_jetht = MakeRDataFrame(file_jetht, branchlist, nsubobjbr, subbranchlist, evtbranchlist, recreatecsv);

	vector<ROOT::RDataFrame> dfs;
	vector<string> names;
	dfs.push_back(df_gjets);
	names.push_back("gjets");
	dfs.push_back(df_jetht);
	names.push_back("jetht");
	
	string ofilename = "GJets_JetHT_photonID";
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += ".root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();

	for(int i = 0; i < names.size(); i++){	
		vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
		vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

		auto df = dfs[i];
		string hextra = names[i];

		//do event level plots
		string file, evtsel;
		if(names[i] == "gjets"){
			file = file_gjets;
			evtsel = "PassGJetsCR";
		}
		else if(names[i] == "jetht"){
			file = file_jetht;
			evtsel = "PassDijetsCR";
		}
		else continue;
		auto df_evt = ROOT::RDataFrame("tree",file);
		auto njets = df_evt.Filter(evtsel,evtsel).Histo1D({("nSelJets_"+hextra).c_str(),("nSelJets_"+hextra).c_str(),10,0,10},"nSelJets");
		hists1d.push_back(njets);


		//FillPUplots(df, hists1d, hists2d);
		
		string label = "";
		if(hextra == "gjets")
			label = "4";
		else if(hextra == "jetht")
			label = "6";
		else continue;

		//filter out endcap photons
		auto df_noEndcap = df.Filter("Photon_EtaCenter < 1.5 && Photon_EtaCenter > -1.5","filterBarrelOnly"); 
		auto dfsig = df_noEndcap.Define("Photon_EtaSig","sqrt(Photon_EtaVar)").Define("Photon_PhiSig","sqrt(Photon_PhiVar)");
		auto df_label = dfsig.Filter("Photon_trueLabel == "+label,"filterLabel == "+label);

		auto pho_energy = df_label.Histo1D({("Photon_Energy_"+hextra).c_str(),("Photon_Energy_"+hextra).c_str(),50,0,650},"Photon_Energy");
		hists1d.push_back(pho_energy);	
		

		auto pho_etaSig = df_label.Histo1D({("Photon_EtaSig_"+hextra).c_str(),("Photon_EtaSig_"+hextra).c_str(),50,0,0.04},"Photon_EtaSig");
		hists1d.push_back(pho_etaSig);	
		auto pho_phiSig = df_label.Histo1D({("Photon_PhiSig_"+hextra).c_str(),("Photon_PhiSig_"+hextra).c_str(),50,0,0.05},"Photon_PhiSig");
		hists1d.push_back(pho_phiSig);


			
		for(auto hist : hists1d) hist->Write();
		for(auto hist : hists2d) hist->Write();
		cout << "Reporting for df " << names[i] << endl;
		df.Report()->Print();
	}		
	ofile->Close();
	cout << "Wrote output to " << ofilename << endl;

}
