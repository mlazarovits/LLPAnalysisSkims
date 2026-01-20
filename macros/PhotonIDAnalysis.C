#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

string GetOutputName(string file){
	string csvname_base = file;
	string eosdir = "root://cmseos.fnal.gov//store/user/mlazarov/LLPMVA_TrainingSamples/LLPSkims_test";
	if(file.find(eosdir) == string::npos){
		csvname_base = csvname_base.substr(csvname_base.find("/")+1);
	}
	else{
		csvname_base = csvname_base.substr(csvname_base.find(eosdir)+eosdir.size()+1);
	}
	csvname_base = csvname_base.substr(0,csvname_base.find(".root"));
	cout << "file " << file << " csvname_base " << csvname_base << endl;
	return csvname_base;
}


ROOT::RDataFrame MakeRDataFrame_CSV(string file, vector<string> branchlist, vector<string> subbranchlist, vector<string> evtbranchlist, bool recreatecsv){
	TTreeInterface TI(file,"tree");
	string csvname_base = GetOutputName(file);
	string csvname = "csv/"+csvname_base+"_PhotonsUnrolled.csv";
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


map<string, vector<string>> file_map;
void MakeFileMaps(){
	string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v46/";
	//data files

	//sms files
	file_map["SMS_GlGl"] = {};
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1300_mN1-1000_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-500_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2150_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2150_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2200_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1000_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2400_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1500_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2300_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2350_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p5"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2350_ct0p5_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1300_mN1-1000_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2200_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1000_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2100_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2000_mN1-1500_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2250_mN1-2150_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2350_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-2200_mN1-2150_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2300_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2400_mN1-2350_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-1000_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2300_mN2-1600_mN1-500_ct0p1_rjrskim.root");
	file_map["SMS_GlGl_0p1"].push_back(eosdir+"SMS_SVIPM100_v31_gogoGZ_AODSIM_mGl-2500_mN2-2450_mN1-2400_ct0p1_rjrskim.root");

	//mc bkg files
	file_map["GJets"] = {};
	file_map["GJets"].push_back(eosdir+"GJets_R18_SVIPM100_v31_GJets_HT-40To100_rjrskim.root");
	file_map["GJets"].push_back(eosdir+"GJets_R18_SVIPM100_v31_GJets_HT-100To200_rjrskim.root");
	file_map["GJets"].push_back(eosdir+"GJets_R18_SVIPM100_v31_GJets_HT-200To400_rjrskim.root");
	file_map["GJets"].push_back(eosdir+"GJets_R18_SVIPM100_v31_GJets_HT-400To600_rjrskim.root");
	file_map["GJets"].push_back(eosdir+"GJets_R18_SVIPM100_v31_GJets_HT-600ToInf_rjrskim.root");

	file_map["QCD"] = {};
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT1000to1500_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT100to200_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT1500to2000_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT2000toInf_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT200to300_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT300to500_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT500to700_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT50to100_rjrskim.root");
	file_map["QCD"].push_back(eosdir+"QCD_R18_SVIPM100_v31_QCD_HT700to1000_rjrskim.root");

}


//void PhotonIDAnalysis(string file, bool recreatecsv = true, string ofilename_extra = ""){
void PhotonIDAnalysis(bool recreatecsv = false, string ofilename_extra = ""){
	//cleaning cuts and preselection
	//MET > 150
	string metcut = "(selCMet > 150)";
	//pTs < 150 - always use index 0 for rjr variables
	string ptscut = "(rjrPTS0 < 150)";
	//triggers applied - 2018
	string triggers = "(Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || Trigger_PFMET120_PFMHT120_IDTight_PFHT60 == 1 || Trigger_PFMET120_PFMHT120_IDTight == 1)";
	//event cleaning filters applied
	string met_filters = "(Flag_BadChargedCandidateFilter == 1 && Flag_BadPFMuonDzFilter == 1 && Flag_BadPFMuonFilter == 1 && Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_ecalBadCalibFilter == 1 && Flag_eeBadScFilter == 1 && Flag_goodVertices == 1 && Flag_hfNoisyHitsFilter == 1)";
	vector<string> branches = {"rjr_Rs","rjr_Ms","selCMet","rjrPTS","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","Trigger_PFMETNoMu120_PFMHTNoMu120_IDTight","Trigger_PFMET120_PFMHT120_IDTight_PFHT60","Trigger_PFMET120_PFMHT120_IDTight","Flag_BadChargedCandidateFilter","Flag_BadPFMuonDzFilter","Flag_BadPFMuonFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_ecalBadCalibFilter","Flag_eeBadScFilter","Flag_goodVertices","Flag_hfNoisyHitsFilter","nSelPhotons","selPhoEta","passNPhoGe1SelectionEndcapNonIsoCR","selPho_nonIsoANNScore"};
	vector<double> msbins = {1000,2000,3000};
	vector<double> rsbins = {0.15,0.3,1};
	int n_msbins = msbins.size()-1;
	int n_rsbins = rsbins.size()-1;

	
	string ofilename = "test_photonID";
	if(ofilename_extra != "") ofilename += "_"+ofilename_extra;
	ofilename += ".root";
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;
	MakeFileMaps();
	auto ch = std::unique_ptr<TChain>(new TChain("kuSkimTree"));
	vector<string> filekeys = {"QCD","GJets"};//,"SMS_GlGl_0p5","SMS_GlGl_0p1"};
	for(int k = 0; k < filekeys.size(); k++){
		cout << "process " << filekeys[k] << endl;
		vector<string> bkgfiles = file_map[filekeys[k]];
		for(auto file : bkgfiles) ch->Add(file.c_str());
		ROOT::RDataFrame df0 = ROOT::RDataFrame(*ch.get(), branches);
		//for skims v45/v46
		auto df1 = df0.Define("passNPhoGe1SelectionEndcapNonIsoCR0","passNPhoGe1SelectionEndcapNonIsoCR[0]").
				Define("passNPhoGe1SelectionBeamHaloCR0","passNPhoGe1SelectionBeamHaloCR[0]").
				Define("nonIsoScoreEndcap","selPho_nonIsoANNScore[selPhoEta > 1.479 || selPhoEta < -1.479]").
				Define("nEndcapPhotons","selPhoEta[selPhoEta > 1.479 || selPhoEta < -1.479].size()").
				Define("nEndcapPhotons_nonIso","selPhoEta[(selPhoEta > 1.479 || selPhoEta < -1.479) && (selPho_nonIsoANNScore > 0.99142313)].size()").
				Define("rjr_Ms0","rjr_Ms[0]").Define("rjr_Rs0","rjr_Rs[0]").Define("rjrPTS0","rjrPTS[0]");

		//look at met and event filters
		string histname = "selCMet_"+filekeys[k];
		auto h_met = df1.Histo1D({histname.c_str(),histname.c_str(),50,0,800},"selCMet");
		hists1d.push_back(h_met);
		
		histname = "rjrPTS_"+filekeys[k];
		auto h_pts = df1.Histo1D({histname.c_str(),histname.c_str(),50,0,800},"rjrPTS0");
		hists1d.push_back(h_pts);

		for(int b = 0; b < branches.size(); b++){
			string branchname = branches[b];
			if(branchname.find("Flag") == string::npos && branchname.find("Trigger") == string::npos) continue;

			auto h = df1.Histo1D({(branchname+"_"+filekeys[k]).c_str(), (branchname+"_"+filekeys[k]).c_str(), 2, 0, 2},branchname);
			hists1d.push_back(h);

		}
		//check preselection efficiencies
		auto df_triggers = df1.Filter(triggers,"triggers");
		auto df_filters = df1.Filter(met_filters,"met_filters");
		auto df_met = df1.Filter(metcut,metcut);
		auto df_ptscut = df1.Filter(ptscut,ptscut);
		string baselinecut = metcut + " && " + ptscut + " && " + triggers+ " && " + met_filters;
		auto df_presel = df1.Filter(baselinecut,"baseline");

	

		//channels
		auto df_ge1pho =   std::make_unique<ROOT::RDF::RNode>(df_presel.Filter("nSelPhotons >= 1","ge1nSelPho"));
		auto df_1pho =   std::make_unique<ROOT::RDF::RNode>(df_presel.Filter("nSelPhotons == 1","1nSelPho"));
		auto df_ge2pho = std::make_unique<ROOT::RDF::RNode>(df_presel.Filter("nSelPhotons >= 2","ge2nSelPho"));

		map<string, std::unique_ptr<ROOT::RDF::RNode>> ch_sels;
		ch_sels["ge1pho"] = std::move(df_ge1pho);
		ch_sels["1pho"] = std::move(df_1pho);
		ch_sels["ge2pho"] = std::move(df_ge2pho);
		
		for(auto ch = ch_sels.begin(); ch != ch_sels.end(); ch++){
			string ch_name = ch->first;
			cout << "channel " << ch_name << endl;

			//noniso score
			histname = "nonIsoPredScore_endcap_"+filekeys[k]+"_"+ch_name;
			auto h1 = ch->second->Histo1D({histname.c_str(),histname.c_str(),50,0,1.05},"nonIsoScoreEndcap");
			hists1d.push_back(h1);

			//check (raw) efficiency of nonIso score cut on photons
			auto df_EEnonIsoCR = ch->second->Filter("passNPhoGe1SelectionEndcapNonIsoCR0 == 1","nonIsoDNNCRCut_"+ch_name);

			//calculate the (weighted) efficiency of the CR definition at the event level
			auto weight_all = ch->second->Sum("evtFillWgt");
			auto weight_pass = df_EEnonIsoCR.Sum("evtFillWgt");	
			cout << "passNPhoGe1SelectionEndcapNonIsoCR " << weight_pass.GetValue() << " all (" << ch_name << ") " << weight_all.GetValue() << " Weighted efficiency : " << (weight_pass.GetValue() / weight_all.GetValue()) * 100 << "%" << endl;
			cout << endl;
			
			auto yields_EEnonIsoCR = df_EEnonIsoCR.
				Histo2D({("yields_EEnonIsoCR"+filekeys[k]).c_str(),("yields_EEnonIsoCR"+filekeys[k]+";Ms;Rs").c_str(),n_msbins,&msbins[0],n_rsbins,&rsbins[0]},"rjr_Ms0","rjr_Rs0");
			hists2d.push_back(yields_EEnonIsoCR);
		}
		df0.Report()->Print();
		cout << endl; cout << endl;
	}
	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();

	ofile->Close();
	cout << "Wrote histograms to " << ofile->GetName() << endl;

}
