#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include "TTreeInterface.h"

//runs on KU skims!!
void DetBkgAnalysis(string file, string tag = "CMS", bool recreatecsv = false, string ofilename_extra = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	string physbkgLabel = "SC_trueLabel_"+tag+" == 1";
	string bhLabel = "SC_trueLabel_"+tag+" == 2";
	string spikeLabel = "SC_trueLabel_"+tag+" == 3";

	//rh maps	
	TTreeInterface TI(file,"kuSkimTree");
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
	skipbranch["SC_seedTime_"+tag] = -999;
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
	

	ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ');

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

	//string predspikes = "SC_predScore_spike_"+tag+" > 0.96000487";
	string earlytimes = "SC_seedTime_"+tag+" < -2";

	auto dfphi = df.Define("SC_PhiSig_"+tag,"sqrt(SC_PhiVar_"+tag+")");
	auto dfeta = dfphi.Define("SC_EtaSig_"+tag,"sqrt(SC_EtaVar_"+tag+")");
	auto df0 = dfeta.Filter("SC_EtaCenter_"+tag+" > -1.5 && SC_EtaCenter_"+tag+" < 1.5","barrelOnly");
	auto df_endcap = dfeta.Filter("!(SC_EtaCenter_"+tag+" > -1.5 && SC_EtaCenter_"+tag+" < 1.5)","endcapOnly");
	auto df1 = df0.Filter(physbkgLabel,"physbkgCR");
	auto df2 = df0.Filter(bhLabel,"bhCR");
	auto df3 = df0.Filter(spikeLabel,"spikeCR");

	auto df_early = df0.Filter(earlytimes,"earlyTimes");
	auto df_bhfilter = df_early.Filter("Flag_globalSuperTightHalo2016Filter == 1","BHFilter");

	//nominal hists
	auto h_nom_tc_ec = df_early. 
				Histo2D({("SC_seedTime_EtaCenter_"+tag).c_str(),("SC_seedTime_EtaCenter_"+tag+";seedTime;EtaCenter;a.u.").c_str(),50,-16,0,50,-1.6,1.6},"SC_seedTime_"+tag,"SC_EtaCenter_"+tag);
	hists2d.push_back(h_nom_tc_ec);
	
	auto h_nom_tc_ec_endcap = df_endcap. 
				Histo2D({("SC_seedTime_EtaCenter_"+tag+"_endcap").c_str(),("SC_seedTime_EtaCenter_"+tag+"_endcap"+";seedTime;EtaCenter;a.u.").c_str(),50,-16,16,50,-3.,3.},"SC_seedTime_"+tag,"SC_EtaCenter_"+tag);
	hists2d.push_back(h_nom_tc_ec_endcap);

	//nominal bh flag applied
	auto h_nom_tc_ec_bhfilter = df_bhfilter. 
				Histo2D({("SC_seedTime_EtaCenter_"+tag+"_bhFilterApplied").c_str(),("SC_seedTime_EtaCenter_"+tag+"_bhFilterApplied"+";seedTime;EtaCenter;a.u.").c_str(),50,-16,0,50,-1.6,1.6},"SC_seedTime_"+tag,"SC_EtaCenter_"+tag);
	hists2d.push_back(h_nom_tc_ec_bhfilter);
	

	//CR hists
	auto h_time_eta_trueBH = df2. 
				Histo2D({("SC_TimeCenter_EtaCenter_trueBH_"+tag).c_str(),("SC_TimeCenter_EtaCenter_trueBH_"+tag+";TimeCenter;EtaCenter;a.u.").c_str(),50,-16,16,50,-1.6,1.6},"SC_seedTime_"+tag,"SC_EtaCenter_"+tag);
	hists2d.push_back(h_time_eta_trueBH);
	auto h_etavar_phivar_trueBH = df2.
				Histo2D({("SC_EtaSig_PhiSig_trueBH_"+tag).c_str(),("SC_EtaSig_PhiSig_trueBH_"+tag+";EtaSig;PhiSig;a.u.").c_str(),50,0.,0.06,50,0,0.3},"SC_EtaSig_"+tag,"SC_PhiSig_"+tag);
	hists2d.push_back(h_etavar_phivar_trueBH);
	auto h_etavar_phivar_truespikes = df3.
				Histo2D({("SC_EtaSig_PhiSig_truespikes_"+tag).c_str(),("SC_EtaSig_PhiSig_truespikes_"+tag+";EtaSig;PhiSig;a.u.").c_str(),50,0.,0.06,50,0,0.3},"SC_EtaSig_"+tag,"SC_PhiSig_"+tag);
	hists2d.push_back(h_etavar_phivar_truespikes);
	auto h_etavar_phivar_truephysBkg = df1.
				Histo2D({("SC_EtaSig_PhiSig_truephysBkg_"+tag).c_str(),("SC_EtaSig_PhiSig_truephysBkg_"+tag+";EtaSig;PhiSig;a.u.").c_str(),50,0.,0.06,50,0,0.3},"SC_EtaSig_"+tag,"SC_PhiSig_"+tag);
	hists2d.push_back(h_etavar_phivar_truephysBkg);
	
	auto h_time_eta_truespikes = df3. 
				Histo2D({("SC_TimeCenter_EtaCenter_truespikes_"+tag).c_str(),("SC_TimeCenter_EtaCenter_truespikes_"+tag+";TimeCenter;EtaCenter;a.u.").c_str(),50,-16,16,50,-1.6,1.6},"SC_seedTime_"+tag,"SC_EtaCenter_"+tag);
	hists2d.push_back(h_time_eta_truespikes);


	string ofilename = csvname_base+"_detbkg";
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
