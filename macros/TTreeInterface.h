//////////////////////////////////////////////////
//Thanks to Justin Anguiano (2024) for the basis 
//for this TTree exploding/unrolling code
//////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <utility>
#include <set>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

typedef std::map<std::string, std::pair<std::string, std::string> > strIdxMap ;
typedef std::map<std::string, std::pair<std::vector<double>*,std::vector<double>*>> vecIdxMap;
typedef std::map<std::string, std::pair<std::vector<std::vector<double>>*,std::vector<double>*>> vecvecIdxMap;

class TTreeInterface{
	public:
	
	TFile* _f{};
	std::vector<std::string> _filelist{};
	TTree* _ttree{0};
	TChain* _tchain{0};
	std::string _filename{};
	std::string _treename{};
	std::set<std::string> _idxSet{}; //set of idx branches
	std::set<std::string> _brSet{}; //set of target branches

	std::string _subobjBranchName = "";
	
	//index chasing var
	strIdxMap _idxMaps{};//stores name associations, keep names around for debug and printing
	strIdxMap _idxSubMaps{};
	vecIdxMap _branchMaps{}; //auto load datastructures with same key mapping
	vecvecIdxMap _branchSubMaps{};	
	
	vector<double>* _nsubobj = nullptr;
	
	//constructors to auto load trees/files from strings
	TTreeInterface(std::string fileName, std::string treeName);
	TTreeInterface(std::vector<std::string> fileList, std::string treeName);
	void SetNSubBranch( std::string branchname);
	
	//do index chasing for quantities in newer versions
	//input idx branchname associated to idx of parent array, output label is csv header for this mapping
	void MapIdx( std::string idxBranchName, std::string targetBranch, std::string outLabel );
	void MapIdxToSubIdx( std::string idxBranchName, std::string targetBranch, std::string outLabel );
	double RetrieveMapValue( std::string mapkey, int dim_idx);
	double RetrieveSubMapValue( std::string mapkey, int dim_idx, int sub_idx);
	bool CheckIdxUniqueMap( std::string idxBranchName);
	bool CheckTargetUniqueMap( std::string targetBranchName);
	std::vector<double>* GetBranchAlreadyAssigned( std::string idxBranchName );
	std::vector<std::vector<double>>* GetSubBranchAlreadyAssigned( std::string idxBranchName );
	std::vector<double>* GetTargetBranchAlreadyAssigned( std::string targetBranchName );
	
	//make a flattened csv
	void CreateFlattenedCSV( std::vector<std::string> branchList, std::vector<string> subBranchList, std::string csvname );
	void CreateFlattenedCSV( std::vector<std::string> branchList, std::vector<string> subBranchList, std::string csvname, std::vector<string> evtbranchlist);
	void CreateFlattenedCSV( std::vector<std::string> branchList, std::string csvname);
	
	void SetSkipBranches(const map<string, double>& skipbranches);
	map<string, double> _skipBranch;
		
};
TTreeInterface::TTreeInterface(std::string fileName, std::string treeName){
	_filename = fileName;
	_treename = treeName;
	_f = TFile::Open(_filename.c_str());
	std::cout<<"Loading "<<_treename<<" from "<< _filename<<"\n";
	_ttree = (TTree*)_f->Get(_treename.c_str());
}

TTreeInterface::TTreeInterface(std::vector<std::string> fileList, std::string treeName){
	_filelist = fileList;
	_treename = treeName;
	_tchain = new TChain(_treename.c_str());
	std::cout<<"Loading "<<_treename<<" from "<<_filelist.size()<<" files \n";
	for( int i=0; i<_filelist.size(); i++){
		_tchain->AddFile(_filelist[i].c_str());
	}
	_ttree = _tchain;
}
bool TTreeInterface::CheckIdxUniqueMap( std::string idxBranchName){
	//can only register the idx array to 1 object, so we need point to the original if it exists already
	bool branchExists=false;
		for (auto& idxstr : _idxSet) {
     
		  if( idxstr == idxBranchName ){
		  	branchExists = true;
		  }	   
		}
	std::cout<<"Checking if branch "<< idxBranchName <<" address is already assigned: branchExists="<<branchExists<<"\n";
	return branchExists;
	
}
bool TTreeInterface::CheckTargetUniqueMap( std::string targetBranchName){
	//can only register the idx array to 1 object, so we need point to the original if it exists already
	bool branchExists=false;
		for (auto& brstr : _brSet) {
     
		  if( brstr == targetBranchName ){
		  	branchExists = true;
		  }	   
		}
	std::cout<<"Checking if branch "<< targetBranchName <<" address is already assigned: branchExists="<<branchExists<<"\n";
	return branchExists;
	
}
std::vector<double>* TTreeInterface::GetTargetBranchAlreadyAssigned( std::string targetBranchName ){
	for(strIdxMap::iterator iter = _idxMaps.begin(); iter != _idxMaps.end(); ++iter)
		{
			std::string key = iter->first;
		  if( _idxMaps[ key].second == targetBranchName ){
		  	return _branchMaps[key].second;
		  }
		}
	return 0;

}
std::vector<double>* TTreeInterface::GetBranchAlreadyAssigned( std::string idxBranchName ){
	for(strIdxMap::iterator iter = _idxMaps.begin(); iter != _idxMaps.end(); ++iter)
		{
			std::string key = iter->first;
		  if( _idxMaps[ key].first == idxBranchName ){
		  	return _branchMaps[key].first;
		  }
		}
	return 0;

}
std::vector<std::vector<double>>* TTreeInterface::GetSubBranchAlreadyAssigned( std::string idxBranchName ){
	for(strIdxMap::iterator iter = _idxSubMaps.begin(); iter != _idxSubMaps.end(); ++iter)
		{
			std::string key = iter->first;
		  if( _idxSubMaps[ key].first == idxBranchName ){
//cout << "found idx branch " << key << " for idxbranchname " << idxBranchName << endl;
		  	return _branchSubMaps[key].first;
		  }
		}
	return 0;

}
void TTreeInterface::MapIdx( std::string idxBranchName, std::string targetBranch, std::string outLabel ){
	std::vector<double>* v1 =0;
	std::vector<double>* v2 =0;
	_branchMaps[ outLabel ] = std::make_pair( v1, v2 );//init pointers to 0
	
	std::cout<<"setting idx unroll branch address of: "<< outLabel <<"\n";
	if( !CheckIdxUniqueMap( idxBranchName ) ){
		_ttree->SetBranchAddress(idxBranchName.c_str(),&(_branchMaps[outLabel].first) );
	}else{
		_branchMaps[ outLabel ].first = GetBranchAlreadyAssigned( idxBranchName );
	}
	if( !CheckTargetUniqueMap(targetBranch) ){
    		_ttree->SetBranchAddress(targetBranch.c_str(), &(_branchMaps[outLabel].second) );
	}
	else{
		_branchMaps[outLabel].second = GetTargetBranchAlreadyAssigned( targetBranch );
	}
     //need to do this after branch setting in case outlabel comes before branches that are already assigned in the map 
    _idxMaps[ outLabel ] = std::make_pair( idxBranchName, targetBranch );
    _idxSet.insert( idxBranchName );
    _brSet.insert( targetBranch );
}

//issue: if targetBranch's address is set to a different output map, it will reset
//if branch address is already set, map _branchSubMaps[outLabel] to the vector it is set to
//similar to what's done for idxBranchName
void TTreeInterface::MapIdxToSubIdx( std::string idxBranchName, std::string targetBranch, std::string outLabel ){
	std::vector<std::vector<double>>* v1 =0;
	std::vector<double>* v2 =0;
	_branchSubMaps[ outLabel ] = std::make_pair( v1, v2 );//init pointers to 0
	
	std::cout<<"setting idx unroll branch address of: "<< outLabel <<"\n";
	if( !CheckIdxUniqueMap( idxBranchName ) ){
		//if idxbranchname doesn't exist in ttree, create that branch and fill it with idx of vector of entries
		_ttree->SetBranchAddress(idxBranchName.c_str(),&(_branchSubMaps[outLabel].first) );
	}else{
		_branchSubMaps[ outLabel ].first = GetSubBranchAlreadyAssigned( idxBranchName );
	}
	if( !CheckTargetUniqueMap(targetBranch) ){
    		_ttree->SetBranchAddress(targetBranch.c_str(), &(_branchSubMaps[outLabel].second) );
	}
	else{
		_branchSubMaps[outLabel].second = GetTargetBranchAlreadyAssigned( targetBranch );
	}
     //need to do this after branch setting in case outlabel comes before branches that are already assigned in the map 
    _idxSubMaps[ outLabel ] = std::make_pair( idxBranchName, targetBranch );
    _idxSet.insert( idxBranchName );
    _brSet.insert( targetBranch );
}

double TTreeInterface::RetrieveMapValue(std::string mapkey, int dim_idx){
	
	double val{};
	double idx{};
	if( (_branchMaps[mapkey].first)->size()==0 || (_branchMaps[mapkey].second)->size()==0){
		return -999;
	}
	
	if( (_branchMaps[mapkey].first)->at(dim_idx) < 0 ){
		if(mapkey.find("ID") != string::npos)
			return -999;
		else
			return -1;
	}
	else{
		idx = (_branchMaps[mapkey].first)->at(dim_idx);
		val = (_branchMaps[mapkey].second)->at(idx);
		return val;
	}
	
}
double TTreeInterface::RetrieveSubMapValue(std::string mapkey, int dim_idx, int sub_idx){
	
	double val{};
	double idx{};
	if( (_branchSubMaps[mapkey].first)->size()==0 || (_branchSubMaps[mapkey].second)->size()==0){
		return -999;
	}
	if( (_branchSubMaps[mapkey].first)->at(dim_idx).at(sub_idx) < 0 ){
		if(mapkey.find("ID") != string::npos)
			return -999;
		else
			return -1;
	}
	else{
		idx = (_branchSubMaps[mapkey].first)->at(dim_idx).at(sub_idx);
		val = (_branchSubMaps[mapkey].second)->at(idx);
		return val;
	}
//typedef std::map<std::string, std::pair<std::vector<std::vector<double>>*,std::vector<std::vector<double>>*>> vecvecIdxMap;
	
}

void TTreeInterface::SetNSubBranch( std::string branchname){
	_ttree->SetBranchAddress(branchname.c_str(),&_nsubobj);
	_subobjBranchName = branchname;
}

void TTreeInterface::SetSkipBranches(const map<string, double>& skipbranches){
	_skipBranch.clear();
	_skipBranch = skipbranches;
}


void TTreeInterface::CreateFlattenedCSV( std::vector<std::string> branchList, std::vector<string> subBranchList, std::string csvname, std::vector<string> evtbranchlist){
	// loop over selected branches
	// expect vector of primitive types and flatten
	// also record a relative event ID

	double evt;	
	_ttree->SetBranchAddress("evt",&evt);
	//prepare space delimited csv file output
	std::ofstream ocsv;
  	ocsv.open(csvname);
  	ocsv << "evtidx ";
	ocsv<<"jetidx ";
	if(_nsubobj != nullptr){
		ocsv<<"subclidx ";
		ocsv<<_subobjBranchName<<" ";
	}
	//evt observables
	std::vector<double> evtobs(evtbranchlist.size(),0);
	for(int i = 0; i < evtbranchlist.size(); i++){
		ocsv << evtbranchlist[i]<<" ";
		std::cout<<"setting branch address of: "<< evtbranchlist[i]<<"\n";
		_ttree->SetBranchAddress(evtbranchlist[i].c_str(),&(evtobs[i]));
	}
	//object observables
	//idx mapping abstraction
	for(strIdxMap::iterator iter = _idxMaps.begin(); iter != _idxMaps.end(); ++iter)
	{
	  std::string key =  iter->first;
	  std::cout<<"found mapping: "<<key<<" "<<_idxMaps[key].first<<":"<<_idxMaps[key].second<<"\n"; 
	  //print headers
	  ocsv<<key<<" ";
	}



	//do for nsubbranch if specified
	if(_subobjBranchName != ""){
		branchList.push_back(_subobjBranchName);
	}
	//create variables for ttree to load into
	std::vector< std::vector<double>* > branchVec(branchList.size(), 0);
	for(int i=0; i<branchList.size(); i++){
		//set dynamic addresses
		std::cout<<"setting branch address of: "<< branchList[i]<<"\n";
		if(branchList[i] == _subobjBranchName){
			branchVec[i] = _nsubobj;
		}
		else
			_ttree->SetBranchAddress(branchList[i].c_str(),&branchVec[i]);
		//write the csv headers
		ocsv << branchList[i] << " ";
	}

	std::vector< std::vector<std::vector<double>>* > subBranchVec(subBranchList.size(), 0);
	if(_nsubobj != nullptr){
		//subobject observables
		for(strIdxMap::iterator iter = _idxSubMaps.begin(); iter != _idxSubMaps.end(); ++iter)
		{
		  std::string key =  iter->first;
		  std::cout<<"found submapping: "<<key<<" "<<_idxSubMaps[key].first<<":"<<_idxSubMaps[key].second<<"\n"; 
		  ocsv<<key<<" ";		
		}
		
		for(int i=0; i<subBranchList.size(); i++){
			//set dynamic addresses
			std::cout<<"setting subBranch address of: "<< subBranchList[i]<<"\n";
			_ttree->SetBranchAddress(subBranchList[i].c_str(),&subBranchVec[i]);
			//write the csv headers
			ocsv << subBranchList[i];
			if(i < subBranchList.size()-1) ocsv << " ";
		}
	}
	ocsv << "\n";

	for(auto it = _skipBranch.begin(); it != _skipBranch.end(); it++){
		cout << "skipping objects with " << it->first << " == " << it->second << endl;
	}
	
	Long64_t nentries = _ttree->GetEntries();
	std::cout<<"Looping over "<<nentries<<" events\n";
	int dim;
	for(Long64_t i=0; i< nentries; i++){
		_ttree->GetEntry(i);
		//cout << "evt " << i << " ttree evt " << evt << endl;
		//do event observables
			
		
		//how many objects to flatten? grab the first one on the list
		dim = (branchVec[0])->size();
		//debug print
		//std::cout<<"found # objs "<<dim<< " from " << branchList[0] << " branch" << " skipBranch has " << branchVec[4]->size() << " objs" << endl;
		//unroll the vector
		for(int j=0; j<dim; j++){
			//cout << "obj # " << j << " has # subobjs " << _nsubobj->at(j) << endl;
			//do skipbranch selection
			bool skip = false;
			for(auto it = _skipBranch.begin(); it != _skipBranch.end(); it++){
//cout << "skipbranch " << it->first << endl;
				auto iit = find(branchList.begin(), branchList.end(), it->first);
				if(iit == branchList.end()) continue;
				//get index for branchvec
				int idx = distance(branchList.begin(), iit);
//cout << "j " << j << " branchvec idx " << idx << " branchvec size " << branchVec[idx]->size() << " branchlist[idx] " << branchList[idx] << endl;
//cout << "branchvec[idx]->at(j) " << branchVec[idx]->at(j) << endl;
				if(branchVec[idx]->at(j) == it->second){
					skip = true;
				}
			}
			if(skip) continue;
			if(_nsubobj == nullptr){
				//cout << "Error: set branch of number of subobjects with SetNSubBranch(branchname)" << endl;
				//return;
				//if newever versions extract the values from the idx mapping
				//event index	
				ocsv<<i<<" ";
				//jet index
				ocsv<<j<<" ";
				for(int o = 0; o < evtobs.size(); o++)
					ocsv<<evtobs[o]<<" ";		
				for(strIdxMap::iterator iter = _idxMaps.begin(); iter != _idxMaps.end(); ++iter)
				{
					//cout << iter->first << " "  << RetrieveMapValue( iter->first, j) << endl;
					ocsv<< RetrieveMapValue( iter->first, j) << " ";
				}
				//unroll objects	
				for(int k=0; k<branchVec.size(); k++){
					//debug print
					//std::cout<<branchList[k]<<" "<<branchVec[k]->at(j)<<" \n";
					//write branch quantities and evtid
					ocsv<<branchVec[k]->at(j)<<" "; 
				}
				ocsv<<"\n";
			}
			else{
			//cout << "obj # " << j << " has # subobjs " << _nsubobj->at(j) << endl;
			//TODO - debug in RunClustering and remove when fixed 
			if(subBranchVec[0]->at(j).size() != _nsubobj->at(j)){
				//cout << "_nsubobj " << _nsubobj->at(j) << " size from subbranchvec[0] " << subBranchVec[0]->at(j).size() << endl;
				continue;
			}
				for(int n = 0; n < _nsubobj->at(j); n++){
					//cout << "subobj #" << n << endl;
					//event index	
					ocsv<<i<<" ";
					//jet index
					ocsv<<j<<" ";
					//subcluster index
					ocsv<<n<<" ";
					ocsv<<_nsubobj->at(j)<<" ";
					for(int o = 0; o < evtobs.size(); o++)
						ocsv<<evtobs[o]<<" ";		
					//if newever versions extract the values from the idx mapping
					for(strIdxMap::iterator iter = _idxMaps.begin(); iter != _idxMaps.end(); ++iter)
					{
						//cout << iter->first << " "  << RetrieveMapValue( iter->first, j) << endl;
						ocsv<< RetrieveMapValue( iter->first, j) << " ";
					}
		
					//unroll objects	
					for(int k=0; k<branchVec.size(); k++){
						//debug print
						//std::cout<<branchList[k]<<" "<<branchVec[k]->at(j)<<" \n";
						//write branch quantities and evtid
						ocsv<<branchVec[k]->at(j)<<" "; 
					}
					//unroll subobjects
					for(strIdxMap::iterator iter = _idxSubMaps.begin(); iter != _idxSubMaps.end(); ++iter)
					{
						//cout << iter->first << " "  << RetrieveSubMapValue( iter->first, j, n) << endl;
						ocsv<< RetrieveSubMapValue( iter->first, j, n) << " ";
					}
					for(int k = 0; k < subBranchVec.size(); k++){
						//debug print
//cout << "subbranchvec " << subBranchList[k] << " size " << subBranchVec[k]->size() << " j " << j << " size " << subBranchVec[k]->at(j).size() << endl;	
						//std::cout<<subBranchList[k]<<" "<<subBranchVec[k]->at(j).at(n)<<" \n";
						//write subBranch quantities and evtid
						ocsv<<subBranchVec[k]->at(j).at(n) << " "; 
					}
					ocsv<<"\n";
				}
			}	
		}
			
	}
	ocsv.close();
	cout << "Writing flattened TTree to " << csvname << endl;
}


/*///////////////////////
Generic CSV global methods
////////////////////////*/
std::vector<std::string> ReadList( std::string fileList){
	std::vector<std::string> filevec{};
	std::ifstream f; 
	f.open(fileList);
	std::string ifile{};
	if( f.is_open()){
		while( f >> ifile ) { 
			std::cout << "adding item: "<<ifile<<"\n";
			filevec.push_back(ifile);
		}
		f.close();
	}		
	else{ std::cout<< "filelist not found\n"; }
	return filevec;
}






