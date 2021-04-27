//This file is used for studying the performance of leptons and jets in semileptonic ttbar events at FCC-ee
//To run: g++ -o PerformanceStudiesNew -std=gnu++17 -I /Users/Julie/ROOT/include PerformanceStudiesNew.cc `root-config --cflags --libs`                                                        
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "iostream"
#include "fstream"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TProfile.h"
#include <unordered_set>
#include <numeric>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "Math/Vector4D.h"

using namespace std;

TFile *outFile;

int main()
{
  Float_t tol = 1.0e-16;

  TChain* chain = new TChain("events");
  chain->Add("outputs/semileptonic/p8_ee_ttbar_semi_ecm365.root");

  std::vector<TString> JetAlgos;
  JetAlgos.emplace_back(TString("kt"));
  JetAlgos.emplace_back(TString("durham"));
  JetAlgos.emplace_back(TString("ee_antikt"));
  JetAlgos.emplace_back(TString("cambridge"));
  JetAlgos.emplace_back(TString("jade"));
  JetAlgos.emplace_back(TString("valencia"));
    
  //Defining relevant branches in the file
  vector<float> *MC_p=0;
  vector<float> *MC_px=0;
  vector<float> *MC_py=0;
  vector<float> *MC_pz=0;
  vector<float> *MC_charge=0;
  vector<float> *MC_mass=0;
  vector<float> *MC_e=0;
  vector<float> *MC_pdg=0;
  vector<float> *MC_status=0;
  
  vector<float> *RP_p=0;
  vector<float> *RP_px=0;
  vector<float> *RP_py=0;
  vector<float> *RP_pz=0;
  vector<float> *RP_charge=0;
  vector<float> *RP_mass=0;
  vector<float> *RP_e=0;
  
  vector<float> *MET_p=0;
  vector<float> *MET_px=0;
  vector<float> *MET_py=0;
  vector<float> *MET_pz=0;
  vector<float> *MET_charge=0;
  vector<float> *MET_mass=0;
  vector<float> *MET_e=0;
  
  vector<float> *RPMC_pdg=0;
  vector<int> *RPMC_index=0;
  vector<int> *RPMC_parentindex=0;
  vector<int> *MC_parent1=0;
  vector<int> *MC_parent2=0;
  vector<int> *MC_daughter1=0;
  vector<int> *MC_daughter2=0;

  vector<int> *RPrest_association=0;
  vector<int> *MCfinal_association=0;
  vector<int> *MCparton_association=0;

  vector<float> *recojets_jetalgo_px=0;
  vector<float> *recojets_jetalgo_py=0;
  vector<float> *recojets_jetalgo_pz=0;
  vector<float> *recojets_jetalgo_e=0;
  vector<vector<int>> *recojetconstituents_jetalgo=0;
  vector<float> *particlejets_jetalgo_px=0;
  vector<float> *particlejets_jetalgo_py=0;
  vector<float> *particlejets_jetalgo_pz=0;
  vector<float> *particlejets_jetalgo_e=0;
  vector<vector<int>> *particlejetconstituents_jetalgo=0;
  vector<float> *partonjets_jetalgo_px=0;
  vector<float> *partonjets_jetalgo_py=0;
  vector<float> *partonjets_jetalgo_pz=0;
  vector<float> *partonjets_jetalgo_e=0;
  vector<vector<int>> *partonjetconstituents_jetalgo=0;


  chain->SetBranchAddress("MC_px", & MC_px);
  chain->SetBranchAddress("MC_py", & MC_py);
  chain->SetBranchAddress("MC_pz", & MC_pz);
  chain->SetBranchAddress("MC_e", & MC_e);
  chain->SetBranchAddress("MC_pdg", & MC_pdg);
  chain->SetBranchAddress("MC_status", & MC_status);
  chain->SetBranchAddress("MC_mass", & MC_mass);

  chain->SetBranchAddress("RP_p", & RP_p);
  chain->SetBranchAddress("RP_px", & RP_px);
  chain->SetBranchAddress("RP_py", & RP_py);
  chain->SetBranchAddress("RP_pz", & RP_pz);
  chain->SetBranchAddress("RP_charge", & RP_charge);
  chain->SetBranchAddress("RP_mass", & RP_mass);
  chain->SetBranchAddress("RP_e", & RP_e);
  
  chain->SetBranchAddress("MET_p", & MET_p);
  chain->SetBranchAddress("MET_px", & MET_px);
  chain->SetBranchAddress("MET_py", & MET_py);
  chain->SetBranchAddress("MET_pz", & MET_pz);
  chain->SetBranchAddress("MET_charge", & MET_charge);
  chain->SetBranchAddress("MET_mass", & MET_mass);
  
  chain->SetBranchAddress("RPMC_pdg", & RPMC_pdg);
  chain->SetBranchAddress("RPMC_index", & RPMC_index);
  chain->SetBranchAddress("RPMC_parentindex", & RPMC_parentindex);

  chain->SetBranchAddress("MC_parent1", & MC_parent1);
  chain->SetBranchAddress("MC_parent2", & MC_parent2);
  chain->SetBranchAddress("MC_daughter1", & MC_daughter1);
  chain->SetBranchAddress("MC_daughter2", & MC_daughter2);

  chain->SetBranchAddress("RPrest_association", & RPrest_association);
  chain->SetBranchAddress("MCfinal_association", & MCfinal_association);
  chain->SetBranchAddress("MCparton_association", & MCparton_association);

  for (auto& jetalgo : JetAlgos) cout << jetalgo << endl;

  for (auto& jetalgo : JetAlgos) {
    outFile = new TFile("p8_analysis_"+jetalgo+".root", "RECREATE");

    chain->SetBranchAddress("recojets_"+jetalgo+"_px", & recojets_jetalgo_px);
    chain->SetBranchAddress("recojets_"+jetalgo+"_py", & recojets_jetalgo_py);
    chain->SetBranchAddress("recojets_"+jetalgo+"_pz", & recojets_jetalgo_pz);
    chain->SetBranchAddress("recojets_"+jetalgo+"_e", & recojets_jetalgo_e);
    chain->SetBranchAddress("recojetconstituents_"+jetalgo, & recojetconstituents_jetalgo);
    
    chain->SetBranchAddress("particlejets_"+jetalgo+"_px", & particlejets_jetalgo_px);
    chain->SetBranchAddress("particlejets_"+jetalgo+"_py", & particlejets_jetalgo_py);
    chain->SetBranchAddress("particlejets_"+jetalgo+"_pz", & particlejets_jetalgo_pz);
    chain->SetBranchAddress("particlejets_"+jetalgo+"_e", & particlejets_jetalgo_e);
    chain->SetBranchAddress("particlejetconstituents_"+jetalgo, & particlejetconstituents_jetalgo);

    chain->SetBranchAddress("partonjets_"+jetalgo+"_px", & partonjets_jetalgo_px);
    chain->SetBranchAddress("partonjets_"+jetalgo+"_py", & partonjets_jetalgo_py);
    chain->SetBranchAddress("partonjets_"+jetalgo+"_pz", & partonjets_jetalgo_pz);
    chain->SetBranchAddress("partonjets_"+jetalgo+"_e", & partonjets_jetalgo_e);
    chain->SetBranchAddress("partonjetconstituents_"+jetalgo, & partonjetconstituents_jetalgo);
  
    auto *hprofLepton  = new TProfile("hprofLepton","Profile for absolute energy difference",100,0,200,0,100);
    TH1F *hJetCosOpeningAngle = new TH1F("hJetCosOpeningAngle", "Cosine to the Matching Angle for jets", 100, 0.01, 1.0);
    auto *hprofJetAbs  = new TProfile("hprofJetAbs","Profile plot for absolute energy difference",100,0,200,0,100);
    auto *hprofJetRel  = new TProfile("hprofJetRel","Profile plot for relative energy difference",100,0,200,0,100,"s"); //Compared to quark level
    TH1F *hJetCosOpeningAngleMC = new TH1F("hJetCosOpeningAngleMC", "Cosine to the Matching Angle for jets", 100, 0.01, 1.0);
    auto *hprofJetAbsMC  = new TProfile("hprofJetAbsMC","Profile plot for absolute energy difference",50,0,200,-100,100);
    auto *hprofJetRelMC  = new TProfile("hprofJetRelMC","Profile plot for relative energy difference",50,0,200,-100,100,"s"); //Compared to MC level jets
    auto *hprofJetRelMCUnique  = new TProfile("hprofJetRelMCUnique","Profile plot for relative energy difference only for unique events",50,0,200,-100,100,"s"); //Compared to MC level jets
    TH1F *hBinnedEntriesMC = new TH1F("hBinnedEntriesMC", "Histogram to count entries in each bin",50,0,200); //Compared to MC level jets
    TH1F *hBinnedUniqueEntries = new TH1F("hBinnedUniqueEntries", "Histogram to count uniquely matched jets in each bin",50,0,200); //Compared to MC level jets
    TH1F *hUniqueEvents = new TH1F("hUniqueEvents", "Histogram to get fraction of unique events",50,0,200); //Compared to MC level jets
    
    TH1F *hConstituentsMatch = new TH1F("hConstituentsMatch","Constituents match in matched jets", 21, 0.0, 1.05);
    TH1F *hConstituentsMatchUnique = new TH1F("hConstituentsMatchUnique","Constituents match in uniquely matched jets", 21, 0.0, 1.05);
    TH1F *hConstituentsMatchNonunique = new TH1F("hConstituentsMatchNonunique","Constituents match in non-uniquely matched jets", 21, 0.0, 1.05);
    
    TH1F *hbhadron2particlejet = new TH1F("hbhadron2particlejet","Distribution of b-hadrons in particle jets", 21, 0.0, 1.05);
    TH1F *hbhadron2recojet = new TH1F("hbhadron2recojet","Distribution of b-hadrons in reco jets", 21, 0.0, 1.05);
    TH1F *hbhadron2recojetUnique = new TH1F("hbhadron2recojetUnique","Distribution of b-hadrons in reco jets", 21, 0.0, 1.05);
    TH1F *hbhadron2recojetNonunique = new TH1F("hbhadron2recojetNonunique","Distribution of b-hadrons in reco jets", 21, 0.0, 1.05);
    //How often do I find at least two reconstructed objects from b-hadron?
    TH1F *hbhadron2RP = new TH1F("hbhadron2RP","Distribution of number of reconstructed b-hadron decay productss", 31, -0.5, 30.5);
    //How many different different b-hadrons from which I have reconstructed particles (to see how many I loose)
    TH1F *hbhadronTypes = new TH1F("hbhadronTypes","Distribution of types b-hadrons from which I see reconstructed particles", 11, -0.5, 10.5);
    TH1F *hDeltaR = new TH1F("hDeltaR", "Distribuation of deltaR between b parton and closest jet", 100,0,10);
    TH1F *hAngle = new TH1F("hAngle", "Distribuation of angle between b parton and closest jet", 100,0,3);
    TH1F *hDeltaRall = new TH1F("hDeltaRall", "Distribuation of deltaR between b parton and all jets", 100,0,10);
    TH1F *hAngleall = new TH1F("hAngleall", "Distribuation of angle between b parton and all jets", 100,0,3);
    
    int Wmatched=0;
    int bmatched=0;
    int b2matched=0;
    int W2matched=0;
    int all=0;
    int itself=0;
    int count1=0;
    int count2=0;
    int countlepton=0;
    int uniquejetmatching=0;
    int MCuniquejetmatching=0;
    int uniguebtaggedparticleJetmatching=0;
    int nonuniguebtaggedparticleJetmatching=0;
    int uniguebtaggedrecoJetmatching=0;
    int nonuniguebtaggedrecoJetmatching=0;
    for (UInt_t j=0; j<chain->GetEntries();j++){
      //if (j>5) break;
      chain->GetEntry(j);
      
      //*********** TESTING WHETHER MC_parent1 IS ALWAYS FILLED FIRST ****************** 
      for (UInt_t i=0; i<MC_parent1->size();i++){
	Int_t parent1_test=MC_parent1->at(i);
	Int_t parent2_test=MC_parent2->at(i);
	if (parent1_test==-999 && parent2_test!=-999) cout << "Ladida not as expected" << endl;
	if (parent1_test==-999 && parent2_test==-999) count1++;
	if (parent1_test==-999) count2++;
      }
      if (count1!=count2) cout << "Ladida not as expected" << endl;
      
      //cout<< "Event no. " <<j << "---------------------------"<<endl;
      
      //*********************** LEPTONS ***********************************************
      //-------------------------------------------------------------------------------
      //Select leptons and get energy resolution from profile histogram
      //Profile histogram shows the spread of the difference in energy between MC and reconstructed particle as a function of energy 
      //Q: only highest energy leptons? or all? (Assume all for now)
      
      for (UInt_t i=0; i<RP_e->size();i++){
	if (abs(RPMC_pdg->at(i))==11 || abs(RPMC_pdg->at(i))==13) {
	  float leptonRP_e = RP_e->at(i); 
	  Int_t idx = RPMC_index->at(i);
	  float leptonMC_e = MC_e->at(idx); 
	  float lepton_diff = leptonRP_e - leptonMC_e; //Difference in energy between reconstructed particle and the MC particle it points to
	  hprofLepton->Fill(leptonRP_e,lepton_diff,1);
	}//end of if lepton loop
      }
      
      //***** Highest energy lepton found from the following loop ******************
      Float_t leptonmax=0;
      Int_t RPMCidx=0;
      for(UInt_t i=0; i < RP_e->size(); i++)  {
	if (abs(RPMC_pdg->at(i))==11 || abs(RPMC_pdg->at(i))==13) {
	  Float_t testmax=RP_e->at(i);
	  if (testmax>leptonmax){
	    leptonmax=testmax;
	    RPMCidx=i;
          }
	}
      }
      if (RPMCidx==0) continue; //Continue if no leptons found (WILL ALSO AFFECT JET STUDIES)
      Int_t idx = RPMC_index->at(RPMCidx);
      all++;
      
      
      //****** Matching highest energy lepton to hard process W boson or b quark ***********
      //------------------------------------------------------------------------------------
      bool running =true;
      vector<int> mothers;
      mothers.push_back(idx);
      while(running){
	if (mothers.size()==0) {
	  cout << "NO MOTHERS FOUND AT EVENT NUMBER " << j << endl;
	  break;
	}
	vector<int> newmothers;
	int counter=0;
	int bcounter=0;
	int Wcounter=0;
	for (auto index : mothers){
	  ///cout << "motherindex = " << index<< endl;
	  int parent1=MC_parent1->at(index);
	  int parent2=MC_parent2->at(index);
	  //if parent1==-999 then MC_parent2==-999
	  if (parent1==-999){
	    //cout << "Parent to lepton not found" << endl;
	    continue;
	  }
	  ///cout<< parent1 << " " << parent2 << endl;
	  newmothers.push_back(parent1);
	  if (parent2!=-999 && parent1!=parent2){
	    //cout << "Two different mothers found" << endl;
	    newmothers.push_back(parent2);
	  }
	  if (abs(MC_pdg->at(index))==24 && MC_status->at(index)==22) {
	    if (MC_pdg->at(index)>0 && MC_pdg->at(idx)>0) countlepton++;//cout << "Lepton comes from opposite W" << endl;
	    if (MC_pdg->at(index)<0 && MC_pdg->at(idx)<0) countlepton++;//cout << "Lepton comes from opposite W" << endl;
	    Wcounter++;
	    running=false;
	  }	
	  if (abs(MC_pdg->at(index))==5 && MC_status->at(index)==23) {
	    bcounter++;
	    running=false;
	  }
	}
	//for (auto index : mothers) cout << "Particle at index = " << index << " with pdg = " << MC_pdg->at(index) <<  " and status code = " << MC_status->at(index) << endl;
	mothers=newmothers;
	
	if (Wcounter>0 && bcounter==0) Wmatched++;
	if (Wcounter==0 && bcounter>0) bmatched++;
	/// cout << "New mothers = ";
	///for (auto index : mothers) cout << index << " "; 
	///cout << endl;
      }
      
      
      //-------------------------------------------------------------------------------     
      //*********************** JETS ************************************************** 
      //------------------------------------------------------------------------------- 

      //******* RP jets vs. hard process quarks ***************************************
      
      //Matching jets to hard process quarks to calculate energy difference and resolution
      //Also calculationg unique macthing 
      if (recojets_jetalgo_px->size()==0) continue;
      if (recojets_jetalgo_px->size()!=4) cout << "Jet size not equal to 4" << endl;
      unordered_set<float> Quark_pdg;  
      bool unique=true;
      for (UInt_t i=0; i<recojets_jetalgo_px->size();i++){
	Float_t pdgmatch;
	Float_t jetsmallestAngle=0;
	Float_t jetpx=recojets_jetalgo_px->at(i);
	Float_t jetpy=recojets_jetalgo_py->at(i);
	Float_t jetpz=recojets_jetalgo_pz->at(i);
	Float_t jete=recojets_jetalgo_e->at(i);
	Float_t jet_diff;
	Float_t jet_diffrel;
	for (UInt_t k=0; k<MC_pdg->size();k++){
	  if (MC_status->at(k)!=23) continue;
	  if (abs(MC_pdg->at(k))>6) continue;
	  Float_t MCpx=MC_px->at(k);
	  Float_t MCpy=MC_py->at(k);
	  Float_t MCpz=MC_pz->at(k);
	  
	  Float_t dot=jetpx*MCpx+jetpy*MCpy+jetpz*MCpz;
	  Float_t norm=sqrt((jetpx*jetpx+jetpy*jetpy+jetpz*jetpz)*(MCpx*MCpx+MCpy*MCpy+MCpz*MCpz));
	  Float_t angle=dot/norm;
	  if (jetsmallestAngle<angle) {
	    jetsmallestAngle=angle;
	    Float_t MCe=MC_e->at(k);
	    jet_diff=jete-MCe;
	    jet_diffrel=jet_diff/MCe;
	    pdgmatch=MC_pdg->at(k);
	    
	  }
	}
	hJetCosOpeningAngle->Fill(jetsmallestAngle);
	hprofJetAbs->Fill(jete,jet_diff,1);
	hprofJetRel->Fill(jete,jet_diffrel,1);
	unordered_set<float>::const_iterator got = Quark_pdg.find (pdgmatch);
	if ( got == Quark_pdg.end() ) {
	  Quark_pdg.insert(pdgmatch);  
	}
	else{
	  unique = false;
	}
	//**cout << pdgmatch << " " ;
      }
      //**if (unique) cout << "Unique" ;
      //**cout << endl;
      if (unique) uniquejetmatching++;
      
      
      //********** Trying to find hard process quark that each jet constituent points too ******
      //----------------------------------------------------------------------------------------
      /*
	vector<vector<float>> signalquarks(recojets_jetalgo_px->size());
	
	for (UInt_t i=0; i<RPrest_association->size() ; i++){
	if (jets1_association->at(i)<0) continue;
	int RPidx=RPrest_association->at(i);
	int MCidx=RPMC_index->at(RPidx);
	bool running =true;
	vector<int> jetmothers;
	jetmothers.push_back(MCidx);
	while(running){
	if (jetmothers.size()==0) {
	running=false;
	}
	vector<int> newjetmothers;
	for (auto index : jetmothers){
	if (MC_status->at(index)==23) {
	//cout << "Signal quark " << MC_pdg->at(index) << " found for jet # " << jets1_association->at(i) << " and RPMC index " << MCidx << " pdg = " << MC_pdg->at(MCidx) <<" px = " << MC_px->at(MCidx) << endl;
	signalquarks[jets1_association->at(i)].push_back(MC_pdg->at(index));
	//running=false;
	continue;
	}
	
	int parent1=MC_parent1->at(index);
	int parent2=MC_parent2->at(index);
	if (parent1==-999){
	//cout << "Parent to particle not found" << endl;                                                                                                                       
	continue;
	}
	newjetmothers.push_back(parent1);
	if (parent2!=-999 && parent1!=parent2){
	newjetmothers.push_back(parent2);
	}
	}
	jetmothers=newjetmothers;
	}
	}
      */
      
      //-------------------------------------------------------------------------------
      //******* Particle jets vs. reco jets *******************************************
      //-------------------------------------------------------------------------------
      
      //Matching jets to hard process quarks to calculate energy difference and resolution
      //Also calculationg unique macthing
      if (recojets_jetalgo_px->size()==0) continue;
      if (recojets_jetalgo_px->size()!=4) cout << "Jet size not equal to 4" << endl;
      vector<int> Jetmatch;
      vector<float> Jetenergy; //to save energy for the four jets and only fill hBinnedUniqueMatch if they are uniquely matched
      vector<float> Jetdiff_Unique; //to save relative energy diff for the four jets and only fill hprofJetRelMCUnique if they are uniquely matched
      unordered_set<int> nMCjet;
      bool MCunique=true;
      for (UInt_t i=0; i<recojets_jetalgo_px->size();i++){
	int jetmatch;
	Float_t jetsmallestAngle=0;
	Float_t jetpx=recojets_jetalgo_px->at(i);
	Float_t jetpy=recojets_jetalgo_py->at(i);
	Float_t jetpz=recojets_jetalgo_pz->at(i);
	Float_t jete=recojets_jetalgo_e->at(i);
	Float_t jet_diff;
	Float_t jet_diffrel;
	for (UInt_t k=0; k<particlejets_jetalgo_px->size();k++){
	  Float_t MCpx=particlejets_jetalgo_px->at(k);
	  Float_t MCpy=particlejets_jetalgo_py->at(k);
	  Float_t MCpz=particlejets_jetalgo_pz->at(k);
	  
	  Float_t dot=jetpx*MCpx+jetpy*MCpy+jetpz*MCpz;
	  Float_t norm=sqrt((jetpx*jetpx+jetpy*jetpy+jetpz*jetpz)*(MCpx*MCpx+MCpy*MCpy+MCpz*MCpz));
	  Float_t angle=dot/norm;
	  if (jetsmallestAngle<angle) {
	    jetsmallestAngle=angle;
	    Float_t MCe=MC_e->at(k);
	    jet_diff=jete-MCe;
	    jet_diffrel=jet_diff/MCe;
	    jetmatch=k;
	  }
	}
	Jetmatch.push_back(jetmatch);
	Jetenergy.push_back(jete);
	Jetdiff_Unique.push_back(jet_diffrel);
	hJetCosOpeningAngleMC->Fill(jetsmallestAngle);
	hprofJetAbsMC->Fill(jete,jet_diff,1);
	hprofJetRelMC->Fill(jete,jet_diffrel,1);
	//Filling histogram with same binning as profile plot to get number of entries in each bin
	//Used for calculating the error on error (assuming Gaussian): sqrt(2*n/std^4)
	hBinnedEntriesMC->Fill(jete); 
	unordered_set<int>::const_iterator got = nMCjet.find (jetmatch);
	if ( got == nMCjet.end() ) {
	  nMCjet.insert(jetmatch);
	}
	else{
	  MCunique = false;
	}
      }
      
      //**if (MCunique) cout << "MCunique" ;
      //**cout << endl;
      if (MCunique) {
	MCuniquejetmatching++;
	for (unsigned i = 0; i < Jetenergy.size(); i++) {
	  float jete=Jetenergy[i];
	  float jet_diffrel=Jetdiff_Unique[i];
	  hBinnedUniqueEntries->Fill(jete);
	  hUniqueEvents->Fill(jete);
	  hprofJetRelMCUnique->Fill(jete,jet_diffrel,1);
	}
      }
      vector<unordered_set<int>> recojet2MCsets;
      vector<unordered_set<int>> particlejet2MCsets;
      
      for (auto constituents : *recojetconstituents_jetalgo){
	unordered_set<int> recojet2MC;
	for (unsigned i = 0; i < constituents.size(); i++) {
	  int recojet2RPindex=RPrest_association->at(constituents[i]);
	  int recojet2MCindex=RPMC_index->at(recojet2RPindex);
	  unordered_set<int>::const_iterator got = recojet2MC.find (recojet2MCindex);
	  if ( got == recojet2MC.end() ) recojet2MC.insert(recojet2MCindex);
	}
	recojet2MCsets.push_back(recojet2MC);
      }
      
      for (auto constituents : *particlejetconstituents_jetalgo){
	unordered_set<int> particlejet2MC;
	for (unsigned i = 0; i < constituents.size(); i++) {
	  int particlejet2MCindex=MCfinal_association->at(constituents[i]);
	  unordered_set<int>::const_iterator got = particlejet2MC.find (particlejet2MCindex);
	  if ( got == particlejet2MC.end() ) particlejet2MC.insert(particlejet2MCindex);
	}
	particlejet2MCsets.push_back(particlejet2MC);
      }
      
      //cout << MCunique << ": " ;
      for (unsigned i = 0; i < recojet2MCsets.size(); i++){
	unordered_set<int> recoset     = recojet2MCsets[i];
	unordered_set<int> particleset = particlejet2MCsets[Jetmatch[i]];
	
	unsigned int found_jets = 0;
	for (const auto& reco: recoset) {
	  unordered_set<int>::const_iterator got = particleset.find (reco);
	  if ( got != particleset.end() ) found_jets++; 
	}
	//cout << float(found_jets)/float(recoset.size()) << ", " ;
	hConstituentsMatch->Fill(float(found_jets)/float(recoset.size()));
	if (MCunique) hConstituentsMatchUnique->Fill(float(found_jets)/float(recoset.size()));
	else hConstituentsMatchNonunique->Fill(float(found_jets)/float(recoset.size()));
	
      }
      //cout << endl;
      
      //-------------------------------------------------------------------------------
      //******* b-hadron vs. particle jets ********************************************
      //-------------------------------------------------------------------------------
      
      /*int countingjets=0;
	for (auto constituents : *particlejetconstituents_jetalgo){
	cout << "Jet no " << countingjets << ": ";
	countingjets++;
	for (unsigned i = 0; i < constituents.size(); i++) {
	int particlejet2MCindex=MCfinal_association->at(constituents[i]);
	cout << MC_pdg->at(particlejet2MCindex) << " ";
	}
	cout << endl;
	}*/
      int countingparticlejets=0;
      unordered_multiset<int> btaggedParticleJets;
      unordered_set<float> allparticlebhadrons;
      vector<unordered_multiset<float>> allparticlebtags;
      for (auto constituents : *particlejetconstituents_jetalgo){
	//cout << "------------------------------------------" << endl;
	//cout << "Jet no " << countingparticlejets << ": " << endl;
	unordered_set<float> bhadrons;
	unordered_multiset<float> btags;
	for (unsigned i = 0; i < constituents.size(); i++) { 
	  //For jet constituent, find corresponding particle in MC collection to check history for
	  int particlejet2MCindex=MCfinal_association->at(constituents[i]);
	  bool running =true;
	  vector<int> jetmothers;
	  jetmothers.push_back(particlejet2MCindex);
	  //cout << "Jet constituent with pdg " << MC_pdg->at(particlejet2MCindex) << ": " <<endl;
	  while(running){
	    if (jetmothers.size()==0) {
	      running=false;
	    }
	    vector<int> newjetmothers;
	    for (auto index : jetmothers){
	      //cout << index <<": " << MC_pdg->at(index) << " "<< MC_status->at(index) << ", " ;
	      if (MC_status->at(index)>10) continue;
	      
	      int parent1=MC_parent1->at(index);
	      int parent2=MC_parent2->at(index);
	      
	      //b-tagging by testing if b-baryon (div(1000)==5) or b-meson (div(100)==5)
	      bool btagged=false;
	      div_t divbaryon= div(MC_pdg->at(index),1000);
	      if (divbaryon.quot==0){
		div_t divmeson= div(MC_pdg->at(index),100);
		if (abs(divmeson.quot)==5) btagged=true;
	      }
	      if (abs(divbaryon.quot)==5) btagged=true;
	      
	      if (btagged && MC_status->at(parent1)>10) {
		unordered_set<float>::const_iterator got = bhadrons.find (MC_pdg->at(index));
		if ( got == bhadrons.end() ) bhadrons.insert(MC_pdg->at(index));
		btags.insert(MC_pdg->at(index));
		continue;
	      }
	      
	      if (parent1==-999){
		//cout << "Parent to particle not found" << endl;
		continue;
	    }
	      newjetmothers.push_back(parent1);
	      if (parent2!=-999 && parent1!=parent2){
		newjetmothers.push_back(parent2);
	      }
	    }
	    jetmothers=newjetmothers;
	    //cout << endl;
	  }
	}
	/*      cout << bhadrons.size() << " " << btags.size() << endl;
		for (auto bpdg : bhadrons) {
		cout << bpdg << ": " << btags.count(bpdg) << endl;
		}*/
	btaggedParticleJets.insert(bhadrons.size());
	for (auto bhadron : bhadrons) allparticlebhadrons.insert(bhadron);
	allparticlebtags.push_back(btags);
	countingparticlejets++;
      }
      unordered_multiset<int> UniquebtaggedJets={1,1,0,0};
      if (btaggedParticleJets==UniquebtaggedJets) uniguebtaggedparticleJetmatching++;//cout << "*********Uniquely b-tagged jets" << endl;
      else nonuniguebtaggedparticleJetmatching++; //cout << "*********Not uniquely b-tagged jets" << endl;
      
      float btagparticlefraction;
      for (auto bhadron : allparticlebhadrons) {
	vector<int> count;
	for (auto btags  : allparticlebtags) {
	  count.push_back(btags.count(bhadron));
	}
	btagparticlefraction=float(*max_element(count.begin(), count.end()))/float(accumulate(count.begin(), count.end(), 0));
	hbhadron2particlejet->Fill(btagparticlefraction);
      }
      

      //-------------------------------------------------------------------------------
      //******* b-hadron vs. reco jets ************************************************
      //-------------------------------------------------------------------------------
      unordered_multiset<int> btaggedRecoJets;
      unordered_set<float> allrecobhadrons;
      vector<unordered_multiset<float>> allrecobtags;
      for (auto constituents : *recojetconstituents_jetalgo){
	unordered_set<float> bhadrons;
	unordered_multiset<float> btags;
	for (unsigned i = 0; i < constituents.size(); i++) { 
	  //For jet constituent, find corresponding particle in MC collection to check history for
	  int recojet2RPindex=RPrest_association->at(constituents[i]);
	  int recojet2MCindex=RPMC_index->at(recojet2RPindex);
	  bool running =true;
	  vector<int> jetmothers;
	  jetmothers.push_back(recojet2MCindex);
	  while(running){
	    if (jetmothers.size()==0) {
	      running=false;
	    }
	    vector<int> newjetmothers;
	    for (auto index : jetmothers){
	      if (MC_status->at(index)>10) continue;
	      int parent1=MC_parent1->at(index);
	      int parent2=MC_parent2->at(index);
	      
	      //b-tagging by testing if b-baryon (div(1000)==5) or b-meson (div(100)==5)
	      bool btagged=false;
	      div_t divbaryon= div(MC_pdg->at(index),1000);
	      if (divbaryon.quot==0){
		div_t divmeson= div(MC_pdg->at(index),100);
		if (abs(divmeson.quot)==5) btagged=true;
	      }
	      if (abs(divbaryon.quot)==5) btagged=true;
	      
	      if (btagged && MC_status->at(parent1)>10) {
		unordered_set<float>::const_iterator got = bhadrons.find (MC_pdg->at(index));
		if ( got == bhadrons.end() ) bhadrons.insert(MC_pdg->at(index));
		btags.insert(MC_pdg->at(index));
		continue;
	      }
	      
	      if (parent1==-999){
		//cout << "Parent to particle not found" << endl;
		continue;
	      }
	      newjetmothers.push_back(parent1);
	      if (parent2!=-999 && parent1!=parent2){
		newjetmothers.push_back(parent2);
	      }
	    }
	    jetmothers=newjetmothers;
	  }
	}
	btaggedRecoJets.insert(bhadrons.size());
	for (auto bhadron : bhadrons) allrecobhadrons.insert(bhadron);
	allrecobtags.push_back(btags);
      }
      if (btaggedRecoJets==UniquebtaggedJets) uniguebtaggedrecoJetmatching++;//cout << "*********Uniquely b-tagged jets" << endl;
      else nonuniguebtaggedrecoJetmatching++; //cout << "*********Not uniquely b-tagged jets" << endl;
      
      float btagrecofraction;
      for (auto bhadron : allrecobhadrons) {
	vector<int> count;
	for (auto btags  : allrecobtags) {
	  count.push_back(btags.count(bhadron));
	}
	btagrecofraction=float(*max_element(count.begin(), count.end()))/float(accumulate(count.begin(), count.end(), 0));
	hbhadron2recojet->Fill(btagrecofraction);
	hbhadron2RP->Fill(float(accumulate(count.begin(), count.end(), 0))); //Count number of bhadron decay products that are RECONSTRUCTED
	if (MCunique) hbhadron2recojetUnique->Fill(btagrecofraction);
	else hbhadron2recojetNonunique->Fill(btagrecofraction);
      }
      
      hbhadronTypes->Fill(allrecobhadrons.size()); //Count number of bhadron that leaves RECONSTRUCTED decay products

      //-------------------------------------------------------------------------------
      //******* b-parton vs. reco jets ********************************************
      //-------------------------------------------------------------------------------
      int loopcount =0;
      for (size_t i = 0; i < MC_status->size(); ++i) {
	if (MC_status->at(i)>80 || MC_status->at(i)<70) continue;
	if (abs(MC_pdg->at(i)) != 5) continue;
	ROOT::Math::PxPyPzMVector partonlv(MC_px->at(i), MC_py->at(i), MC_pz->at(i), MC_mass->at(i));
	float deltaRmin=100;
	float anglemin=100;
	loopcount++;
	//cout << loopcount << " :    ";
	for (size_t j = 0; j < recojets_jetalgo_px->size(); ++j) {
	  ROOT::Math::PxPyPzEVector jetlv(recojets_jetalgo_px->at(j), recojets_jetalgo_py->at(j), recojets_jetalgo_pz->at(j), recojets_jetalgo_e->at(j));
	  float dEta = partonlv.Eta() - jetlv.eta();
	  float dPhi = partonlv.Phi() - jetlv.phi();
	  float deltaR = sqrt(dEta*dEta+dPhi*dPhi);
	  //cout << j << ": " << deltaR << ";    ";
	  hDeltaRall->Fill(deltaR);
	  if (deltaR<deltaRmin) deltaRmin=deltaR;

	  Float_t dot = jetlv.px()*partonlv.px()+jetlv.py()*partonlv.py()+jetlv.pz()*partonlv.pz();
	  Float_t lenSq1 = jetlv.px()*jetlv.px()+jetlv.py()*jetlv.py()+jetlv.pz()*jetlv.pz();
	  Float_t lenSq2 = partonlv.px()*partonlv.px()+partonlv.py()*partonlv.py()+partonlv.pz()*partonlv.pz();
	  Float_t norm = sqrt(lenSq1*lenSq2);
	  Float_t angle = acos(dot/norm);
	  //cout << j << ": " << angle << ";    ";
	  hAngleall->Fill(angle);
	  if (angle<anglemin) anglemin=angle;
	}
	//cout << endl;
	hDeltaR->Fill(deltaRmin);
	hAngle->Fill(anglemin);
      }
      

    }
    
    //**************** COUT STATEMENTS, WRITING TO FILE AND MAKING PLOTS *****************************
    cout << "------------------- "+jetalgo+" -------------------" << endl;

    cout << "Number of wrongly matched leptons = " << countlepton << endl;
    
    cout << "Number of highest energy leptons matched to W boson: " << Wmatched << "/" << all << "=" << float(Wmatched)/float(all) <<endl;   
    cout << "Number of highest energy leptons matched to b quark: " <<bmatched << "/" << all << "=" << float(bmatched)/float(all) <<endl;   
  
    cout << "Unique matching between quark and reco jet: " << uniquejetmatching << "/" << all << "=" << float(uniquejetmatching)/float(all) <<endl;
    cout << "Unique matching between particle and reco jet: " << MCuniquejetmatching << "/" << all << "=" << float(MCuniquejetmatching)/float(all) <<endl;
    
    cout << "Events with full separation of b-hadron decay products into particle jets: " << uniguebtaggedparticleJetmatching  << "/" << all << "=" << float(uniguebtaggedparticleJetmatching)/float(all) <<endl;
    cout << "Events with full separation of b-hadron decay products into reco jets: " << uniguebtaggedrecoJetmatching << "/" << all << "=" << float(uniguebtaggedrecoJetmatching)/float(all) <<endl; 
    
    hUniqueEvents->Divide(hBinnedEntriesMC);

    hprofLepton->SetDirectory(outFile);
    hJetCosOpeningAngle->SetDirectory(outFile);
    hprofJetRel->SetDirectory(outFile);
    hprofJetAbs->SetDirectory(outFile);
    hJetCosOpeningAngleMC->SetDirectory(outFile);
    hprofJetRelMC->SetDirectory(outFile);
    hprofJetAbsMC->SetDirectory(outFile);
    hprofJetRelMCUnique->SetDirectory(outFile);
    hBinnedEntriesMC->SetDirectory(outFile);
    hBinnedUniqueEntries->SetDirectory(outFile);
    hUniqueEvents->SetDirectory(outFile);
    hConstituentsMatch->SetDirectory(outFile);  
    hConstituentsMatchUnique->SetDirectory(outFile);
    hConstituentsMatchNonunique->SetDirectory(outFile);
    hbhadron2particlejet->SetDirectory(outFile);
    hbhadron2recojet->SetDirectory(outFile);
    hbhadron2recojetUnique->SetDirectory(outFile);
    hbhadron2recojetNonunique->SetDirectory(outFile);
    hbhadronTypes->SetDirectory(outFile);
    hbhadron2RP->SetDirectory(outFile);
    hDeltaR->SetDirectory(outFile);
    hAngle->SetDirectory(outFile);
    hDeltaRall->SetDirectory(outFile);
    hAngleall->SetDirectory(outFile);
    // kT, kT, eeantikt, eeCam
    
    TCanvas* c1 = new TCanvas("cst","stacked hists",10,10,700,500);
    hprofJetAbs->GetXaxis()->SetTitle("Energy [GeV]");
    hprofJetAbs->GetYaxis()->SetTitle("Mean on absolute difference [GeV]");
    hprofJetAbs->GetXaxis()->SetTitleSize(0.05);
    hprofJetAbs->GetYaxis()->SetTitleSize(0.05);
    hprofJetAbs->GetXaxis()->SetTitleOffset(0.9);
    hprofJetAbs->GetYaxis()->SetTitleOffset(0.9);
    hprofJetAbs->SetStats(0);
    hprofJetAbs->SetMarkerSize(0.7);
    hprofJetAbs->SetMarkerStyle(9);
    hprofJetAbs->Draw("P");
    c1->SaveAs("hprofJetAbs_"+jetalgo+".png");
    delete c1;
    
    
    TCanvas* c2 = new TCanvas("cst","stacked hists",10,10,700,500);
    c2->SetLogy();
    hJetCosOpeningAngle->GetXaxis()->SetTitle("Cos(matching angle)");
    hJetCosOpeningAngle->GetYaxis()->SetTitle("Event rate per 0.01");
    hJetCosOpeningAngle->GetXaxis()->SetTitleSize(0.05);
    hJetCosOpeningAngle->GetYaxis()->SetTitleSize(0.05);
    hJetCosOpeningAngle->GetXaxis()->SetTitleOffset(0.9);
    hJetCosOpeningAngle->GetYaxis()->SetTitleOffset(0.9);
    hJetCosOpeningAngle->SetStats(0);
    hJetCosOpeningAngle->SetMarkerSize(16);
    hJetCosOpeningAngle->Draw();
    c2->SaveAs("hJetCosOpeningAngle_"+jetalgo+".png");
    delete c2;
    
    
    TCanvas* c3 = new TCanvas("cst","stacked hists",10,10,700,500);
    TGraph* gJetRel = new TGraph();
    gJetRel->GetXaxis()->SetTitle("Energy [GeV]");
    gJetRel->GetYaxis()->SetTitle("Standard deviation [GeV]");
    gJetRel->GetXaxis()->SetTitleSize(0.05);
    gJetRel->GetYaxis()->SetTitleSize(0.05);
    gJetRel->GetXaxis()->SetTitleOffset(0.9);
    gJetRel->GetYaxis()->SetTitleOffset(0.9);
    gJetRel->SetTitle("Energy resolution of jets");
    //gJetRel->SetTitleSize(0.1);                                                                                                                                  
    for (int i=0 ; i<hprofJetRel->GetNbinsX();i++){
      float Std = hprofJetRel->GetBinError(i);
      if (Std==0) continue;
      float energy = hprofJetRel->GetXaxis()->GetBinCenter(i);
      gJetRel->SetPoint(i,energy,Std);
    }
    gJetRel->SetMarkerColor(kBlack);
    gJetRel->SetMarkerSize(0.7);
    gJetRel->SetMarkerStyle(8);
    gJetRel->Draw("AP");
    c3->SaveAs("gJetRel_"+jetalgo+".png");
    delete c3;

    //************ Resolution plot for RP vs. MC jets **********************
    //----------------------------------------------------------------------
    
    TCanvas* c4 = new TCanvas("cst","stacked hists",10,10,700,500);
    //TGraph* gJetRelMC = new TGraph();
    TGraphErrors* gJetRelMC = new TGraphErrors();
    gJetRelMC->GetXaxis()->SetTitle("Energy [GeV]");
    gJetRelMC->GetYaxis()->SetTitle("Standard deviation [GeV]");
    gJetRelMC->GetXaxis()->SetTitleSize(0.05);
    gJetRelMC->GetYaxis()->SetTitleSize(0.05);
    gJetRelMC->GetXaxis()->SetTitleOffset(0.9);
    gJetRelMC->GetYaxis()->SetTitleOffset(0.9);
    gJetRelMC->SetTitle("Energy resolution of RP vs. MC jets");
    //gJetRelMC->SetTitleSize(0.1);                                                                                                                                                    
    for (int i=0 ; i<hprofJetRelMC->GetNbinsX();i++){
      double n = hBinnedEntriesMC->GetBinContent(i);
      if (n<2) continue;
      double Std = hprofJetRelMC->GetBinError(i);
      double energy = hprofJetRelMC->GetXaxis()->GetBinCenter(i);
      //double error = sqrt(2/n)*(Std*Std);
      //double error = sqrt((3-(n-3)/(n-1))/n)*(Std*Std);
      double error = Std/sqrt(2*n);
      //cout << Std << "^2*sqrt(2/" << n << ")=" << error << endl;
      //cout << i <<  ": " << n << ", " << energy << ", " << Std << ", " << error << endl; 
      gJetRelMC->SetPoint(i,energy,Std);
      gJetRelMC->SetPointError(i,0,error);
      
    }
    gJetRelMC->SetMarkerColor(kBlack);
    gJetRelMC->SetMarkerSize(0.7);
    gJetRelMC->SetMarkerStyle(8);
    gJetRelMC->Draw("AP");
    c4->SaveAs("gJetRelMC_"+jetalgo+".png");
    delete c4;
    
    
    
    outFile->Write();
    outFile->Close();
    
    delete outFile;
    
  }

  return 0;
    
}
