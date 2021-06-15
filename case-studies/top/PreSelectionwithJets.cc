//The order of the selection cuts are based on a significance-optimized selection cut strategy

//To run: g++ -o PreSelectionwithJets -std=gnu++17 -I /Users/Julie/ROOT/include PreSelectionwithJets.cc `root-config --cflags --libs`
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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

using namespace std;

TFile *outFile;

struct descending
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};

float invMass(vector< array<float,4> > ParticleSet) {
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  for (auto P : ParticleSet){
    E+=P[0];
    px+=P[1];
    py+=P[2];
    pz+=P[3];
  }
  float invMass_squared = pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2); 
  return sqrt(invMass_squared);
}

int main()
{
  std::vector<TString> process;
  /*process.emplace_back(TString("ttbar_semi_ecm365"));
  process.emplace_back(TString("mumu_ecm365"));
  process.emplace_back(TString("qqbar_ecm365"));
  process.emplace_back(TString("bbbar_ecm365"));
  process.emplace_back(TString("gmZ_ecm365"));
  process.emplace_back(TString("WW_ecm365"));
  process.emplace_back(TString("ZZ_ecm365"));
  process.emplace_back(TString("ZWW_ecm365_v2"));
  process.emplace_back(TString("ZZZ_ecm365_v2"));
  process.emplace_back(TString("singleTop_ecm365_v3"));*/

  process.emplace_back(TString("ttbar_semi_pos_el"));
  process.emplace_back(TString("ttbar_semi_neg_el"));
  process.emplace_back(TString("ttbar_semi_pos_mu"));
  process.emplace_back(TString("ttbar_semi_neg_mu"));
  process.emplace_back(TString("ttbar_semi_pos_ta"));
  process.emplace_back(TString("ttbar_semi_neg_ta"));
  process.emplace_back(TString("mumu"));
  process.emplace_back(TString("tautau"));
  process.emplace_back(TString("qqbar"));
  process.emplace_back(TString("bbbar"));
  process.emplace_back(TString("gmZ"));
  process.emplace_back(TString("WW"));
  process.emplace_back(TString("ZZ"));
  process.emplace_back(TString("ZH"));
  process.emplace_back(TString("ZWW"));
  process.emplace_back(TString("ZZZ"));
  process.emplace_back(TString("singletop"));

  
  std::vector<float> sigma; //cross sections in pb                               
  //sigma.push_back(float(0.4868)); //Sigma for total ttbar                                              
  //sigma.push_back(float(0.1315)); //Sigma for semi-leptonic ttbar only with l={e,mu}
  //sigma.push_back(float(0.1933)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.1933/6)); //Sigma for semi-leptonic ttbar only with l={e,mu,tau}  
  sigma.push_back(float(0.7942));
  sigma.push_back(float(0.7937));
  sigma.push_back(float(4.143));
  sigma.push_back(float(0.7448));
  sigma.push_back(float(3.386));
  sigma.push_back(float(10.72));
  sigma.push_back(float(0.6428));
  sigma.push_back(float(0.1173));
  sigma.push_back(float(0.01591));
  sigma.push_back(float(0.0007633));
  sigma.push_back(float(0.003337));

  std::vector<TString> modepath;
  //modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_semi_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_electron_pos_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_electron_neg_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_muon_pos_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_muon_neg_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_tau_pos_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ttbar_tau_neg_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_mumu_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_tautau_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_qqbar_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_bbbar_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_gmZ_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_WW_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_ZZ_ecm365.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_ZH_ecm365_v1.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_ZWW_ecm365_v2.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov/p8_ee_ZZZ_ecm365_v2.root"));
  modepath.emplace_back(TString("ttbar_IDEAtrkCov_new/p8_ee_singleTop_no_t_ecm365_run_03.root"));
  
  Float_t tol = 1.0e-16;
  Float_t MW = 80.4; //mass of W boson
  Float_t Mtop = 173; //mass of top quark
  
  int ix=0;
  vector<float> Passed;
  vector<float> Sanity;
  vector<vector<UInt_t>> cutflow;
  vector<vector<UInt_t>> flowerror;
  vector<float> scalefactors;
  for(auto& mode : process) {

    outFile = new TFile("EventSelection/p8_eventselection_"+mode+".root", "RECREATE");

    TChain* chain = new TChain("events");
    chain->Add("/eos/user/j/jutornda/FCCee/"+modepath[ix]);

    //Defining relevant branches in the file                                     
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

    vector<float> *jets_px=0;
    vector<float> *jets_py=0;
    vector<float> *jets_pz=0;
    vector<float> *jets_e=0;
    vector<vector<int>> *jetconstituents=0;
    double jets_dmin;
    vector<int> *jets_btag=0;
    vector<int> *RPrest_association=0; //To get the correct association between jet constituents and RP
    vector<float> *jetsignificance=0;

    float event_thrust_val=0;
    float eventRest_thrust_val=0;

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
    chain->SetBranchAddress("MET_e", & MET_e);

    chain->SetBranchAddress("RPMC_pdg", & RPMC_pdg);
    chain->SetBranchAddress("RPMC_index", & RPMC_index);

    chain->SetBranchAddress("jets_px", &jets_px);
    chain->SetBranchAddress("jets_py", &jets_py);
    chain->SetBranchAddress("jets_pz", &jets_pz);
    chain->SetBranchAddress("jets_e", &jets_e);
    chain->SetBranchAddress("jetconstituents", &jetconstituents);
    chain->SetBranchAddress("jets_dmin", &jets_dmin);
    chain->SetBranchAddress("jets_btag", &jets_btag);
    chain->SetBranchAddress("RPrest_association", &RPrest_association);
    chain->SetBranchAddress("jetsignificance", &jetsignificance);
    
    chain->SetBranchAddress("EVT_thrust_val", & event_thrust_val);
    chain->SetBranchAddress("EVTrest_thrust_val", & eventRest_thrust_val);
    
    
    Int_t n=1000;//1000
    TH1F *hHighestEnergyLepton = new TH1F("hHighestEnergyLepton", "Highest Energy Lepton", 70*n, 0.0, 140.0);
    TH1F *h2ndHighestEnergyLepton = new TH1F("h2ndHighestEnergyLepton", "2nd Highest Energy Lepton", 70*n, 0.0, 140.0);//
    TH1F *hLeptonMomentum = new TH1F("hLeptonMomentum", "Lepton Momentum", 70*n, 0.0, 140.0);
    TH1F *hLeptonMomentumRest = new TH1F("hLeptonMomentumRest", "Lepton Momentum excluding highest energy lepton", 70*n, 0.0, 140.0);   //
    TH1F *hDiffLepton = new TH1F("hDiffLepton", "Momentum difference", 70*n, 0.0, 140.0);
    TH1F *hMissingMomentum = new TH1F("hMissingMomentum", "Missing Momentum", 200*n, 0.0, 400.0);
    TH1F *hMissingMomentum2 = new TH1F("hMissingMomentum2", "Missing Momentum", 200*n, 0.0, 400.0); 
    TH1F *hInvariantMassLepton = new TH1F("hInvariantMassLepton", "Invariant mass of leptonic pair", 100*n, 0.0, 200.0);
    TH1F *hInvariantMassHighestLeptons = new TH1F("hInvariantMassHighestLeptons", "Invariant mass of highest and second highest lepton", 100*n, 0.0, 200.0); //
    TH1F *hInvariantMassRest = new TH1F("hInvariantMassRest", "Invariant mass of the rest", 100*n, 0.0, 400.0);
    TH1F *hInvariantMassAllRP = new TH1F("hInvariantMassAllRP", "Invariant mass of the allRP", 100*n, 0.0, 400.0);
    TH1F *hThrustAllRP = new TH1F("hThrustAllRP", "Thrust of entire event", 100*n, 0.0, 1.0);
    TH1F *hThrustRest = new TH1F("hThrustRest","Thrust of event excluding highest energy lepton", 100*n, 0.0, 1.0);

    TH1F *hJetMass = new TH1F("hJetMass", "Mass of jet", 120*n,0.0,60);
    TH1F *hJetEnergy = new TH1F("hJetEnergy", "Energy of jet", 100*n,0.0,200);
    TH1F *hbtagNumber = new TH1F("hbtagNumber", "Number of b-tagged jets per event", 5, 0.0, 5.0);
    TH1F *hSignificance = new TH1F("hSignificance", "Significance distribution", 100*n,0,0.1);
    TH1F *hdmin = new TH1F("hdmin", "dmin for each event", 250*n,0,500);
    TH1F *hDeltaWdecay = new TH1F("hDeltaWdecay", "Delta(m_i-m_W/2)", 25*n,0,50);
    TH1F *hDeltaW = new TH1F("hDeltaW", "Delta(m_ij-m_W)", 50*n,0,100);
    TH1F *hDeltaTop = new TH1F("hDeltaTop", "Delta(m_ijk-m_t)", 100*n,0,200);
    TH1F *hDeltaToplep = new TH1F("hDeltaToplep", "Delta(m_{l,nu,i}-m_t)", 100*n,0,200);
      
    UInt_t missinglepton=0;
    UInt_t zerojets=0;
    UInt_t eventcount=0;

    int interval_1=5000;
    int interval_2=50000;

    UInt_t count=0;
    vector<UInt_t> cutcount(10,0);   
    cout << "--------------------------- "+mode+" loop -----------------------------" << endl;
    cout << "Number of entries = " << chain->GetEntries() << endl;
    for (UInt_t j=0; j<chain->GetEntries();j++){
      chain->GetEntry(j);
      count++;
      cutcount[0]++;
      //if (count>1000) break; //for testing
      if (count>200000 && mode!=TString("WW")) break;
      if (count>1000000) break;
      if (count%interval_1 == 0) cout << "."; cout.flush();
      if (count%interval_2 == 0) cout << " entries " << count << endl;

      //cout<< "Event no. " <<j << "---------------------------"<<endl;
      UInt_t nLepton=0;
      vector<Int_t> Lepton_idx;
      vector<Float_t> Lepton_p;
      vector<Float_t> Lepton_px;
      vector<Float_t> Lepton_py;
      vector<Float_t> Lepton_pz;
      vector<Float_t> Lepton_pdg;
      Float_t ME_p=0;
      for (UInt_t i=0; i<RP_p->size();i++){
	//---------- MISSING MOMENTUM
	ME_p+=RP_p->at(i);
	//---------- LEPTON MOMENTUM
	if (abs(RPMC_pdg->at(i))==11 || abs(RPMC_pdg->at(i))==13) {
          Int_t idx = RPMC_index->at(i);
	  Lepton_idx.push_back(i);
	  Lepton_p.push_back(RP_p->at(i));
	  Lepton_px.push_back(RP_px->at(i));
	  Lepton_py.push_back(RP_py->at(i));
	  Lepton_pz.push_back(RP_pz->at(i));
	  Lepton_pdg.push_back(RPMC_pdg->at(i));
	  nLepton++;                         
	}
      }

      // 1. cut ------- SANITY CHECK THAT THERE IS AT LEAST ONE RECONSTRUCTED LEPTON IN THE EVENT
      if (nLepton==0) {
	missinglepton+=1;
	//cout << "No leptons found at " << j << endl;
	continue;     
	}
      cutcount[1]++;
      // 2. cut ------- UPPER CUT ON THRUST
      if (event_thrust_val>0.85) continue;
      cutcount[2]++;
      
      //---------- Sorting LEPTON MOMENTUM 
      vector<int> sortedLepton(Lepton_p.size()); //Saving index for sorting the momentum vector in descinding order to also change for other properties
      size_t n(0);
      generate(begin(sortedLepton), end(sortedLepton), [&]{ return n++; });
      
      sort(  begin(sortedLepton), 
	     end(sortedLepton),
	     [&](int i1, int i2) { return Lepton_p[i1] > Lepton_p[i2]; } );
      
      int iL=sortedLepton[0]; //index for highest energy lepton

      // ---------- INVARIANT MASS for the rest of RP collection
      vector < array<float,4> > rest;
      vector < array<float,4> > allRP;
      for (UInt_t i=0; i<RP_p->size();i++){
	allRP.push_back({RP_e->at(i), RP_px->at(i), RP_py->at(i), RP_pz->at(i)});
	if (i==Lepton_idx[iL]) continue;
	rest.push_back({RP_e->at(i), RP_px->at(i), RP_py->at(i), RP_pz->at(i)});
      }
      float invMass_rest=invMass(rest);
      // 3. cut ------- LOWER CUT ON INVARIANT MASS OF EVENT EXCLUDING HIGHEST ENERGY LEPTON
      if (invMass_rest < 160) continue;
      cutcount[3]++;
      //if (invMass_rest > 300) continue;

      // 4. cut ---------- LOWER CUT ON INVARIANT MASS for lepton and missing ET
      vector < array<float,4> > leptonset;
      leptonset.push_back({Lepton_p[iL], Lepton_px[iL], Lepton_py[iL], Lepton_pz[iL]});
      leptonset.push_back({MET_e->at(0), MET_px->at(0), MET_py->at(0), MET_pz->at(0)});
      float invMass_lepton=invMass(leptonset);
      if (invMass_lepton<50) continue;
      //if (invMass_lepton<60 || invMass_lepton>110) continue;
      cutcount[4]++;
      // 5. cut ------- UPPER CUT ON HIGHEST ENERGY LEPTON
      if (Lepton_p[iL]>100) continue;
      cutcount[5]++;
      // 6. cut ------- LOWER CUT ON HIGHEST ENERGY LEPTON  
      if (Lepton_p[iL]<15) continue;
      cutcount[6]++;

      if (nLepton>1){ 
	int i2L=sortedLepton[1]; //index for second highest energy lepton
	//7. cut --- UPPER CUT ON 2ND HIGHEST ENERGY LEPTON/LEPTON MOMENTUM EXCLUDING HIGHEST ENERGY LEPTON  
	//Because the lepton momentums are sorted, it's the same thing
	if (Lepton_p[i2L]>40) continue;
      }
      cutcount[7]++;

      // Cut -------- SANITY CHECK THAT THERE IS ARE RECONSTRUCTED JETS
      if (jets_px->size()==0) {
	zerojets++;
	continue;
      }
      if (jets_px->size()!=4) cout << "Jet size not equal to 4" << endl;
      cutcount[8]++;
      // Cut -------- B-TAGGING
      int Nbtags=0;
      for (UInt_t i=0; i<jets_btag->size();i++){
        if (jets_btag->at(i)==1) Nbtags++;
      }
      if (Nbtags==0) continue;
      cutcount[9]++;


      if (nLepton>1){
	int i2L=sortedLepton[1]; //index for second highest energy lepton
	
	//------ MOMENTUM DIFFERENCE BETWEEN HIGHEST AND SECOND HIGHEST ENERGY LEPTON
	hDiffLepton->Fill(Lepton_p[iL]-Lepton_p[i2L]);
	// ------ MOMENTUM FOR SECOND HIGHEST ENERGY LEPTON
	h2ndHighestEnergyLepton->Fill(Lepton_p[i2L]); 
	//cout  << " Highest = " <<Lepton_p[0] << " Second = " <<Lepton_p[1] <<" Diff = " <<Lepton_p[0]-Lepton_p[1] << endl;
	//--------- Check if highest energy lepton has exactly the same momentum as second highest
	if (Lepton_p[iL]-Lepton_p[i2L]==0) {
	  cout << "Exact lepton momentum match at event= " <<j << endl; 
	  //count1+=1;
	}
	//---------- INVARIANT MASS OF HIGHEST AND SECOND HIGHEST ENERGY LEPTON 
	if (Lepton_pdg[iL]+Lepton_pdg[i2L]==0) { //same type opposite charge
	  vector < array<float,4> > HighestEnergyLeptonSet;
	  HighestEnergyLeptonSet.push_back({Lepton_p[iL], Lepton_px[iL], Lepton_py[iL], Lepton_pz[iL]});
	  HighestEnergyLeptonSet.push_back({Lepton_p[i2L], Lepton_px[i2L], Lepton_py[i2L], Lepton_pz[i2L]});
	  float invMass_highLepton=invMass(HighestEnergyLeptonSet);
	  hInvariantMassHighestLeptons->Fill(invMass_highLepton);
	}

	//---------- LEPTON MOMENTUM EXCLUDING HIGHEST ENERGY LEPTON
	for (UInt_t ilep=1; ilep<Lepton_p.size(); ilep++){
	  int iLepton=sortedLepton[ilep];
	  hLeptonMomentumRest->Fill(Lepton_p[iLepton]);
	}
      }

      for (UInt_t i=0; i<Lepton_p.size();i++){
	hLeptonMomentum->Fill(Lepton_p[i]);	
      }
      
      float invMass_allRP=invMass(allRP);
      hInvariantMassRest->Fill(invMass_rest);
      hInvariantMassAllRP->Fill(invMass_allRP);      
      hInvariantMassLepton->Fill(invMass_lepton);
      hHighestEnergyLepton->Fill(Lepton_p[iL]);      


      //---------- MISSING MOMENTUM
  
      ME_p=365.0-ME_p;
      hMissingMomentum->Fill(ME_p);
      if (MET_p->size()!=1) cout << MET_p->size() << endl;
      hMissingMomentum2->Fill(MET_p->at(0)); 

      //---------- THRUST
      hThrustAllRP->Fill(event_thrust_val);
      hThrustRest->Fill(eventRest_thrust_val);
    
      //---------- B-TAGGING
      hbtagNumber->Fill(Nbtags);
      for (UInt_t i=0; i<jetsignificance->size(); i++) hSignificance->Fill(jetsignificance->at(i));
      
      //---------- dmin
      hdmin->Fill(jets_dmin);
      
      //---------- INVARIANT MASS for each jet
      /*Int_t jetcount=0;
      for (auto constituents : *jetconstituents){
	vector < array<float,4> > constituentset;
	vector < array<float,4> > singlejet;
	singlejet.push_back({jets_e->at(jetcount), jets_px->at(jetcount), jets_py->at(jetcount), jets_pz->at(jetcount)});
	float jetMass=invMass(singlejet);
	hJetMass->Fill(jetMass);
	hJetEnergy->Fill(jets_e->at(jetcount));
	for (unsigned i = 0; i < constituents.size(); i++) {
          //For jet constituent, find corresponding particle in RP collection to get the correct E/px/py/pz
          int idx=RPrest_association->at(constituents[i]); //recojet2RPindex
	  constituentset.push_back({RP_e->at(idx), RP_px->at(idx), RP_py->at(idx), RP_pz->at(idx)});
	}
	float invMass_jet=invMass(constituentset);
	hInvariantMassJet->Fill(invMass_jet);
	hDiffJetMass->Fill(invMass_jet-jetMass);
	jetcount++;
	}*/
      
      //---------- INVARIANT MASS for jet systems
      float minDeltaW = 999;
      float minDeltaTop = 999;
      float minDeltaToplep = 999;
      
      for (UInt_t i=0; i<jets_px->size(); i++){
	//Jet from W decay ~ W/2 mass (does not account for momentum)
	vector < array<float,4> > monojet;
	monojet.push_back({jets_e->at(i), jets_px->at(i), jets_py->at(i), jets_pz->at(i)});
	float monojetMass = invMass(monojet);
	float jetMass=invMass(monojet);
        hJetMass->Fill(jetMass);
	hJetEnergy->Fill(jets_e->at(i));
	float deltamonoJet = abs(monojetMass-MW/2);
	hDeltaWdecay->Fill(deltamonoJet);

	//Jet + lepton + neutrino ~ top mass
	vector < array<float,4> > semijet = leptonset;
	semijet.push_back({jets_e->at(i), jets_px->at(i), jets_py->at(i), jets_pz->at(i)});
	float semijetMass = invMass(semijet);
	float deltasemijet = abs(semijetMass-Mtop);
	if (deltasemijet < minDeltaToplep) minDeltaToplep = deltasemijet;
	
	for (UInt_t k=1; k<jets_px->size(); k++) {
	  //Dijet with ~ W mass
	  vector < array<float,4> > dijet=monojet;
	  dijet.push_back({jets_e->at(k), jets_px->at(k), jets_py->at(k), jets_pz->at(k)});
	  float dijetMass  = invMass(dijet);
	  float deltadijet = abs(dijetMass-MW);
	  if (deltadijet < minDeltaW ) minDeltaW = deltadijet;
	  
	  for (UInt_t l=2; l<jets_px->size(); l++) {
	    //Trijet with ~ top mass
	    vector < array<float,4> > trijet=dijet;
	    trijet.push_back({jets_e->at(k), jets_px->at(k), jets_py->at(k), jets_pz->at(k)});
	    float trijetMass  = invMass(trijet);
	    float deltatrijet = abs(trijetMass-Mtop);
	    if (deltatrijet < minDeltaTop) minDeltaTop = deltatrijet;
	  }
	}
      }
      hDeltaW->Fill(minDeltaW);
      hDeltaTop->Fill(minDeltaTop);
      hDeltaToplep->Fill(minDeltaToplep);
      
      eventcount++;      
    }


    cout << "Number of events with zero leptons " << missinglepton << endl;
    cout << "Number of events with zero jets " << zerojets << endl;
    cout << "Event count " << count << endl;
    cout << ix << " " << sigma[ix] << endl;
    int nEvent = count;//chain->GetEntries(); //changed to count so that the events are scaled to the expected number of events before selection
    float lumi = 1.36*1000000.0; //Integrated luminosity in ab^-1 --> pb^-1     
    float scalefactor = sigma[ix]*lumi/float(nEvent);

    hHighestEnergyLepton->Sumw2();
    h2ndHighestEnergyLepton->Sumw2();
    hLeptonMomentum->Sumw2();
    hLeptonMomentumRest->Sumw2();
    hDiffLepton->Sumw2();
    hMissingMomentum->Sumw2();
    hMissingMomentum2->Sumw2();
    hInvariantMassLepton->Sumw2();
    hInvariantMassHighestLeptons->Sumw2();
    hInvariantMassRest->Sumw2();
    hInvariantMassAllRP->Sumw2();
    hThrustAllRP->Sumw2();
    hThrustRest->Sumw2();
    
    hJetMass->Sumw2();
    hJetEnergy->Sumw2();
    hbtagNumber->Sumw2();
    hSignificance->Sumw2();
    hdmin->Sumw2();
    hDeltaWdecay->Sumw2();
    hDeltaW->Sumw2();
    hDeltaTop->Sumw2();
    hDeltaToplep->Sumw2();
    
    hHighestEnergyLepton->Scale(scalefactor);
    h2ndHighestEnergyLepton->Scale(scalefactor);
    hLeptonMomentum->Scale(scalefactor);
    hLeptonMomentumRest->Scale(scalefactor);
    hDiffLepton->Scale(scalefactor);
    hMissingMomentum->Scale(scalefactor);
    hMissingMomentum2->Scale(scalefactor);
    hInvariantMassLepton->Scale(scalefactor);
    hInvariantMassHighestLeptons->Scale(scalefactor);
    hInvariantMassRest->Scale(scalefactor);
    hInvariantMassAllRP->Scale(scalefactor);
    hThrustAllRP->Scale(scalefactor);
    hThrustRest->Scale(scalefactor);
    
    hJetMass->Scale(scalefactor);
    hJetEnergy->Scale(scalefactor);
    hbtagNumber->Scale(scalefactor);
    hSignificance->Scale(scalefactor);
    hdmin->Scale(scalefactor);
    hDeltaWdecay->Scale(scalefactor);
    hDeltaW->Scale(scalefactor);
    hDeltaTop->Scale(scalefactor);
    hDeltaToplep->Scale(scalefactor);
    
    hHighestEnergyLepton->SetDirectory(outFile);
    h2ndHighestEnergyLepton->SetDirectory(outFile);
    hLeptonMomentum->SetDirectory(outFile);
    hLeptonMomentumRest->SetDirectory(outFile);
    hDiffLepton->SetDirectory(outFile);
    hMissingMomentum->SetDirectory(outFile);
    hMissingMomentum2->SetDirectory(outFile);
    hInvariantMassLepton->SetDirectory(outFile);
    hInvariantMassHighestLeptons->SetDirectory(outFile);
    hInvariantMassRest->SetDirectory(outFile);
    hInvariantMassAllRP->SetDirectory(outFile);
    hThrustAllRP->SetDirectory(outFile);
    hThrustRest->SetDirectory(outFile);
    
    hJetMass->SetDirectory(outFile);
    hJetEnergy->SetDirectory(outFile);
    hbtagNumber->SetDirectory(outFile);
    hSignificance->SetDirectory(outFile);
    hdmin->SetDirectory(outFile);
    hDeltaWdecay->SetDirectory(outFile);
    hDeltaW->SetDirectory(outFile);
    hDeltaTop->SetDirectory(outFile);
    hDeltaToplep->SetDirectory(outFile);
    
    /*
    TCanvas* c1 = new TCanvas();                                                  
    hHighestEnergyLepton->Draw("HIST");
    hHighestEnergyLepton->GetXaxis()->SetTitle("Momentum [GeV]");
    hHighestEnergyLepton->GetYaxis()->SetTitle("Events per 20 GeV"); 
    //hHighestEnergyLepton->SetStats(0); 
    TLegend* legend1 = new TLegend(0.8,0.4,0.95,0.6);
    legend1->SetHeader("","C"); // option "C" allows to center the header
    legend1->AddEntry(hHighestEnergyLepton,mode,"lep");
    legend1->Draw();
    c1->SaveAs("hHighestEnergyLepton"+mode+".pdf");
    delete c1;
    delete legend1;

    TCanvas* c2 = new TCanvas();
    hLeptonMomentum->Draw("HIST");
    hLeptonMomentum->GetXaxis()->SetTitle("Momentum [GeV]");
    hLeptonMomentum->GetYaxis()->SetTitle("Events per 20 GeV");
    //hLeptonMomentum->SetStats(0);
    TLegend* legend2 = new TLegend(0.8,0.4,0.95,0.6);
    legend2->SetHeader("","C"); // option "C" allows to center the header 
    legend2->AddEntry(hLeptonMomentum,mode,"lep");
    legend2->Draw();
    c2->SaveAs("hLeptonMomentum"+mode+".pdf");
    delete c2;
    delete legend2;


    TCanvas* c3 = new TCanvas();
    hDiffLepton->Draw("HIST");
    hDiffLepton->GetXaxis()->SetTitle("#Delta Momentum [GeV]");
    hDiffLepton->GetYaxis()->SetTitle("Events per 20 GeV");
    //hDiffLepton->SetStats(0);
    TLegend* legend3 = new TLegend(0.8,0.4,0.95,0.6);
    legend3->SetHeader("","C"); // option "C" allows to center the header                                                                                         
    legend3->AddEntry(hDiffLepton,mode,"lep");
    legend3->Draw();
    c3->SaveAs("hDiffLepton"+mode+".pdf");
    delete c3;
    delete legend3;

    TCanvas* c4 = new TCanvas();
    hMissingMomentum->Draw("HIST");
    hMissingMomentum->GetXaxis()->SetTitle("Missing Momentum [GeV]");
    hMissingMomentum->GetYaxis()->SetTitle("Events per 20 GeV");
    //hMissingMomentum->SetStats(0);
    TLegend* legend4 = new TLegend(0.8,0.4,0.95,0.6);
    legend4->SetHeader("","C"); // option "C" allows to center the header                                                                                         
    legend4->AddEntry(hMissingMomentum,mode,"lep");
    legend4->Draw();
    c4->SaveAs("hMissingMomentum"+mode+".pdf");
    delete c4;
    delete legend4;

    TCanvas* c5 = new TCanvas();
    hMissingMomentum2->Draw("HIST");
    hMissingMomentum2->GetXaxis()->SetTitle("Missing Momentum [GeV]");
    hMissingMomentum2->GetYaxis()->SetTitle("Events per 20 GeV");
    //hMissingMomentum2->SetStats(0);
    TLegend* legend5 = new TLegend(0.8,0.4,0.95,0.6);
    legend5->SetHeader("","C"); // option "C" allows to center the header                                                                                         
    legend5->AddEntry(hMissingMomentum2,mode,"lep");
    legend5->Draw();
    c5->SaveAs("hMissingMomentum2"+mode+".pdf");
    delete c5;
    delete legend5;

    TCanvas* c6 = new TCanvas();
    hInvariantMassLepton->Draw("HIST");
    hInvariantMassLepton->GetXaxis()->SetTitle("Missing Momentum [GeV]");
    hInvariantMassLepton->GetYaxis()->SetTitle("Events per 20 GeV");
    //hInvariantMassLepton->SetStats(0);
    TLegend* legend6 = new TLegend(0.8,0.5,0.95,0.6);
    legend6->SetHeader("","C"); // option "C" allows to center the header                                                                                         
    legend6->AddEntry(hInvariantMassLepton,mode,"lep");
    legend6->Draw();
    c6->SaveAs("hInvariantMassLepton"+mode+".pdf");
    delete c6;
    delete legend6;

    TCanvas* c7 = new TCanvas();
    hInvariantMassRest->Draw("HIST");
    hInvariantMassRest->GetXaxis()->SetTitle("Missing Momentum [GeV]");
    hInvariantMassRest->GetYaxis()->SetTitle("Events per 20 GeV");
    //hInvariantMassRest->SetStats(0);
    TLegend* legend7 = new TLegend(0.8,0.5,0.95,0.6);
    legend7->SetHeader("","C"); // option "C" allows to center the header                                                                                         
    legend7->AddEntry(hInvariantMassRest,mode,"lep");
    legend7->Draw();
    c7->SaveAs("hInvariantMassRest"+mode+".pdf");
    delete c7;
    delete legend7;

    */
    outFile->Write();
    outFile->Close();

    delete outFile;

    //cout << count1 << endl;
    //break; //only run ttbar_semi

    Passed.push_back(float(eventcount)*scalefactor);
    Sanity.push_back((nEvent-float(missinglepton))*scalefactor);
    //vector<UInt_t> counterror = cutcount;
    //transform(cutcount.begin(), cutcount.end(), cutcount.begin(), [&scalefactor](auto& c){return c*scalefactor;});
    //transform(counterror.begin(), counterror.end(), counterror.begin(), [&scalefactor](auto& c){return scalefactor/sqrt(c);});
    cutflow.push_back(cutcount);
    scalefactors.push_back(scalefactor);
    ix++;
    cout << "Number of events that passed the event selection for " << mode << " = " << eventcount <<endl;

  }

  cout << "----------- Preselection ----------" << endl;
  float Passed_signal=0;
  float Passed_totalbkg=0;
  for (UInt_t i=0;i<Passed.size();i++) {
    if (i<6) Passed_signal+=Passed[i];
    else Passed_totalbkg+=Passed[i];
  }
  for (UInt_t i=5;i<Passed.size();i++) {
    if (i==5){
      cout << "Significance for all backgrounds = " << Passed_signal/sqrt(Passed_signal+Passed_totalbkg) << endl;
      continue;
    }
    cout << "Significance for background " << process[i] << " = " << Passed_signal/sqrt(Passed_signal+Passed[i]) << endl;
  }

  cout << "----------- just sanity check -------- " << endl;
  float Sanity_signal=0;
  float Sanity_totalbkg=0;
  for (UInt_t i=0;i<Sanity.size();i++) {
    if (i<6) Sanity_signal+=Sanity[i];
    else Sanity_totalbkg+=Sanity[i];
  }
  for (UInt_t i=5;i<Sanity.size();i++) {
    if (i==5){
      cout << "Significance for all backgrounds = " << Sanity_signal/sqrt(Sanity_signal+Sanity_totalbkg) << endl;
      continue;
    }
    cout << "Significance for background " << process[i] << " = " << Sanity_signal/sqrt(Sanity_signal+Sanity[i]) << endl;
  }

  
  cout << "----------- Cut flow ----------" << endl;
  UInt_t iMode=0;
  vector<vector<float>> CF;
  vector<vector<float>>	errorCF;
  for (auto& cutcounts : cutflow){
    cout << process[iMode] << endl;
    float scalefactor=scalefactors[iMode];
    vector<float> cf;
    vector<float> errorcf;
    for (UInt_t i=0;i<cutcounts.size();i++){
      cout << "Cut no. " << i << " = " << float(cutcounts[i])*scalefactor << " \\pm " << scalefactor*sqrt(float(cutcounts[i])) << endl;
      cf.push_back(float(cutcounts[i])*scalefactor);
      errorcf.push_back(scalefactor*sqrt(float(cutcounts[i])));
    }
    CF.push_back(cf);
    errorCF.push_back(errorcf);
    iMode++;
  }

  ofstream myfile;
  myfile.open ("Cutflow.txt");
  for (UInt_t i=0; i<CF.size(); i++){
    myfile << process.at(i) << ":";
    for (UInt_t j=0; j<CF[i].size(); j++){
      myfile << CF[i][j] << "," << errorCF[i][j];
      if (j<CF[i].size()-1) myfile << ";" ;
    }
    myfile << endl;
  }
  myfile.close();


  cout << CF[0][0] << " " << CF[1][0]<< endl;
  cout << CF[0][1] << " " << CF[1][1]<< endl;
  
  cout << "----------- Cut flow on signal ----------" << endl;
  cout << "\\toprule " << endl;
  cout << "Selection & $t\\bar{t}, (\\ell=e)$ & $t\\bar{t}, (\\ell=\\mu)$ & $t\\bar{t}, (\\ell=\\tau)$ \\\\ \\midrule" << endl;
  cout << "Initial & " << CF[0][0] + CF[1][0]<< " $\\pm$ " << sqrt(errorCF[0][0]*errorCF[0][0] + errorCF[1][0]*errorCF[1][0]) << " & " << CF[2][0] + CF[3][0]<< " $\\pm$ " << sqrt(errorCF[2][0]*errorCF[2][0] + errorCF[3][0]*errorCF[3][0]) << " & " << CF[4][0] + CF[4][0]<< " $\\pm$ " << sqrt(errorCF[4][0]*errorCF[4][0] + errorCF[5][0]*errorCF[5][0]) << "\\\\" << endl;
  cout << "At least 1 lepton & " << CF[0][1] + CF[1][1]<< " $\\pm$ " << sqrt(errorCF[0][1]*errorCF[0][1] + errorCF[1][1]*errorCF[1][1]) << " & " << CF[2][1] + CF[3][1]<< " $\\pm$ " << sqrt(errorCF[2][1]*errorCF[2][1] + errorCF[3][1]*errorCF[3][1]) << " & " << CF[4][1] + CF[4][1]<< " $\\pm$ " << sqrt(errorCF[4][1]*errorCF[4][1] + errorCF[5][1]*errorCF[5][1]) << "\\\\" << endl;
  cout << "Thrust $<$ 0 .85 & " << CF[0][2] + CF[1][2]<< " $\\pm$ " << sqrt(errorCF[0][2]*errorCF[0][2] + errorCF[1][2]*errorCF[1][2]) << " & " << CF[2][2] + CF[3][2]<< " $\\pm$ " << sqrt(errorCF[2][2]*errorCF[2][2] + errorCF[3][2]*errorCF[3][2]) << " & " << CF[4][2] + CF[4][2]<< " $\\pm$ " << sqrt(errorCF[4][2]*errorCF[4][2] + errorCF[5][2]*errorCF[5][2]) << "\\\\" << endl;
  cout << "$M_(rest)$ $>$ 160 GeV & " << CF[0][3] + CF[1][3]<< " $\\pm$ " << sqrt(errorCF[0][3]*errorCF[0][3] + errorCF[1][3]*errorCF[1][3]) << " & " << CF[2][3] + CF[3][3]<< " $\\pm$ " << sqrt(errorCF[2][3]*errorCF[2][3] + errorCF[3][3]*errorCF[3][3]) << " & " << CF[4][3] + CF[4][3]<< " $\\pm$ " << sqrt(errorCF[4][3]*errorCF[4][3] + errorCF[5][3]*errorCF[5][3]) << "\\\\" << endl;
  cout << "$M_(\\ell_{HE},\\cancel{E})$ $>$ 50 GeV & " << CF[0][4] + CF[1][4]<< " $\\pm$ " << sqrt(errorCF[0][4]*errorCF[0][4] + errorCF[1][4]*errorCF[1][4]) << " & " << CF[2][4] + CF[3][4]<< " $\\pm$ " << sqrt(errorCF[2][4]*errorCF[2][4] + errorCF[3][4]*errorCF[3][4]) << " & " << CF[4][4] + CF[4][4]<< " $\\pm$ " << sqrt(errorCF[4][4]*errorCF[4][4] + errorCF[5][4]*errorCF[5][4]) << "\\\\" << endl;
  cout << "$p_{\\ell_{HE}}$ $<$ 100 GeV & " << CF[0][5] + CF[1][5]<< " $\\pm$ " << sqrt(errorCF[0][5]*errorCF[0][5] + errorCF[1][5]*errorCF[1][5]) << " & " << CF[2][5] + CF[3][5]<< " $\\pm$ " << sqrt(errorCF[2][5]*errorCF[2][5] + errorCF[3][5]*errorCF[3][5]) << " & " << CF[4][5] + CF[4][5]<< " $\\pm$ " << sqrt(errorCF[4][5]*errorCF[4][5] + errorCF[5][5]*errorCF[5][5]) << "\\\\" << endl;
  cout << "$p_{\\ell_{HE}}$ $>$ 15 GeV & " << CF[0][6] + CF[1][6]<< " $\\pm$ " << sqrt(errorCF[0][6]*errorCF[0][6] + errorCF[1][6]*errorCF[1][6]) << " & " << CF[2][6] + CF[3][6]<< " $\\pm$ " << sqrt(errorCF[2][6]*errorCF[2][6] + errorCF[3][6]*errorCF[3][6]) << " & " << CF[4][6] + CF[4][6]<< " $\\pm$ " << sqrt(errorCF[4][6]*errorCF[4][6] + errorCF[5][6]*errorCF[5][6]) << "\\\\" << endl;
  cout << "$p_{\\ell_{2^{nd}HE}}$ $<$ 40 GeV & " << CF[0][7] + CF[1][7]<< " $\\pm$ " << sqrt(errorCF[0][7]*errorCF[0][7] + errorCF[1][7]*errorCF[1][7]) << " & " << CF[2][7] + CF[3][7]<< " $\\pm$ " << sqrt(errorCF[2][7]*errorCF[2][7] + errorCF[3][7]*errorCF[3][7]) << " & " << CF[4][7] + CF[4][7]<< " $\\pm$ " << sqrt(errorCF[4][7]*errorCF[4][7] + errorCF[5][7]*errorCF[5][7]) << "\\\\" << endl;
  cout << "Exactly 4 jets & " << CF[0][8] + CF[1][8]<< " $\\pm$ " << sqrt(errorCF[0][8]*errorCF[0][8] + errorCF[1][8]*errorCF[1][8]) << " & " << CF[2][8] + CF[3][8]<< " $\\pm$ " << sqrt(errorCF[2][8]*errorCF[2][8] + errorCF[3][8]*errorCF[3][8]) << " & " << CF[4][8] + CF[4][8]<< " $\\pm$ " << sqrt(errorCF[4][8]*errorCF[4][8] + errorCF[5][8]*errorCF[5][8]) << "\\\\" << endl;
  cout << "At least 1 b-tag & " << CF[0][9] + CF[1][9]<< " $\\pm$ " << sqrt(errorCF[0][9]*errorCF[0][9] + errorCF[1][9]*errorCF[1][9]) << " & " << CF[2][9] + CF[3][9]<< " $\\pm$ " << sqrt(errorCF[2][9]*errorCF[2][9] + errorCF[3][9]*errorCF[3][9]) << " & " << CF[4][9] + CF[4][9]<< " $\\pm$ " << sqrt(errorCF[4][9]*errorCF[4][9] + errorCF[5][9]*errorCF[5][9]) << "\\\\ \\bottomrule" << endl;


  
	  //Initial & \\ \midrule
	  //At least 1 lepton & \\ 
	  //Thrust $<$ 0 .85 & \\ 
	  //$M_(rest)$ $>$ 160 GeV  & \\ 
	  //$M_(\ell_{HE},\cancel{E})$ $>$ 50 GeV  & \\ 
	  //$p_{\ell_{HE}}$ $<$ 100 GeV & \\ 
	  //$p_{\ell_{HE}}$ $>$ 15 GeV & \\ 
	  //$p_{\ell_{2^{nd}HE}}$ $<$ 40 GeV & \\ 
	  //Exactly 4 jets& \\ 
	  //At least 1 b-tag & \\ \bottomrule
			  
  return 0;
}
