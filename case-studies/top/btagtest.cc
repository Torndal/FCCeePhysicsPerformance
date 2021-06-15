//To run: g++ -o btagtest -std=gnu++17 -I /Users/Julie/ROOT/include btagtest.cc `root-config --cflags --libs`

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

#include "TAxis.h"
#include "TGaxis.h"
#include <cmath>
#include <cassert>

using namespace std;

TFile *outFile;


template <typename T>
inline std::vector<T> range(T a, T b, T dist) {
  assert(b >= a);
  int n = (int)std::ceil((b-a)/dist);
  std::vector<T> ret(n);
  T newElement = a + dist;
  ret[0] = a;
  for(unsigned int i = 1; i < n; i++) {
    ret[i] = newElement;
    newElement += dist;
  }
  return ret;
}

int main()
{
  Float_t tol = 1.0e-16;

  TChain* chain = new TChain("events");
  chain->Add("outputs/semileptonic/btag/p8_ee_ttbar_semi_ecm365_IDEAtrkCov.root");

  outFile = new TFile("btagtest.root", "RECREATE");

  vector<float> *jets_btag=0;
  vector<float> *jets_truebtag=0;
  vector<vector<int>> *jetconstituents=0;
  vector<float> *jetsignificance=0;
  vector<float> *jetdistance=0;
  vector<float> *jetsigma=0;
  vector<float> *jetchi2=0;
  
  
  vector<float> *MC_pdg=0;
  vector<float> *MC_status=0;
  vector<int> *RPrest_association=0;
  vector<int> *RPMC_index=0;
  vector<int> *MC_parent1=0;
  vector<int> *MC_parent2=0;
  
  chain->SetBranchAddress("jets_btag", & jets_btag);
  chain->SetBranchAddress("jets_truebtag", & jets_truebtag);
  chain->SetBranchAddress("jetconstituents", & jetconstituents);
  chain->SetBranchAddress("jetsignificance", & jetsignificance);
  chain->SetBranchAddress("jetdistance", & jetdistance);
  chain->SetBranchAddress("jetsigma", & jetsigma);
  chain->SetBranchAddress("jetchi2", & jetchi2);
  
  chain->SetBranchAddress("MC_pdg", & MC_pdg);
  chain->SetBranchAddress("MC_status", & MC_status);
  chain->SetBranchAddress("RPrest_association", & RPrest_association);
  chain->SetBranchAddress("RPMC_index", & RPMC_index);
  chain->SetBranchAddress("MC_parent1", & MC_parent1);
  chain->SetBranchAddress("MC_parent2", & MC_parent2);
  
  TH1F *hbtagNumber = new TH1F("hbtagNumber", "Number of b-tags per event", 5,0,5);
  TH1F *htruebtagNumber = new TH1F("htruebtagNumber", "Number of b-tags per event", 5,0,5);
  TH1F *hSignificance = new TH1F("hSignificance", "Significance distribution", 50,0,5000);
  TH1F *hbtagSignificance = new TH1F("hbtagSignificance", "Signifance distibution for b-tagged jets", 50,0,5000);
  TH1F *hnotagSignificance = new TH1F("hnotagSignificance", "Signifance distibution for not b-tagged jets", 50,0,5000);
  TH1F *hAllDecaysSignificance = new TH1F("hAllDecaysSignificance", "Signifance distibution for jets with all b-hadron decay products", 50,0,5000);
  TH1F *hNoDecaysSignificance = new TH1F("hNoDecaysSignificance", "Signifance distibution for jets with no b-hadron decay products", 50,0,5000);

  TH1F *hNumberOfTracks = new TH1F("hNumberOfTracks", "Number of Constituents in a b-tagged jet", 50,0,50);

  TH1F *hFlightDistance = new TH1F("hFlightDistance", "Flight Distance", 100,0,1000);
  TH1F *hFlightDistErr  = new TH1F("hFlightDistErr",  "Error on flight distance", 100,0,1);
  TH1F *hChi2 = new TH1F("hChi2", "Chi2/ndf",40,0,20);

  TH1F *hbtagFlightDistance = new TH1F("hbtagFlightDistance", "Flight Distance for b-tagged jets", 100,0,1000);
  TH1F *hbtagFlightDistErr  = new TH1F("hbtagFlightDistErr",  "Error on flight distance for b-tagged jets", 100,0,1);
  TH1F *hbtagChi2 =	new TH1F("hbtagChi2", "Chi2/ndf for b-tagged jets",40,0,20);

  //TH1F *hnotagFlightDistance = new TH1F("hnotagFlightDistance", "Flight Distance for not b-tagged jets", 100,0,1000);
  //TH1F *hnotagFlightDistErr  = new TH1F("hnotagFlightDistErr",  "Error on flight distance for not b-tagged jets", 100,0,1);
  //TH1F *hnotagChi2 =	new TH1F("hnotagChi2", "Chi2/ndf for not b-tagged jets",40,0,20);

  TH1F *hAllDecaysFlightDistance = new TH1F("hAllDecaysFlightDistance", "Flight Distance for jets with all b-hadron decay products", 100,0,1000);
  TH1F *hAllDecaysFlightDistErr  = new TH1F("hAllDecaysFlightDistErr",  "Error on flight distance for jets with all b-hadron decay products", 100,0,1);
  TH1F *hAllDecaysChi2 =	new TH1F("hAllDecaysChi2", "Chi2/ndf for jets with all b-hadron decay products",40,0,20);

  TH1F *hNoDecaysFlightDistance = new TH1F("hNoDecaysFlightDistance", "Flight Distance for jets with no b-hadron decay products", 100,0,1000);
  TH1F *hNoDecaysFlightDistErr  = new TH1F("hNoDecaysFlightDistErr",  "Error on flight distance for jets with no b-hadron decay products", 100,0,1);
  TH1F *hNoDecaysChi2 =	new TH1F("hNoDecaysChi2", "Chi2/ndf for jets with no b-hadron decay products",40,0,20);

  
  for (UInt_t j=0; j<chain->GetEntries();j++){
    chain->GetEntry(j);
    int count=0;
    int truecount=0;
    for (UInt_t i=0; i<jets_btag->size();i++){
      if (jets_btag->at(i)==1) count++;
      if (jets_truebtag->at(i)==1) truecount++;
    }
    hbtagNumber->Fill(count);
    htruebtagNumber->Fill(truecount);
    
    //unordered_multiset<int> btaggedRecoJets;
    vector<unordered_multiset<float>> allrecobtags;
    for (auto constituents : *jetconstituents){
      hNumberOfTracks->Fill(constituents.size());
      ///**unordered_set<float> bhadrons;
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
	      ///**unordered_set<float>::const_iterator got = bhadrons.find (MC_pdg->at(index));
	      ///**if ( got == bhadrons.end() ) bhadrons.insert(MC_pdg->at(index));
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
      ///**btaggedRecoJets.push_back(bhadrons.size());
      ///**for (auto bhadron : bhadrons) allrecobhadrons.insert(bhadron);
      allrecobtags.push_back(btags);
    }

    unordered_multiset<float> totbtags;
    for (auto btags : allrecobtags) totbtags.insert(btags.begin(), btags.end());
    unordered_set<float> totbhadrons(totbtags.begin(), totbtags.end());

    for (UInt_t i=0; i<jetsignificance->size(); i++) {
      hSignificance->Fill(jetsignificance->at(i));
      hFlightDistance->Fill(jetdistance->at(i));
      hFlightDistErr->Fill(jetsigma->at(i));
      hChi2->Fill(jetchi2->at(i));
      
      if (jets_truebtag->at(i)==1){
	hbtagSignificance->Fill(jetsignificance->at(i));
	hbtagFlightDistance->Fill(jetdistance->at(i));
	hbtagFlightDistErr->Fill(jetsigma->at(i));
	hbtagChi2->Fill(jetchi2->at(i));
      }
      auto btags =allrecobtags[i];
      unordered_set<float> bhadrons(btags.begin(), btags.end());
      if (bhadrons.size() == 1) {
	for (auto bhadron: bhadrons) {
	  if (btags.count(bhadron) == totbtags.count(bhadron)) {
	    hAllDecaysSignificance->Fill(jetsignificance->at(i));
	    hAllDecaysFlightDistance->Fill(jetdistance->at(i));
	    hAllDecaysFlightDistErr->Fill(jetsigma->at(i));
	    hAllDecaysChi2->Fill(jetchi2->at(i));
	  }
	}
      }
      //if (bhadrons.size()==1 && btags.count(bhadrons.begin())==totbtags.count(bhadrons.begin())) hAllDecaysSignificance->Fill(jetsignificance->at(i));
      if (bhadrons.size()==0) {
	hNoDecaysSignificance->Fill(jetsignificance->at(i));
	hNoDecaysFlightDistance->Fill(jetdistance->at(i));
	hNoDecaysFlightDistErr->Fill(jetsigma->at(i));
	hNoDecaysChi2->Fill(jetchi2->at(i));
      }
    }
    
    
  }
  
  hbtagNumber->SetDirectory(outFile);
  htruebtagNumber->SetDirectory(outFile);

  hSignificance->SetDirectory(outFile);
  hbtagSignificance->SetDirectory(outFile);
  hnotagSignificance->SetDirectory(outFile);
  hAllDecaysSignificance->SetDirectory(outFile);
  hNoDecaysSignificance->SetDirectory(outFile);

  hFlightDistance->SetDirectory(outFile);
  hbtagFlightDistance->SetDirectory(outFile);
  hAllDecaysFlightDistance->SetDirectory(outFile);
  hNoDecaysFlightDistance->SetDirectory(outFile);
  
  hFlightDistErr->SetDirectory(outFile);
  hbtagFlightDistErr->SetDirectory(outFile);
  hAllDecaysFlightDistErr->SetDirectory(outFile);
  hNoDecaysFlightDistErr->SetDirectory(outFile);
  
  hChi2->SetDirectory(outFile);
  hbtagChi2->SetDirectory(outFile);
  hAllDecaysChi2->SetDirectory(outFile);
  hNoDecaysChi2->SetDirectory(outFile);
  
  
  
  TCanvas* c1 = new TCanvas("cst", "b-tags", 10,10,700,500);
  htruebtagNumber->GetXaxis()->SetTitle("#b-tags");
  htruebtagNumber->GetYaxis()->SetTitle("Event rate per 1");
  htruebtagNumber->GetXaxis()->SetTitleSize(0.05);
  htruebtagNumber->GetYaxis()->SetTitleSize(0.05);
  htruebtagNumber->GetXaxis()->SetTitleOffset(0.9);
  htruebtagNumber->GetYaxis()->SetTitleOffset(0.9);
  htruebtagNumber->SetMarkerColor(kRed);
  htruebtagNumber->SetMarkerStyle(1);
  htruebtagNumber->SetLineColor(kRed);
  hbtagNumber->SetMarkerColor(kBlack);
  hbtagNumber->SetMarkerStyle(1);
  hbtagNumber->SetLineColor(kBlack);
  htruebtagNumber->Draw();
  hbtagNumber->Draw("same");
  
  TLegend *legend = new TLegend(0.15,0.65,0.35,0.85);
  legend->SetHeader("","C"); // option "C" allows to center the header                                                                                                           
  legend->AddEntry(htruebtagNumber,"100 % efficiency","l");
  legend->AddEntry(hbtagNumber,"80 % efficiency","l");
  legend->SetTextSize(0.03);
  legend->Draw();
  
  c1->SaveAs("bTagPlots/hbtags.png");
  delete c1;

  TCanvas* cTrack = new TCanvas("cst", "b-tags", 10,10,700,500);
  hNumberOfTracks->GetXaxis()->SetTitle("#b-tags");
  hNumberOfTracks->GetYaxis()->SetTitle("Event rate per 1");
  hNumberOfTracks->GetXaxis()->SetTitleSize(0.05);
  hNumberOfTracks->GetYaxis()->SetTitleSize(0.05);
  hNumberOfTracks->GetXaxis()->SetTitleOffset(0.9);
  hNumberOfTracks->GetYaxis()->SetTitleOffset(0.9);
  hNumberOfTracks->SetMarkerColor(kRed);
  hNumberOfTracks->SetMarkerStyle(1);
  hNumberOfTracks->SetLineColor(kRed);
  hNumberOfTracks->Draw();

  TLegend *legendTrack = new TLegend(0.15,0.65,0.35,0.85);
  legendTrack->SetHeader("","C"); // option "C" allows to center the header                                                                                                                                                                                                          
  legendTrack->AddEntry(hNumberOfTracks,"Number of tracks associated with a jet","l");
  legendTrack->SetTextSize(0.03);
  legendTrack->Draw();

  cTrack->SaveAs("bTagPlots/hNumberOfTracks.png");
  delete cTrack;
  
  TCanvas* c2 = new TCanvas("cst","stacked hists",10,10,700,500);
  c2->SetLogy();
  //c2->SetLogx();
  hSignificance->GetXaxis()->SetTitle("Significance");
  hSignificance->GetYaxis()->SetTitle("Event rate per 100");
  hSignificance->GetXaxis()->SetTitleSize(0.05);
  hSignificance->GetYaxis()->SetTitleSize(0.05);
  hSignificance->GetXaxis()->SetTitleOffset(0.9);
  hSignificance->GetYaxis()->SetTitleOffset(0.9);
  hSignificance->SetLineColor(kBlack);
  hbtagSignificance->SetLineColor(kBlue);
  hAllDecaysSignificance->SetLineColor(kRed);
  hNoDecaysSignificance->SetLineColor(kGray+1);
  hSignificance->Draw();
  hbtagSignificance->Draw("same");
  hAllDecaysSignificance->Draw("same");
  hNoDecaysSignificance->Draw("same");
  
  TLegend *legend2 = new TLegend(0.15,0.65,0.75,0.85);
  legend2->SetHeader("","C"); // option "C" allows to center the header                                                                                                                                                        
  legend2->AddEntry(hSignificance,"Significance for all jets","l");
  legend2->AddEntry(hbtagSignificance,"Signifince for jets with b-tag","l");
  legend2->AddEntry(hAllDecaysSignificance,"Signifince for jets with all b-hadron decay products","l");
  legend2->AddEntry(hNoDecaysSignificance,"Signifince for jets with no b-hadron decay products","l");
  legend2->SetTextSize(0.03);
  legend2->Draw();
  
  c2->SaveAs("bTagPlots/hSignificance.png");
  delete c2;


  TCanvas* cDist = new TCanvas("cst","Flight Distance",10,10,700,500);
  cDist->SetLogy();
  //cDist->SetLogx();
  hFlightDistance->GetXaxis()->SetTitle("Flight Distance");
  hFlightDistance->GetYaxis()->SetTitle("Event rate per 100");
  hFlightDistance->GetXaxis()->SetTitleSize(0.05);
  hFlightDistance->GetYaxis()->SetTitleSize(0.05);
  hFlightDistance->GetXaxis()->SetTitleOffset(0.9);
  hFlightDistance->GetYaxis()->SetTitleOffset(0.9);
  hFlightDistance->SetLineColor(kBlack);
  hbtagFlightDistance->SetLineColor(kBlue);
  hAllDecaysFlightDistance->SetLineColor(kRed);
  hNoDecaysFlightDistance->SetLineColor(kGray+1);
  hFlightDistance->Draw();
  hbtagFlightDistance->Draw("same");
  hAllDecaysFlightDistance->Draw("same");
  hNoDecaysFlightDistance->Draw("same");

  TLegend *legendDist = new TLegend(0.15,0.65,0.75,0.85);
  legendDist->SetHeader("","C"); // option "C" allows to center the header
  legendDist->AddEntry(hFlightDistance,"Flight distance for all jets","l");
  legendDist->AddEntry(hbtagFlightDistance,"Flight distance for jets with b-tag","l");
  legendDist->AddEntry(hAllDecaysFlightDistance,"Flight distance for jets with all b-hadron decay products","l");
  legendDist->AddEntry(hNoDecaysFlightDistance,"Flight distance for jets with no b-hadron decay products","l");
  legendDist->SetTextSize(0.03);
  legendDist->Draw();

  cDist->SaveAs("bTagPlots/hFlightDistance.png");
  delete cDist;


  TCanvas* cSigma = new TCanvas("cst","Error on Flight Distance",10,10,700,500);
  cSigma->SetLogy();
  //cSigma->SetLogx();
  hFlightDistErr->GetXaxis()->SetTitle("Error on Flight Distance");
  hFlightDistErr->GetYaxis()->SetTitle("Event rate per 100");
  hFlightDistErr->GetXaxis()->SetTitleSize(0.05);
  hFlightDistErr->GetYaxis()->SetTitleSize(0.05);
  hFlightDistErr->GetXaxis()->SetTitleOffset(0.9);
  hFlightDistErr->GetYaxis()->SetTitleOffset(0.9);
  hFlightDistErr->SetLineColor(kBlack);
  hbtagFlightDistErr->SetLineColor(kBlue);
  hAllDecaysFlightDistErr->SetLineColor(kRed);
  hNoDecaysFlightDistErr->SetLineColor(kGray+1);
  hFlightDistErr->Draw();
  hbtagFlightDistErr->Draw("same");
  hAllDecaysFlightDistErr->Draw("same");
  hNoDecaysFlightDistErr->Draw("same");

  TLegend *legendSigma = new TLegend(0.15,0.65,0.75,0.85);
  legendSigma->SetHeader("","C"); // option "C" allows to center the header
  legendSigma->AddEntry(hFlightDistErr,"Error for all jets","l");
  legendSigma->AddEntry(hbtagFlightDistErr,"Error for jets with b-tag","l");
  legendSigma->AddEntry(hAllDecaysFlightDistErr,"Error for jets with all b-hadron decay products","l");
  legendSigma->AddEntry(hNoDecaysFlightDistErr,"Error for jets with no b-hadron decay products","l");
  legendSigma->SetTextSize(0.03);
  legendSigma->Draw();

  cSigma->SaveAs("bTagPlots/hFlightDistErr.png");
  delete cSigma;

  TCanvas* cChi2 = new TCanvas("cst","Chi2/ndf",10,10,700,500);
  cChi2->SetLogy();
  //cChi2->SetLogx();
  hChi2->GetXaxis()->SetTitle("Chi2/ndf");
  hChi2->GetYaxis()->SetTitle("Event rate per 100");
  hChi2->GetXaxis()->SetTitleSize(0.05);
  hChi2->GetYaxis()->SetTitleSize(0.05);
  hChi2->GetXaxis()->SetTitleOffset(0.9);
  hChi2->GetYaxis()->SetTitleOffset(0.9);
  hChi2->SetLineColor(kBlack);
  hbtagChi2->SetLineColor(kBlue);
  hAllDecaysChi2->SetLineColor(kRed);
  hNoDecaysChi2->SetLineColor(kGray+1);
  hChi2->Draw();
  hbtagChi2->Draw("same");
  hAllDecaysChi2->Draw("same");
  hNoDecaysChi2->Draw("same");

  TLegend *legendChi2 = new TLegend(0.15,0.65,0.75,0.85);
  legendChi2->SetHeader("","C"); // option "C" allows to center the header
  legendChi2->AddEntry(hChi2,"Error for all jets","l");
  legendChi2->AddEntry(hbtagChi2,"Error for jets with b-tag","l");
  legendChi2->AddEntry(hAllDecaysChi2,"Error for jets with all b-hadron decay products","l");
  legendChi2->AddEntry(hNoDecaysChi2,"Error for jets with no b-hadron decay products","l");
  legendChi2->SetTextSize(0.03);
  legendChi2->Draw();

  cChi2->SaveAs("bTagPlots/hChi2.png");
  delete cChi2;

  


  
  Float_t tot_sig=hbtagSignificance->Integral(); //Get integral of entire signal distribution for efficiency
  TString GraphTitle=(TString)hbtagSignificance->GetTitle();
  TAxis *xaxis = hbtagSignificance->GetXaxis();

  vector<Int_t> bins;
  vector<float> cuts = range((float)xaxis->GetXmin(), (float)xaxis->GetXmax(), (float)xaxis->GetXmax()/100);

  for(auto const & cut : cuts) {
    bins.push_back(xaxis->FindBin(cut));
  }
    int n=cuts.size();
    vector<Float_t> sig_notag; //integral for significance without b-tag
    vector<Float_t> sig_btag; //integral for significance with b-tag
    
    for (int i = 0; i < n; ++i) {
      Int_t bincut=bins[i];
      sig_notag.push_back(hnotagSignificance->Integral(bincut,-1));  //Integral from lower cutvalue (-1 gives the last bin)
      sig_btag.push_back(hbtagSignificance->Integral(bincut,-1));  //Integral from lower cutvalue (-1 gives the last bin)
    }
    

    vector<Float_t> cuts_low;
    vector<Float_t> Eff_low;
    vector<Float_t> Puri_low;
    vector<Float_t> Significance_low;
    vector<Float_t> PE_low;

    vector<Float_t> cuts_up;
    vector<Float_t> Eff_up;
    vector<Float_t> Puri_up;
    vector<Float_t> Significance_up;
    vector<Float_t> PE_up;

    Float_t sigMax_low=0;
    for (int i = 0; i < n; ++i) {
      if (sig_btag[i]+sig_notag[i]>tol) {
        cuts_low.push_back(cuts[i]);
        Float_t t_Eff_low = sig_btag[i]/tot_sig;
        Float_t t_Puri_low = sig_btag[i]/(sig_btag[i]+sig_notag[i]);
        Float_t t_Significance_low =sig_btag[i]/sqrt(sig_btag[i]+sig_notag[i]);
        Eff_low.push_back(t_Eff_low);
        Puri_low.push_back(t_Puri_low);
        Significance_low.push_back(t_Significance_low);
        PE_low.push_back(t_Eff_low*t_Puri_low);

        if (sigMax_low<t_Significance_low) sigMax_low=t_Significance_low;
      }
    }

    TCanvas *c3 = new TCanvas("c3","Lower cut",200,10,500,300);
    TGraph* gr1 = new TGraph();
    for (UInt_t i=0; i<PE_low.size();i++) gr1->SetPoint(i,cuts_low[i],PE_low[i]);
    gr1->GetXaxis()->SetTitle("cut");
    gr1->GetYaxis()->SetTitle("Efficiency");
    gr1->SetTitle("Lower cut on "+GraphTitle);
    gr1->SetMarkerColor(kBlack);
    gr1->SetMarkerSize(0.7);
    gr1->SetMarkerStyle(8);
    TGraph* gr1p = new TGraph();
    for (UInt_t i=0; i<PE_low.size();i++) gr1p->SetPoint(i,cuts_low[i],Puri_low[i]);
    gr1p->SetLineColor(kBlue);
    gr1p->SetLineWidth(5);
    TGraph* gr1e = new TGraph();
    for (UInt_t i=0; i<PE_low.size();i++) gr1e->SetPoint(i,cuts_low[i],Eff_low[i]);
    gr1e->SetLineColor(kRed);
    gr1e->SetLineWidth(5);
    TGraph* gr1s = new TGraph();
    for (UInt_t i=0; i<PE_low.size();i++) gr1s->SetPoint(i,cuts_low[i],Significance_low[i]);
    gr1s->SetLineColor(kMagenta);
    gr1s->SetLineWidth(5);


    auto legend3 = new TLegend(0.7,0.4,0.9,0.6);
    legend3->SetHeader("","C"); // option "C" allows to center the header                                                                                                                                                        
    legend3->AddEntry(gr1,"PE","lep");
    legend3->AddEntry(gr1p,"Purity","lep");
    legend3->AddEntry(gr1e,"Efficiency","lep");
    legend3->AddEntry(gr1s,"Significance","lep");
    legend3->SetTextSize(0.04);

    gr1->GetHistogram()->SetMaximum(1.1);
    gr1->GetXaxis()->SetTitleSize(0.05);
    gr1->GetYaxis()->SetTitleSize(0.05);
    gr1->GetXaxis()->SetTitleOffset(0.9);
    gr1->GetYaxis()->SetTitleOffset(0.9);

    gr1->Draw("AP");
    gr1e->Draw("C SAME");
    gr1p->Draw("C SAME");
    legend3->Draw();

// scale hint1 to the pad coordinates                                                                                                                                                                                       
    float rightmax = 1.1*TMath::MaxElement(gr1s->GetN(),gr1s->GetY());
    float scale = 1.1*gPad->GetUymax()/rightmax;
    //cout << "ladida" << rightmax << " " << scale << endl;                                                                                                                                                                     
    for (int i=0;i<gr1s->GetN();i++) {
      gr1s->GetY()[i] *= scale; //gr1s->Scale(scale) (TH1) equivalent for TGraph                                                                                                                                                
      //cout << cuts_low[i] << " " << gr1s->GetY()[i] << endl;                                                                                                                                                                  
    }
    gr1s->Draw("C SAME");
    c3->Update();
    // draw an axis on the right side                                                                                                                                                                                           
    auto axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,50510,"+L");
    axis->SetLineColor(kMagenta);
    axis->SetTextColor(kMagenta);
    axis->Draw();

    c3->SaveAs("bTagPlots/SigofSig.pdf");
    delete c3;

    cout << " Maximum significance for lower cut" << sigMax_low << endl;
  
    
  outFile->Write();
  outFile->Close();
  delete outFile;
  return 0;
}
