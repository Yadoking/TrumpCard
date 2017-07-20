#define Yosi_cxx
#include "Yosi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Yosi::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Yosi.C
//      root> Yosi t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TFile* fout = new TFile("temp.root", "recreate");
   TTree* tree = new TTree("namu", "namu");

   TH1F* hdR = new TH1F("dR", "dR", 100, 0, 5);
   TH1F* hdijetM = new TH1F("dijet_Mass", "dijet Mass", 300, 0, 300);
   TH1F* hNJet = new TH1F("NJet", "NJet", 10, 0, 10);
   TH1F* hBJetPt_H = new TH1F("BJetPt", "BJetPt", 150, 0, 300); 
   TH1F* hm3 = new TH1F("M3", "M3", 1000, 0, 2000);


   int b_njet;
   int b_nbjet;

   float b_lepton_pt;
   float b_lepton_eta;
   float b_lepton_phi;
   float b_lepton_m;

   float b_jet_pt;
   float b_jet_eta;
   float b_jet_phi;
   float b_jet_m;

   float b_deltaR;
   float b_diJetMass;
   float b_maxJetPt;
   
   float b_m3;

   tree->Branch("njet", &b_njet, "njet/I");
   tree->Branch("nbjet", &b_nbjet, "nbjet/I");
   
   tree->Branch("lepton_pt", &b_lepton_pt, "lepton_pt/F");
   tree->Branch("lepton_eta", &b_lepton_eta, "lepton_eta/F");
   tree->Branch("lepton_phi", &b_lepton_phi, "lepton_phi/F");
   tree->Branch("lepton_m", &b_lepton_m, "lepton_m/F");

   tree->Branch("jet_pt", &b_jet_pt, "jet_pt/F");
   tree->Branch("jet_eta", &b_jet_eta, "jet_eta/F");
   tree->Branch("jet_phi", &b_jet_phi, "jet_phi/F");
   tree->Branch("jet_m", &b_jet_m, "jet_m/F");

   tree->Branch("deltaR", &b_deltaR, "deltaR/F");
   tree->Branch("diJetMass", &b_diJetMass, "diJetMass/F");
   tree->Branch("maxJetPt", &b_maxJetPt, "maxJetPt");

   tree->Branch("M3", &b_m3, "M3/F");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      //
      // Yosi!Programming Season!
      
      int njets = 0;
      int nbjets = 0;
      TLorentzVector leptonP4;

      if ( muons_n >= 1 and muons_pt[0] > 30 and std::abs(muons_eta[0]) < 2.5)   leptonP4.SetPtEtaPhiM(muons_pt[0], muons_eta[0], muons_phi[0], muons_m[0]);
      else if ( electrons_n >= 1 and electrons_pt[0] and std::abs(electrons_eta[0]) < 2.5)   leptonP4.SetPtEtaPhiM(electrons_pt[0], electrons_eta[0], electrons_phi[0], electrons_m[0]);
      else continue;
      
      b_lepton_pt = leptonP4.Pt();
      b_lepton_eta = leptonP4.Eta();
      b_lepton_phi = leptonP4.Phi();
      b_lepton_m = leptonP4.M();

      vector<TLorentzVector> p4s;

      for (int i = 0; i < jets_n; ++i)
      {
         if (jets_pt[i] < 30 or std::abs(jets_eta[i]) > 2.5) continue;
         njets++;
         if (jets_bTag[i] == 1)
         {
            ++nbjets;
            TLorentzVector p4;
            p4.SetPtEtaPhiM(jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
            p4s.push_back(p4);
            b_jet_pt = p4.Pt();
            b_jet_eta = p4.Eta();
            b_jet_phi = p4.Phi();
            b_jet_m = p4.M();
         }
      }
      b_njet = njets;
      b_nbjet = nbjets;
      // m3 variable
      vector<TLorentzVector> keepm3;
      if (njets > 2)
      {
         double m3pt = 0;
         for (int i = 2; i < njets; ++i)
         {
            for (int j = 1; j < i; ++j)
            {
               for (int k = 0; k < j; ++k)
               {
                  TLorentzVector j1, j2, j3;
                  j1.SetPtEtaPhiM(jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
                  j2.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
                  j3.SetPtEtaPhiM(jets_pt[k], jets_eta[k], jets_phi[k], jets_m[k]);
                  const double M3pt = (j1+j2+j3).Pt();
		  if (M3pt >= m3pt){
                     m3pt = M3pt;
                     keepm3 = { j1, j2, j3 };
                  }
               }
            }
         }
         b_m3 = (keepm3[0]+keepm3[1]+keepm3[2]).M();
      }
      hm3->Fill(b_m3);
      // minimum delta R
      if (nbjets > 1)
      {
         double minDR = 999, minM = -1;
         double maxPt = 0;
         for (int i = 1; i < nbjets; ++i)
         {
            const auto& Bjet1 = p4s.at(i);
            for (int j = 0; j < i; ++j)
            {
               const auto& Bjet2 = p4s.at(j);
               const double dR = Bjet1.DeltaR(Bjet2);
               const double PT = Bjet1.Pt() >= Bjet2.Pt() ? Bjet1.Pt() : Bjet2.Pt();
               if ( dR < minDR ) {
                 // using Lambda
                 /*
                 vector<pair<double, double> > mins;
                 std::sort(mins.begin(), mins.end(),
                 [&](  pair<double, double>a, pair<double, double> b)
                 {
                    return a.first < b.first;
                 }
                 );
                 */

                 minDR = dR;
                 minM = (Bjet1+Bjet2).M();
                 maxPt = PT;
               }
            }
         }
         b_deltaR = minDR;
         b_diJetMass = minM;
         b_maxJetPt = maxPt;

	 hdR->Fill(minDR);
	 hdijetM->Fill(minM);
         hNJet->Fill(njets);
         hBJetPt_H->Fill(maxPt);
      }
      tree->Fill();
   }
   fout->Write();
   /*
   TCanvas* c = new TCanvas("cDr", "delta R", 500, 500);
   hdR->Draw();
   c = new TCanvas("cdiM", "dijet Mass", 500, 500);
   hdijetM->Draw();
   c = new TCanvas("cNJ", "N Jet", 500, 500);
   hNJet->Draw();
   c = new TCanvas("cPt", "Max Pt", 500, 500);
   hBJetPt_H->Draw();
   */
}
