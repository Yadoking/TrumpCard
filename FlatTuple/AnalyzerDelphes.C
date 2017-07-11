#define AnalyzerDelphes_cxx
#include "AnalyzerDelphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TTLJKinFit.h"
#include <TLorentzVector.h>

void AnalyzerDelphes::Loop(const double hadW_m)
{
//   In a ROOT session, you can do:
//      root> .L A.C
//      root> A t
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

  TH1F* hLW_m = new TH1F("hLW_m", "Leptonic W Mass", 100, 0, 300);
  TH1F* hLT_m = new TH1F("hLT_m", "Leptonic Top Mass", 100, 0, 300);
  TH1F* hHW_m = new TH1F("hHW_m", "Hadronic W Mass", 100, 0, 300);
  TH1F* hHW_dR = new TH1F("hHW_dR", "Hadronic W #DeltaR", 100, 0, 5);
  TH1F* hHT_m = new TH1F("hHT_m", "Hadronic Top Mass", 100, 0, 300);
  TH1F* hJES = new TH1F("hJES", "Residual JES", 100, 0, 2);
  TH1F* hHadJJ_dR = new TH1F("hHadJJ_dR", "hHadJJ_dR", 100, 0, 5);
  TH1F* hAddJJ_dR = new TH1F("hAddJJ_dR", "hAddJJ_dR", 100, 0, 5);
  TH1F* hAddJJ_m = new TH1F("hAddJJ_m", "hAddJJ_m", 100, 0, 300);

  if (fChain == 0) return;

  TTLJKinFit fit(hadW_m);

  //Long64_t nentries = 10000;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    cout << jentry << '/' << nentries << '\r';

    TLorentzVector leptonP4;
    if ( muons_n >= 1 ) leptonP4.SetPtEtaPhiM(muons_pt[0], muons_eta[0], muons_phi[0], muons_m[0]);
    else if ( electrons_n >= 1 ) leptonP4.SetPtEtaPhiM(electrons_pt[0], electrons_eta[0], electrons_phi[0], electrons_m[0]);
    else continue;

    if ( leptonP4.Pt() < 30 or std::abs(leptonP4.Eta()) > 2.1 ) continue;

    TLorentzVector metP4;
    metP4.SetPtEtaPhiM(met_pt, 0, met_phi, 0);

    std::vector<size_t> jetIdxs;
    for ( size_t j=0; j<jets_n; ++j ) {
      if ( jets_pt[j] < 30 or std::abs(jets_eta[j]) > 2.5 ) continue;
      jetIdxs.push_back(j);
    }

    double bestChi2 = 1e9;
    std::vector<size_t> bestIdxs;
    TLorentzVector jetP4s[4];
    for ( size_t j1 : jetIdxs ) {
      jetP4s[0].SetPtEtaPhiM(jets_pt[j1], jets_eta[j1], jets_phi[j1], jets_m[j1]);
      for ( size_t j2 : jetIdxs ) {
        if ( j2 == j1 ) continue;
        jetP4s[1].SetPtEtaPhiM(jets_pt[j2], jets_eta[j2], jets_phi[j2], jets_m[j2]);
        for ( size_t j3 : jetIdxs ) {
          if ( j2 >= j3 ) continue;
          if ( j3 == j1 ) continue;
          jetP4s[2].SetPtEtaPhiM(jets_pt[j3], jets_eta[j3], jets_phi[j3], jets_m[j3]);
          for ( size_t j4 : jetIdxs ) {
            if ( j4 == j3 or j4 == j2 or j4 == j1 ) continue;
            jetP4s[3].SetPtEtaPhiM(jets_pt[j4], jets_eta[j4], jets_phi[j4], jets_m[j4]);
            const double chi2 = fit.compute(metP4, leptonP4, jetP4s[0], jetP4s[1], jetP4s[2], jetP4s[3]);

            if ( chi2 < bestChi2 ) {
              bestChi2 = chi2;
              bestIdxs = {j1, j2, j3, j4};
            }
          }
        }
      }
    }
    if ( !bestIdxs.empty() ) {
      for ( size_t i=0; i<4; ++i ) {
        const size_t j = bestIdxs[i];
        jetP4s[i].SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
      }
      const double chi2 = fit.compute(metP4, leptonP4, jetP4s[0], jetP4s[1], jetP4s[2], jetP4s[3]);
      const std::vector<TLorentzVector> solution = fit.getSolution();
      const auto& sol_nuP4 = solution[0], sol_lepP4 = solution[1], sol_ljP4 = solution[2];
      const auto& sol_wj1P4 = solution[3], sol_wj2P4 = solution[4], sol_hbP4 = solution[5];

      const double hadJJ_dR = sol_wj1P4.DeltaR(sol_wj2P4);
      hHadJJ_dR->Fill(hadJJ_dR);
      hLW_m->Fill( (sol_lepP4+sol_nuP4).M() );
      hLT_m->Fill( (sol_lepP4+sol_nuP4+sol_ljP4).M() );
      hHW_m->Fill( (sol_wj1P4+sol_wj2P4).M() );
      hHW_dR->Fill( sol_wj1P4.DeltaR(sol_wj2P4) );
      hHT_m->Fill( (sol_wj1P4+sol_wj2P4+sol_hbP4).M() );

      hJES->Fill(fit.min_->X()[0]);

      std::vector<size_t> addJetIdxs;
      for ( size_t j : jetIdxs ) {
        bool isTaken = false;
        for ( size_t bestIdx : bestIdxs ) {
          if ( bestIdx == j ) { isTaken = true; break; }
        }
        if ( isTaken ) continue;

        addJetIdxs.push_back(j);
      }
      if ( addJetIdxs.size() >= 2 ) {
        std::sort(addJetIdxs.begin(), addJetIdxs.end(), [&](size_t a, size_t b) { return jets_pt[a] > jets_pt[b]; });
        std::stable_sort(addJetIdxs.begin(), addJetIdxs.end(), [&](size_t a, size_t b) { return jets_bTag[a]==1; });
        const size_t j1 = addJetIdxs[0], j2 = addJetIdxs[1];
        TLorentzVector addJet1, addJet2;
        addJet1.SetPtEtaPhiM(jets_pt[j1], jets_eta[j1], jets_phi[j1], jets_m[j1]);
        addJet1.SetPtEtaPhiM(jets_pt[j2], jets_eta[j2], jets_phi[j2], jets_m[j2]);
        const double dR = addJet1.DeltaR(addJet2);
        const double m = (addJet1+addJet2).M();
        hAddJJ_dR->Fill(dR);
        hAddJJ_m->Fill(m);
      }
    }
  }

  TFile* fout = new TFile("hist.root", "recreate");
  hLW_m->SetDirectory(fout);
  hLT_m->SetDirectory(fout);
  hHW_m->SetDirectory(fout);
  hHW_dR->SetDirectory(fout);
  hHT_m->SetDirectory(fout);
  hJES->SetDirectory(fout);
  fout->Write();

  TCanvas* c = 0;
  c = new TCanvas("cHadJJ_dR", "hadJJ dR", 500, 500); hHadJJ_dR->Draw();
  c = new TCanvas("cAddJJ_dR", "addJJ dR", 500, 500); hAddJJ_dR->Draw();
  c = new TCanvas("cHadJJ_m", "hadJJ m", 500, 500); hAddJJ_m->Draw();
  c = new TCanvas("cLW_m", "MW lep", 500, 500); hLW_m->Draw();
  c = new TCanvas("cLT_m", "MTop lep", 500, 500); hLT_m->Draw();
  c = new TCanvas("cHW_m", "MW had", 500, 500); hHW_m->Draw();
  c = new TCanvas("cHW_dR", "dR W had", 500, 500); hHW_dR->Draw();
  c = new TCanvas("cHT_m", "MTop had", 500, 500); hHT_m->Draw();
  c = new TCanvas("cJES", "JES", 500, 500); hJES->Draw();
}
