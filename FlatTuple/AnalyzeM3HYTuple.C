#define AnalyzeM3HYTuple_cxx

#include "AnalyzeM3HYTuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>

void AnalyzeM3HYTuple::Loop(const string outFileName)
{
  //   In a ROOT session, you can do:
  //      root> .L AnalyzeM3HYTuple.C
  //      root> AnalyzeM3HYTuple t
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

  const double CSVM = 0.8484;
  const double CSVT = 0.9535;

  TFile* fout = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");
  tree->SetDirectory(fout);

  TH1F* hLW_m = new TH1F("hLW_m", "Leptonic W Mass", 250, 0, 500);
  TH1F* hLT_m = new TH1F("hLT_m", "Leptonic Top Mass", 250, 0, 500);
  TH1F* hHW_m = new TH1F("hHW_m", "Hadronic W Mass", 250, 0, 500);
  TH1F* hHW_dR = new TH1F("hHW_dR", "Hadronic W DeltaR", 100, 0, 5);
  TH1F* hHT_m = new TH1F("hHT_m", "Hadronic Top Mass", 250, 0, 500);
  TH1F* hAddJJ_dR = new TH1F("hAddJJ_dR", "hAddJJ_dR", 100, 0, 5);
  TH1F* hAddJJ_m = new TH1F("hAddJJ_m", "hAddJJ_m", 250, 0, 500);

  int b_run, b_event, b_vertex_n;
  int b_jets_n, b_bjets_n;
  float b_weight_gen;
  float b_lepton_pt, b_lepton_eta, b_lepton_phi;
  float b_met_pt, b_met_phi;
  int b_m3_bjetcode;
  float b_m3_lepB_pt, b_m3_lepB_eta, b_m3_lepB_phi, b_m3_lepB_m;
  float b_m3_lepW_pt, b_m3_lepW_eta, b_m3_lepW_phi, b_m3_lepW_m;
  float b_m3_lepT_pt, b_m3_lepT_eta, b_m3_lepT_phi, b_m3_lepT_m;
  float b_m3_hadJ1_pt, b_m3_hadJ1_eta, b_m3_hadJ1_phi, b_m3_hadJ1_m;
  float b_m3_hadJ2_pt, b_m3_hadJ2_eta, b_m3_hadJ2_phi, b_m3_hadJ2_m;
  float b_m3_hadB_pt, b_m3_hadB_eta, b_m3_hadB_phi, b_m3_hadB_m;
  float b_m3_hadW12_pt, b_m3_hadW12_eta, b_m3_hadW12_phi, b_m3_hadW12_m, b_m3_hadW12_dR;
  float b_m3_hadW23_pt, b_m3_hadW23_eta, b_m3_hadW23_phi, b_m3_hadW23_m, b_m3_hadW23_dR;
  float b_m3_hadT_pt, b_m3_hadT_eta, b_m3_hadT_phi, b_m3_hadT_m;
  float b_m3_theta1, b_m3_theta2;
  float b_m3_lepB_CSV, b_m3_hadB_CSV, b_m3_hadJ1_CSV, b_m3_hadJ2_CSV;
  float b_m3_lepB_CvsB, b_m3_hadB_CvsB, b_m3_hadJ1_CvsB, b_m3_hadJ2_CvsB;
  float b_m3_lepB_CvsL, b_m3_hadB_CvsL, b_m3_hadJ1_CvsL, b_m3_hadJ2_CvsL;
  float b_m3_addJetByPt1_pt, b_m3_addJetByPt1_CSV;
  float b_m3_addJetByPt2_pt, b_m3_addJetByPt2_CSV;
  float b_m3_addJetByCSV1_pt, b_m3_addJetByCSV1_CSV;
  float b_m3_addJetByCSV2_pt, b_m3_addJetByCSV2_CSV;
  float b_m3_addJetsByPt_m, b_m3_addJetsByPt_dR;
  float b_m3_addJetsByCSV_m, b_m3_addJetsByCSV_dR;

  tree->Branch("run", &b_run, "run/I");
  tree->Branch("event", &b_event, "event/I");
  tree->Branch("vertex_n", &b_vertex_n, "vertex_n/I");
  tree->Branch("weight_gen", &b_weight_gen, "weight_gen/F");
  tree->Branch("lepton_pt", &b_lepton_pt, "lepton_pt/F");
  tree->Branch("lepton_eta", &b_lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &b_lepton_phi, "lepton_phi/F");
  tree->Branch("met_pt", &b_met_pt, "met_pt/F");
  tree->Branch("met_phi", &b_met_phi, "met_phi/F");

  tree->Branch("jets_n", &b_jets_n, "jets_n/I");
  tree->Branch("bjets_n", &b_bjets_n, "bjets_n/I");
  tree->Branch("m3_lepW_m", &b_m3_lepW_m, "m3_lepW_m/F");
  tree->Branch("m3_lepT_m", &b_m3_lepT_m, "m3_lepT_m/F");
  tree->Branch("m3_hadW12_m", &b_m3_hadW12_m, "m3_hadW12_m/F");
  tree->Branch("m3_hadW23_m", &b_m3_hadW23_m, "m3_hadW23_m/F");
  tree->Branch("m3_hadT_m", &b_m3_hadT_m, "m3_hadT_m/F");

  tree->Branch("m3_bjetcode", &b_m3_bjetcode, "m3_bjetcode/I"); // b-jet contribution "code". Format=[nbjetInLepT][nbjetInHadT]
  tree->Branch("m3_lepB_pt", &b_m3_lepB_pt, "m3_lepB_pt/F");
  tree->Branch("m3_lepB_eta", &b_m3_lepB_eta, "m3_lepB_eta/F");
  tree->Branch("m3_lepB_phi", &b_m3_lepB_phi, "m3_lepB_phi/F");
  tree->Branch("m3_lepB_m", &b_m3_lepB_m, "m3_lepB_m/F");
  tree->Branch("m3_lepW_pt", &b_m3_lepW_pt, "m3_lepW_pt/F");
  tree->Branch("m3_lepW_eta", &b_m3_lepW_eta, "m3_lepW_eta/F");
  tree->Branch("m3_lepW_phi", &b_m3_lepW_phi, "m3_lepW_phi/F");
  tree->Branch("m3_lepT_pt", &b_m3_lepT_pt, "m3_lepT_pt/F");
  tree->Branch("m3_lepT_eta", &b_m3_lepT_eta, "m3_lepT_eta/F");
  tree->Branch("m3_lepT_phi", &b_m3_lepT_phi, "m3_lepT_phi/F");
  tree->Branch("m3_hadJ1_pt", &b_m3_hadJ1_pt, "m3_hadJ1_pt/F");
  tree->Branch("m3_hadJ1_eta", &b_m3_hadJ1_eta, "m3_hadJ1_eta/F");
  tree->Branch("m3_hadJ1_phi", &b_m3_hadJ1_phi, "m3_hadJ1_phi/F");
  tree->Branch("m3_hadJ1_m", &b_m3_hadJ1_m, "m3_hadJ1_m/F");
  tree->Branch("m3_hadJ2_pt", &b_m3_hadJ2_pt, "m3_hadJ2_pt/F");
  tree->Branch("m3_hadJ2_eta", &b_m3_hadJ2_eta, "m3_hadJ2_eta/F");
  tree->Branch("m3_hadJ2_phi", &b_m3_hadJ2_phi, "m3_hadJ2_phi/F");
  tree->Branch("m3_hadJ2_m", &b_m3_hadJ2_m, "m3_hadJ2_m/F");
  tree->Branch("m3_hadB_pt", &b_m3_hadB_pt, "m3_hadB_pt/F");
  tree->Branch("m3_hadB_eta", &b_m3_hadB_eta, "m3_hadB_eta/F");
  tree->Branch("m3_hadB_phi", &b_m3_hadB_phi, "m3_hadB_phi/F");
  tree->Branch("m3_hadB_m", &b_m3_hadB_m, "m3_hadB_m/F");
  tree->Branch("m3_hadW12_pt", &b_m3_hadW12_pt, "m3_hadW12_pt/F");
  tree->Branch("m3_hadW12_eta", &b_m3_hadW12_eta, "m3_hadW12_eta/F");
  tree->Branch("m3_hadW12_phi", &b_m3_hadW12_phi, "m3_hadW12_phi/F");
  tree->Branch("m3_hadW12_dR", &b_m3_hadW12_dR, "m3_hadW12_dR/F");
  tree->Branch("m3_hadW23_pt", &b_m3_hadW23_pt, "m3_hadW23_pt/F");
  tree->Branch("m3_hadW23_eta", &b_m3_hadW23_eta, "m3_hadW23_eta/F");
  tree->Branch("m3_hadW23_phi", &b_m3_hadW23_phi, "m3_hadW23_phi/F");
  tree->Branch("m3_hadW23_dR", &b_m3_hadW23_dR, "m3_hadW23_dR/F");
  tree->Branch("m3_hadT_pt", &b_m3_hadT_pt, "m3_hadT_pt/F");
  tree->Branch("m3_hadT_eta", &b_m3_hadT_eta, "m3_hadT_eta/F");
  tree->Branch("m3_hadT_phi", &b_m3_hadT_phi, "m3_hadT_phi/F");

  tree->Branch("m3_theta1", &b_m3_theta1, "m3_theta1/F"); // Angle between top and b
  tree->Branch("m3_theta2", &b_m3_theta2, "m3_theta2/F"); // Angle between t-b and w->jj plane

  tree->Branch("m3_lepB_CSV", &b_m3_lepB_CSV, "m3_lepB_CSV/F");
  tree->Branch("m3_hadB_CSV", &b_m3_hadB_CSV, "m3_hadB_CSV/F");
  tree->Branch("m3_hadJ1_CSV", &b_m3_hadJ1_CSV, "m3_hadJ1_CSV/F");
  tree->Branch("m3_hadJ2_CSV", &b_m3_hadJ2_CSV, "m3_hadJ2_CSV/F");

  tree->Branch("m3_lepB_CvsB", &b_m3_lepB_CvsB, "m3_lepB_CvsB/F");
  tree->Branch("m3_hadB_CvsB", &b_m3_hadB_CvsB, "m3_hadB_CvsB/F");
  tree->Branch("m3_hadJ1_CvsB", &b_m3_hadJ1_CvsB, "m3_hadJ1_CvsB/F");
  tree->Branch("m3_hadJ2_CvsB", &b_m3_hadJ2_CvsB, "m3_hadJ2_CvsB/F");

  tree->Branch("m3_lepB_CvsL", &b_m3_lepB_CvsL, "m3_lepB_CvsL/F");
  tree->Branch("m3_hadB_CvsL", &b_m3_hadB_CvsL, "m3_hadB_CvsL/F");
  tree->Branch("m3_hadJ1_CvsL", &b_m3_hadJ1_CvsL, "m3_hadJ1_CvsL/F");
  tree->Branch("m3_hadJ2_CvsL", &b_m3_hadJ2_CvsL, "m3_hadJ2_CvsL/F");

  tree->Branch("m3_addJetsByPt_m", &b_m3_addJetsByPt_m, "m3_addJetsByPt_m/F");
  tree->Branch("m3_addJetsByPt_dR", &b_m3_addJetsByPt_dR, "m3_addJetsByPt_dR/F");
  tree->Branch("m3_addJetsByCSV_m", &b_m3_addJetsByCSV_m, "m3_addJetsByCSV_m/F");
  tree->Branch("m3_addJetsByCSV_dR", &b_m3_addJetsByCSV_dR, "m3_addJetsByCSV_dR/F");

  tree->Branch("m3_addJetByPt1_pt", &b_m3_addJetByPt1_pt, "m3_addJetByPt1_pt/F");
  tree->Branch("m3_addJetByPt2_pt", &b_m3_addJetByPt2_pt, "m3_addJetByPt2_pt/F");
  tree->Branch("m3_addJetByPt1_CSV", &b_m3_addJetByPt1_CSV, "m3_addJetByPt1_CSV/F");
  tree->Branch("m3_addJetByPt2_CSV", &b_m3_addJetByPt2_CSV, "m3_addJetByPt2_CSV/F");

  tree->Branch("m3_addJetByCSV1_pt", &b_m3_addJetByCSV1_pt, "m3_addJetByCSV1_pt/F");
  tree->Branch("m3_addJetByCSV2_pt", &b_m3_addJetByCSV2_pt, "m3_addJetByCSV2_pt/F");
  tree->Branch("m3_addJetByCSV1_CSV", &b_m3_addJetByCSV1_CSV, "m3_addJetByCSV1_CSV/F");
  tree->Branch("m3_addJetByCSV2_CSV", &b_m3_addJetByCSV2_CSV, "m3_addJetByCSV2_CSV/F");

  if (fChain == 0) return;

  //Long64_t nentries = 10000;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    cout << jentry << '/' << nentries << '\r';

    b_run = run;
    b_event = event;
    b_vertex_n = GoodPV;
    b_weight_gen = genweight;

    TLorentzVector leptonP4;
    leptonP4.SetPtEtaPhiE(lepton_pT, lepton_eta, lepton_phi, lepton_E);
    if ( leptonP4.Pt() < 30 or std::abs(leptonP4.Eta()) > 2.1 ) continue;

    b_lepton_pt = lepton_pT;
    b_lepton_eta = lepton_eta;
    b_lepton_phi = lepton_phi;

    TLorentzVector metP4;
    metP4.SetPtEtaPhiM(MET, 0, MET_phi, 0);

    b_met_pt = MET;
    b_met_phi = MET_phi;

    std::vector<size_t> jetIdxs;
    const auto jets_pt = *jet_pT, jets_eta = *jet_eta, jets_phi = *jet_phi, jets_e = *jet_E;
    const auto jets_CSV = *jet_CSV, jets_CvsB = *jet_CvsB, jets_CvsL = *jet_CvsL;
    b_bjets_n = 0;
    for ( int j=0; j<jet_number; ++j ) {
      if ( jets_pt[j] < 30 or std::abs(jets_eta[j]) > 2.5 ) continue;
      if ( jets_CSV[j] > CSVM ) ++b_bjets_n;
      jetIdxs.push_back(j);
    }
    if ( jetIdxs.size() < 4 ) continue;
    b_jets_n = jetIdxs.size();

    double maxM3Pt = 0;
    std::vector<size_t> bestIdxs;
    TLorentzVector jetP4s[4];
    for ( auto ii1 = jetIdxs.begin(); ii1 != jetIdxs.end(); ++ii1 ) {
      jetP4s[1].SetPtEtaPhiE(jets_pt[*ii1], jets_eta[*ii1], jets_phi[*ii1], jets_e[*ii1]);
      for ( auto ii2 = ii1+1; ii2 != jetIdxs.end(); ++ii2 ) {
        jetP4s[2].SetPtEtaPhiE(jets_pt[*ii2], jets_eta[*ii2], jets_phi[*ii2], jets_e[*ii2]);
        for ( auto ii3 = ii2+1; ii3 != jetIdxs.end(); ++ii3 ) {
          jetP4s[3].SetPtEtaPhiE(jets_pt[*ii3], jets_eta[*ii3], jets_phi[*ii3], jets_e[*ii3]);

          const double m3Pt = (jetP4s[1]+jetP4s[2]+jetP4s[3]).Pt();
          if ( m3Pt > maxM3Pt ) {
            maxM3Pt = m3Pt;
            bestIdxs = {size_t(b_jets_n), *ii1, *ii2, *ii3};
          }
        }
      }
    }
    if ( bestIdxs.empty() ) continue;

    // Sort by CSV in increasing order - j1+j2 will prefer 80GeV for SM top, j2+j3 will prefer 125GeV for FCNC top
    // Keep j0 as is which is the jet from leptonically decaying top
    std::sort(std::next(bestIdxs.begin()), bestIdxs.end(),
              [&](size_t a, size_t b){ return jets_CSV[a] < jets_CSV[b]; });

    for ( auto i : jetIdxs ) {
      if ( i == bestIdxs[1] or i == bestIdxs[2] or i == bestIdxs[3] ) continue;
      if ( bestIdxs[0] == size_t(b_jets_n) or jets_pt[bestIdxs[0]] < jets_pt[i] ) {
        bestIdxs[0] = i;
      }
    }
    b_m3_bjetcode = 0;
    for ( size_t i=0; i<4; ++i ) {
      const size_t j = bestIdxs[i];
      jetP4s[i].SetPtEtaPhiE(jets_pt[j], jets_eta[j], jets_phi[j], jets_e[j]);
      if ( jets_CSV[j] > CSVM ) {
        if ( i == 0 ) b_m3_bjetcode = 10;
        else b_m3_bjetcode += 1;
      }
    }

    hLW_m->Fill( (leptonP4+metP4).M() );
    hLT_m->Fill( (leptonP4+metP4+jetP4s[0]).M() );
    hHW_m->Fill( (jetP4s[1]+jetP4s[2]).M() );
    hHW_dR->Fill( jetP4s[1].DeltaR(jetP4s[2]) );
    hHT_m->Fill( (jetP4s[1]+jetP4s[2]+jetP4s[3]).M() );

    b_m3_lepB_pt = jetP4s[0].Pt(); b_m3_lepB_eta = jetP4s[0].Eta(); b_m3_lepB_phi = jetP4s[0].Phi(); b_m3_lepB_m = jetP4s[0].M(); 
    b_m3_hadJ1_pt = jetP4s[1].Pt(); b_m3_hadJ1_eta = jetP4s[1].Eta(); b_m3_hadJ1_phi = jetP4s[1].Phi(); b_m3_hadJ1_m = jetP4s[1].M(); 
    b_m3_hadJ2_pt = jetP4s[1].Pt(); b_m3_hadJ2_eta = jetP4s[1].Eta(); b_m3_hadJ2_phi = jetP4s[1].Phi(); b_m3_hadJ2_m = jetP4s[1].M(); 
    b_m3_hadB_pt = jetP4s[2].Pt(); b_m3_hadB_eta = jetP4s[2].Eta(); b_m3_hadB_phi = jetP4s[2].Phi(); b_m3_hadB_m = jetP4s[2].M(); 

    const auto lepW = leptonP4+metP4;
    const auto lepT = lepW+jetP4s[0];
    b_m3_lepW_pt = lepW.Pt(); b_m3_lepW_eta = lepW.Eta(); b_m3_lepW_phi = lepW.Phi(); b_m3_lepW_m = lepW.M(); 
    b_m3_lepT_pt = lepT.Pt(); b_m3_lepT_eta = lepT.Eta(); b_m3_lepT_phi = lepT.Phi(); b_m3_lepT_m = lepT.M(); 

    const auto hadW12 = jetP4s[1]+jetP4s[2];
    const auto hadW23 = jetP4s[2]+jetP4s[3];
    const auto hadT = hadW12+jetP4s[3];
    b_m3_hadW12_pt = hadW12.Pt(); b_m3_hadW12_eta = hadW12.Eta(); b_m3_hadW12_phi = hadW12.Phi(); b_m3_hadW12_m = hadW12.M(); 
    b_m3_hadW23_pt = hadW23.Pt(); b_m3_hadW23_eta = hadW23.Eta(); b_m3_hadW23_phi = hadW23.Phi(); b_m3_hadW23_m = hadW23.M(); 
    b_m3_hadW12_dR = jetP4s[1].DeltaR(jetP4s[2]);
    b_m3_hadW23_dR = jetP4s[2].DeltaR(jetP4s[3]);
    b_m3_hadT_pt = hadT.Pt(); b_m3_hadT_eta = hadT.Eta(); b_m3_hadT_phi = hadT.Phi(); b_m3_hadT_m = hadT.M(); 

    TLorentzVector cm_hb = jetP4s[3], cm_hj1 = jetP4s[1], cm_hj2 = jetP4s[2];
    TLorentzVector cm_top = cm_hb+cm_hj1+cm_hj2;
    cm_hb.Boost(-cm_top.BoostVector());
    cm_hj1.Boost(-cm_top.BoostVector());
    cm_hj2.Boost(-cm_top.BoostVector());
    cm_hb *= 1./cm_hb.P();
    cm_hj1 *= 1./cm_hj1.P();
    cm_hj2 *= 1./cm_hj2.P();
    cm_top *= 1./cm_top.P();
    b_m3_theta1 = cm_hb.Vect().Dot(cm_top.Vect());
    b_m3_theta2 = cm_hb.Vect().Cross(cm_hj1.Vect()).Dot(cm_hb.Vect().Cross(cm_top.Vect()));

    b_m3_lepB_CSV = jets_CSV[bestIdxs[0]];
    b_m3_hadJ1_CSV = jets_CSV[bestIdxs[1]];
    b_m3_hadJ2_CSV = jets_CSV[bestIdxs[2]];
    b_m3_hadB_CSV = jets_CSV[bestIdxs[3]];

    b_m3_lepB_CvsB = jets_CvsB[bestIdxs[0]];
    b_m3_hadJ1_CvsB = jets_CvsB[bestIdxs[1]];
    b_m3_hadJ2_CvsB = jets_CvsB[bestIdxs[2]];
    b_m3_hadB_CvsB = jets_CvsB[bestIdxs[3]];

    b_m3_lepB_CvsL = jets_CvsL[bestIdxs[0]];
    b_m3_hadJ1_CvsL = jets_CvsL[bestIdxs[1]];
    b_m3_hadJ2_CvsL = jets_CvsL[bestIdxs[2]];
    b_m3_hadB_CvsL = jets_CvsL[bestIdxs[3]];

    std::vector<size_t> addJetIdxs;
    for ( size_t j : jetIdxs ) {
      bool isTaken = false;
      for ( size_t bestIdx : bestIdxs ) {
        if ( bestIdx == j ) { isTaken = true; break; }
      }
      if ( isTaken ) continue;

      addJetIdxs.push_back(j);
    }
    if ( addJetIdxs.size() < 2 ) {
      b_m3_addJetByPt1_pt = b_m3_addJetByPt2_pt = b_m3_addJetByCSV1_pt = b_m3_addJetByCSV2_pt = 0;
      b_m3_addJetByPt1_CSV = b_m3_addJetByPt2_CSV = b_m3_addJetByCSV1_CSV = b_m3_addJetByCSV2_CSV = -10;
      b_m3_addJetsByPt_m = b_m3_addJetsByCSV_m = 0;
      b_m3_addJetsByPt_dR = b_m3_addJetsByCSV_dR = 0;
    }
    else {
      auto addJetIdxsByPt = addJetIdxs;
      std::sort(addJetIdxsByPt.begin(), addJetIdxsByPt.end(), [&](size_t a, size_t b) { return jets_pt[a] > jets_pt[b]; });
      std::sort(addJetIdxs.begin(), addJetIdxs.end(), [&](size_t a, size_t b) { return jets_CSV[a] > jets_CSV[b]; });
      const size_t jByPt1 = addJetIdxsByPt[0], jByPt2 = addJetIdxsByPt[1];
      const size_t jByCSV1 = addJetIdxs[0], jByCSV2 = addJetIdxs[1];
      TLorentzVector addJetByCSV1, addJetByCSV2, addJetByPt1, addJetByPt2;
      addJetByCSV1.SetPtEtaPhiE(jets_pt[jByPt1], jets_eta[jByPt1], jets_phi[jByPt1], jets_e[jByPt1]);
      addJetByCSV1.SetPtEtaPhiE(jets_pt[jByPt2], jets_eta[jByPt2], jets_phi[jByPt2], jets_e[jByPt2]);
      addJetByPt1.SetPtEtaPhiE(jets_pt[jByCSV1], jets_eta[jByCSV1], jets_phi[jByCSV1], jets_e[jByCSV1]);
      addJetByPt1.SetPtEtaPhiE(jets_pt[jByCSV2], jets_eta[jByCSV2], jets_phi[jByCSV2], jets_e[jByCSV2]);

      b_m3_addJetByPt1_pt = jets_pt[jByPt1];
      b_m3_addJetByPt2_pt = jets_pt[jByPt2];
      b_m3_addJetByCSV1_pt = jets_pt[jByCSV1];
      b_m3_addJetByCSV2_pt = jets_pt[jByCSV2];
      b_m3_addJetByPt1_CSV = jets_CSV[jByPt1];
      b_m3_addJetByPt2_CSV = jets_CSV[jByPt2];
      b_m3_addJetByCSV1_CSV = jets_CSV[jByCSV1];
      b_m3_addJetByCSV2_CSV = jets_CSV[jByCSV2];
      b_m3_addJetsByPt_dR = addJetByPt1.DeltaR(addJetByPt2);
      b_m3_addJetsByCSV_dR = addJetByCSV1.DeltaR(addJetByCSV2);
      b_m3_addJetsByPt_m = (addJetByPt1+addJetByPt2).M();
      b_m3_addJetsByCSV_m = (addJetByCSV1+addJetByCSV2).M();

      hAddJJ_dR->Fill(b_m3_addJetsByCSV_dR);
      hAddJJ_m->Fill(b_m3_addJetsByCSV_m);

    }

    tree->Fill();
  }

  fout->Write();

  TCanvas* c = 0;
  c = new TCanvas("cAddJJ_dR", "addJJ dR", 500, 500); hAddJJ_dR->Draw();
  c = new TCanvas("cLW_m", "MW lep", 500, 500); hLW_m->Draw();
  c = new TCanvas("cLT_m", "MTop lep", 500, 500); hLT_m->Draw();
  c = new TCanvas("cHW_m", "MW had", 500, 500); hHW_m->Draw();
  c = new TCanvas("cHW_dR", "dRW had", 500, 500); hHW_dR->Draw();
  c = new TCanvas("cHT_m", "MTop had", 500, 500); hHT_m->Draw();
}
