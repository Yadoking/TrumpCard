#define AnalyzeM3Delphes_cxx

#include "AnalyzeM3Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>

void AnalyzeM3Delphes::Loop(const string outFileName)
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

  const double CSVM = 0.8484;
  const double CSVT = 0.9535;

  TFile* fout = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");
  tree->SetDirectory(fout);

  TH1F* hLW_m = new TH1F("hLW_m", "Leptonic W Mass", 250, 0, 500);
  TH1F* hLT_m = new TH1F("hLT_m", "Leptonic Top Mass", 250, 0, 500);
  TH1F* hHW_m = new TH1F("hHW_m", "Hadronic W Mass", 250, 0, 500);
  TH1F* hHW_dR = new TH1F("hHW_dR", "Hadronic W #DeltaR", 100, 0, 5);
  TH1F* hHT_m = new TH1F("hHT_m", "Hadronic Top Mass", 250, 0, 500);
  TH1F* hAddJJ_dR = new TH1F("hAddJJ_dR", "hAddJJ_dR", 100, 0, 5);
  TH1F* hAddJJ_m = new TH1F("hAddJJ_m", "hAddJJ_m", 250, 0, 500);

  int b_run, b_event, b_vertex_n;
  int b_jets_n, b_bjets_n;
  float b_weight_gen;
  float b_lepton_pt, b_lepton_eta, b_lepton_phi;
  float b_met_pt, b_met_phi;
  int b_kin_bjetcode;
  float b_kin_lepB_pt, b_kin_lepB_eta, b_kin_lepB_phi, b_kin_lepB_m;
  float b_kin_lepW_pt, b_kin_lepW_eta, b_kin_lepW_phi, b_kin_lepW_m;
  float b_kin_lepT_pt, b_kin_lepT_eta, b_kin_lepT_phi, b_kin_lepT_m;
  float b_kin_hadJ1_pt, b_kin_hadJ1_eta, b_kin_hadJ1_phi, b_kin_hadJ1_m;
  float b_kin_hadJ2_pt, b_kin_hadJ2_eta, b_kin_hadJ2_phi, b_kin_hadJ2_m;
  float b_kin_hadB_pt, b_kin_hadB_eta, b_kin_hadB_phi, b_kin_hadB_m;
  float b_kin_hadW12_pt, b_kin_hadW12_eta, b_kin_hadW12_phi, b_kin_hadW12_m, b_kin_hadW12_dR;
  float b_kin_hadW23_pt, b_kin_hadW23_eta, b_kin_hadW23_phi, b_kin_hadW23_m, b_kin_hadW23_dR;
  float b_kin_hadW13_pt, b_kin_hadW13_eta, b_kin_hadW13_phi, b_kin_hadW13_m, b_kin_hadW13_dR;
  float b_kin_hadT_pt, b_kin_hadT_eta, b_kin_hadT_phi, b_kin_hadT_m;
  float b_kin_theta1, b_kin_theta2;
  float b_kin_lepB_CSV, b_kin_hadB_CSV, b_kin_hadJ1_CSV, b_kin_hadJ2_CSV;
  float b_kin_lepB_CvsB, b_kin_hadB_CvsB, b_kin_hadJ1_CvsB, b_kin_hadJ2_CvsB;
  float b_kin_lepB_CvsL, b_kin_hadB_CvsL, b_kin_hadJ1_CvsL, b_kin_hadJ2_CvsL;
  float b_kin_addJetByPt1_pt, b_kin_addJetByPt1_CSV;
  float b_kin_addJetByPt2_pt, b_kin_addJetByPt2_CSV;
  float b_kin_addJetByCSV1_pt, b_kin_addJetByCSV1_CSV;
  float b_kin_addJetByCSV2_pt, b_kin_addJetByCSV2_CSV;
  float b_kin_addJetsByPt_m, b_kin_addJetsByPt_dR;
  float b_kin_addJetsByCSV_m, b_kin_addJetsByCSV_dR;

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
  tree->Branch("kin_lepW_m", &b_kin_lepW_m, "kin_lepW_m/F");
  tree->Branch("kin_lepT_m", &b_kin_lepT_m, "kin_lepT_m/F");
  tree->Branch("kin_hadW12_m", &b_kin_hadW12_m, "kin_hadW12_m/F");
  tree->Branch("kin_hadW23_m", &b_kin_hadW23_m, "kin_hadW23_m/F");
  tree->Branch("kin_hadW13_m", &b_kin_hadW13_m, "kin_hadW13_m/F");
  tree->Branch("kin_hadT_m", &b_kin_hadT_m, "kin_hadT_m/F");

  tree->Branch("kin_bjetcode", &b_kin_bjetcode, "kin_bjetcode/I"); // b-jet contribution "code". Format=[nbjetInLepT][nbjetInHadT]
  tree->Branch("kin_lepB_pt", &b_kin_lepB_pt, "kin_lepB_pt/F");
  tree->Branch("kin_lepB_eta", &b_kin_lepB_eta, "kin_lepB_eta/F");
  tree->Branch("kin_lepB_phi", &b_kin_lepB_phi, "kin_lepB_phi/F");
  tree->Branch("kin_lepB_m", &b_kin_lepB_m, "kin_lepB_m/F");
  tree->Branch("kin_lepW_pt", &b_kin_lepW_pt, "kin_lepW_pt/F");
  tree->Branch("kin_lepW_eta", &b_kin_lepW_eta, "kin_lepW_eta/F");
  tree->Branch("kin_lepW_phi", &b_kin_lepW_phi, "kin_lepW_phi/F");
  tree->Branch("kin_lepT_pt", &b_kin_lepT_pt, "kin_lepT_pt/F");
  tree->Branch("kin_lepT_eta", &b_kin_lepT_eta, "kin_lepT_eta/F");
  tree->Branch("kin_lepT_phi", &b_kin_lepT_phi, "kin_lepT_phi/F");
  tree->Branch("kin_hadJ1_pt", &b_kin_hadJ1_pt, "kin_hadJ1_pt/F");
  tree->Branch("kin_hadJ1_eta", &b_kin_hadJ1_eta, "kin_hadJ1_eta/F");
  tree->Branch("kin_hadJ1_phi", &b_kin_hadJ1_phi, "kin_hadJ1_phi/F");
  tree->Branch("kin_hadJ1_m", &b_kin_hadJ1_m, "kin_hadJ1_m/F");
  tree->Branch("kin_hadJ2_pt", &b_kin_hadJ2_pt, "kin_hadJ2_pt/F");
  tree->Branch("kin_hadJ2_eta", &b_kin_hadJ2_eta, "kin_hadJ2_eta/F");
  tree->Branch("kin_hadJ2_phi", &b_kin_hadJ2_phi, "kin_hadJ2_phi/F");
  tree->Branch("kin_hadJ2_m", &b_kin_hadJ2_m, "kin_hadJ2_m/F");
  tree->Branch("kin_hadB_pt", &b_kin_hadB_pt, "kin_hadB_pt/F");
  tree->Branch("kin_hadB_eta", &b_kin_hadB_eta, "kin_hadB_eta/F");
  tree->Branch("kin_hadB_phi", &b_kin_hadB_phi, "kin_hadB_phi/F");
  tree->Branch("kin_hadB_m", &b_kin_hadB_m, "kin_hadB_m/F");
  tree->Branch("kin_hadW12_pt", &b_kin_hadW12_pt, "kin_hadW12_pt/F");
  tree->Branch("kin_hadW12_eta", &b_kin_hadW12_eta, "kin_hadW12_eta/F");
  tree->Branch("kin_hadW12_phi", &b_kin_hadW12_phi, "kin_hadW12_phi/F");
  tree->Branch("kin_hadW12_dR", &b_kin_hadW12_dR, "kin_hadW12_dR/F");
  tree->Branch("kin_hadW23_pt", &b_kin_hadW23_pt, "kin_hadW23_pt/F");
  tree->Branch("kin_hadW23_eta", &b_kin_hadW23_eta, "kin_hadW23_eta/F");
  tree->Branch("kin_hadW23_phi", &b_kin_hadW23_phi, "kin_hadW23_phi/F");
  tree->Branch("kin_hadW23_dR", &b_kin_hadW23_dR, "kin_hadW23_dR/F");
  tree->Branch("kin_hadW13_pt", &b_kin_hadW13_pt, "kin_hadW13_pt/F");
  tree->Branch("kin_hadW13_eta", &b_kin_hadW13_eta, "kin_hadW13_eta/F");
  tree->Branch("kin_hadW13_phi", &b_kin_hadW13_phi, "kin_hadW13_phi/F");
  tree->Branch("kin_hadW13_dR", &b_kin_hadW13_dR, "kin_hadW13_dR/F");
  tree->Branch("kin_hadT_pt", &b_kin_hadT_pt, "kin_hadT_pt/F");
  tree->Branch("kin_hadT_eta", &b_kin_hadT_eta, "kin_hadT_eta/F");
  tree->Branch("kin_hadT_phi", &b_kin_hadT_phi, "kin_hadT_phi/F");

  tree->Branch("kin_theta1", &b_kin_theta1, "kin_theta1/F"); // Angle between top and b
  tree->Branch("kin_theta2", &b_kin_theta2, "kin_theta2/F"); // Angle between t-b and w->jj plane

  tree->Branch("kin_lepB_CSV", &b_kin_lepB_CSV, "kin_lepB_CSV/F");
  tree->Branch("kin_hadB_CSV", &b_kin_hadB_CSV, "kin_hadB_CSV/F");
  tree->Branch("kin_hadJ1_CSV", &b_kin_hadJ1_CSV, "kin_hadJ1_CSV/F");
  tree->Branch("kin_hadJ2_CSV", &b_kin_hadJ2_CSV, "kin_hadJ2_CSV/F");

  tree->Branch("kin_lepB_CvsB", &b_kin_lepB_CvsB, "kin_lepB_CvsB/F");
  tree->Branch("kin_hadB_CvsB", &b_kin_hadB_CvsB, "kin_hadB_CvsB/F");
  tree->Branch("kin_hadJ1_CvsB", &b_kin_hadJ1_CvsB, "kin_hadJ1_CvsB/F");
  tree->Branch("kin_hadJ2_CvsB", &b_kin_hadJ2_CvsB, "kin_hadJ2_CvsB/F");

  tree->Branch("kin_lepB_CvsL", &b_kin_lepB_CvsL, "kin_lepB_CvsL/F");
  tree->Branch("kin_hadB_CvsL", &b_kin_hadB_CvsL, "kin_hadB_CvsL/F");
  tree->Branch("kin_hadJ1_CvsL", &b_kin_hadJ1_CvsL, "kin_hadJ1_CvsL/F");
  tree->Branch("kin_hadJ2_CvsL", &b_kin_hadJ2_CvsL, "kin_hadJ2_CvsL/F");

  tree->Branch("kin_addJetsByPt_m", &b_kin_addJetsByPt_m, "kin_addJetsByPt_m/F");
  tree->Branch("kin_addJetsByPt_dR", &b_kin_addJetsByPt_dR, "kin_addJetsByPt_dR/F");
  tree->Branch("kin_addJetsByCSV_m", &b_kin_addJetsByCSV_m, "kin_addJetsByCSV_m/F");
  tree->Branch("kin_addJetsByCSV_dR", &b_kin_addJetsByCSV_dR, "kin_addJetsByCSV_dR/F");

  tree->Branch("kin_addJetByPt1_pt", &b_kin_addJetByPt1_pt, "kin_addJetByPt1_pt/F");
  tree->Branch("kin_addJetByPt2_pt", &b_kin_addJetByPt2_pt, "kin_addJetByPt2_pt/F");
  tree->Branch("kin_addJetByPt1_CSV", &b_kin_addJetByPt1_CSV, "kin_addJetByPt1_CSV/F");
  tree->Branch("kin_addJetByPt2_CSV", &b_kin_addJetByPt2_CSV, "kin_addJetByPt2_CSV/F");

  tree->Branch("kin_addJetByCSV1_pt", &b_kin_addJetByCSV1_pt, "kin_addJetByCSV1_pt/F");
  tree->Branch("kin_addJetByCSV2_pt", &b_kin_addJetByCSV2_pt, "kin_addJetByCSV2_pt/F");
  tree->Branch("kin_addJetByCSV1_CSV", &b_kin_addJetByCSV1_CSV, "kin_addJetByCSV1_CSV/F");
  tree->Branch("kin_addJetByCSV2_CSV", &b_kin_addJetByCSV2_CSV, "kin_addJetByCSV2_CSV/F");

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
    b_vertex_n = 1;
    b_weight_gen = 1;

    TLorentzVector leptonP4;
    if ( muons_n >= 1 and muons_pt[0] > 30 and std::abs(muons_eta[0]) < 2.1 ) {
      leptonP4.SetPtEtaPhiM(muons_pt[0], muons_eta[0], muons_phi[0], muons_m[0]);
    }
    else if ( electrons_n >= 1 and electrons_pt[0] > 30 and std::abs(electrons_eta[0]) < 2.1 ) {
      leptonP4.SetPtEtaPhiM(electrons_pt[0], electrons_eta[0], electrons_phi[0], electrons_m[0]);
    }
    else continue;

    b_lepton_pt = leptonP4.Pt();
    b_lepton_eta = leptonP4.Eta();
    b_lepton_phi = leptonP4.Phi();

    TLorentzVector metP4;
    metP4.SetPtEtaPhiM(met_pt, 0, met_phi, 0);

    std::vector<size_t> jetIdxs;
    b_bjets_n = 0;
    for ( size_t j=0; j<jets_n; ++j ) {
      if ( jets_pt[j] < 30 or std::abs(jets_eta[j]) > 2.5 ) continue;
      if ( jets_bTag[j] > CSVM ) ++b_bjets_n;
      jetIdxs.push_back(j);
    }
    if ( jetIdxs.size() < 4 ) continue;
    b_jets_n = jetIdxs.size();

    double maxM3Pt = 0;
    std::vector<size_t> bestIdxs;
    TLorentzVector jetP4s[4];
    for ( auto ii1 = jetIdxs.begin(); ii1 != jetIdxs.end(); ++ii1 ) {
      jetP4s[1].SetPtEtaPhiM(jets_pt[*ii1], jets_eta[*ii1], jets_phi[*ii1], jets_m[*ii1]);
      for ( auto ii2 = ii1+1; ii2 != jetIdxs.end(); ++ii2 ) {
        jetP4s[2].SetPtEtaPhiM(jets_pt[*ii2], jets_eta[*ii2], jets_phi[*ii2], jets_m[*ii2]);
        for ( auto ii3 = ii2+1; ii3 != jetIdxs.end(); ++ii3 ) {
          jetP4s[3].SetPtEtaPhiM(jets_pt[*ii3], jets_eta[*ii3], jets_phi[*ii3], jets_m[*ii3]);

          const double m3Pt = (jetP4s[1]+jetP4s[2]+jetP4s[3]).Pt();
          if ( m3Pt > maxM3Pt ) {
            maxM3Pt = m3Pt;
            bestIdxs = {jets_n, *ii1, *ii2, *ii3};
          }
        }
      }
    }
    if ( bestIdxs.empty() ) continue;

    // Sort by CSV in increasing order - j1+j2 will prefer 80GeV for SM top, j2+j3 will prefer 125GeV for FCNC top
    // Keep j0 as is which is the jet from leptonically decaying top
    std::sort(std::next(bestIdxs.begin()), bestIdxs.end(),
              [&](size_t a, size_t b){ return jets_bTag[a] < jets_bTag[b]; });
    std::stable_sort(std::next(bestIdxs.begin()), bestIdxs.end(),
                     [&](size_t a, size_t b){ return jets_pt[a] > jets_pt[b]; });

    for ( auto i : jetIdxs ) {
      if ( i == bestIdxs[1] or i == bestIdxs[2] or i == bestIdxs[3] ) continue;
      if ( bestIdxs[0] == jets_n or jets_pt[bestIdxs[0]] < jets_pt[i] ) {
        bestIdxs[0] = i;
      }
    }
    b_kin_bjetcode = 0;
    for ( size_t i=0; i<4; ++i ) {
      const size_t j = bestIdxs[i];
      jetP4s[i].SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
      if ( jets_bTag[j] > CSVM ) {
        if ( i == 0 ) b_kin_bjetcode = 10;
        else b_kin_bjetcode += 1;
      }
    }

    hLW_m->Fill( (leptonP4+metP4).M() );
    hLT_m->Fill( (leptonP4+metP4+jetP4s[0]).M() );
    hHW_m->Fill( (jetP4s[1]+jetP4s[2]).M() );
    hHW_dR->Fill( jetP4s[1].DeltaR(jetP4s[2]) );
    hHT_m->Fill( (jetP4s[1]+jetP4s[2]+jetP4s[3]).M() );

    b_kin_lepB_pt = jetP4s[0].Pt(); b_kin_lepB_eta = jetP4s[0].Eta(); b_kin_lepB_phi = jetP4s[0].Phi(); b_kin_lepB_m = jetP4s[0].M();
    b_kin_hadJ1_pt = jetP4s[1].Pt(); b_kin_hadJ1_eta = jetP4s[1].Eta(); b_kin_hadJ1_phi = jetP4s[1].Phi(); b_kin_hadJ1_m = jetP4s[1].M();
    b_kin_hadJ2_pt = jetP4s[2].Pt(); b_kin_hadJ2_eta = jetP4s[2].Eta(); b_kin_hadJ2_phi = jetP4s[2].Phi(); b_kin_hadJ2_m = jetP4s[2].M();
    b_kin_hadB_pt = jetP4s[3].Pt(); b_kin_hadB_eta = jetP4s[3].Eta(); b_kin_hadB_phi = jetP4s[3].Phi(); b_kin_hadB_m = jetP4s[3].M();

    const auto lepW = leptonP4+metP4;
    const auto lepT = lepW+jetP4s[0];
    b_kin_lepW_pt = lepW.Pt(); b_kin_lepW_eta = lepW.Eta(); b_kin_lepW_phi = lepW.Phi(); b_kin_lepW_m = lepW.M();
    b_kin_lepT_pt = lepT.Pt(); b_kin_lepT_eta = lepT.Eta(); b_kin_lepT_phi = lepT.Phi(); b_kin_lepT_m = lepT.M();

    const auto hadW12 = jetP4s[1]+jetP4s[2];
    const auto hadW23 = jetP4s[2]+jetP4s[3];
    const auto hadW13 = jetP4s[1]+jetP4s[3];
    const auto hadT = hadW12+jetP4s[3];
    b_kin_hadW12_pt = hadW12.Pt(); b_kin_hadW12_eta = hadW12.Eta(); b_kin_hadW12_phi = hadW12.Phi(); b_kin_hadW12_m = hadW12.M();
    b_kin_hadW23_pt = hadW23.Pt(); b_kin_hadW23_eta = hadW23.Eta(); b_kin_hadW23_phi = hadW23.Phi(); b_kin_hadW23_m = hadW23.M();
    b_kin_hadW13_pt = hadW13.Pt(); b_kin_hadW13_eta = hadW13.Eta(); b_kin_hadW13_phi = hadW13.Phi(); b_kin_hadW13_m = hadW13.M();
    b_kin_hadW12_dR = jetP4s[1].DeltaR(jetP4s[2]);
    b_kin_hadW23_dR = jetP4s[2].DeltaR(jetP4s[3]);
    b_kin_hadW13_dR = jetP4s[1].DeltaR(jetP4s[3]);
    b_kin_hadT_pt = hadT.Pt(); b_kin_hadT_eta = hadT.Eta(); b_kin_hadT_phi = hadT.Phi(); b_kin_hadT_m = hadT.M();

    TLorentzVector cm_hb = jetP4s[3], cm_hj1 = jetP4s[1], cm_hj2 = jetP4s[2];
    TLorentzVector cm_top = cm_hb+cm_hj1+cm_hj2;
    cm_hb.Boost(-cm_top.BoostVector());
    cm_hj1.Boost(-cm_top.BoostVector());
    cm_hj2.Boost(-cm_top.BoostVector());
    cm_hb *= 1./cm_hb.P();
    cm_hj1 *= 1./cm_hj1.P();
    cm_hj2 *= 1./cm_hj2.P();
    cm_top *= 1./cm_top.P();
    b_kin_theta1 = cm_hb.Vect().Dot(cm_top.Vect());
    b_kin_theta2 = cm_hb.Vect().Cross(cm_hj1.Vect()).Dot(cm_hb.Vect().Cross(cm_top.Vect()));

    b_kin_lepB_CSV = jets_bTag[bestIdxs[0]];
    b_kin_hadJ1_CSV = jets_bTag[bestIdxs[1]];
    b_kin_hadJ2_CSV = jets_bTag[bestIdxs[2]];
    b_kin_hadB_CSV = jets_bTag[bestIdxs[3]];

    b_kin_lepB_CvsB = 0;
    b_kin_hadJ1_CvsB = 0;
    b_kin_hadJ2_CvsB = 0;
    b_kin_hadB_CvsB = 0;

    b_kin_lepB_CvsL = 0;
    b_kin_hadJ1_CvsL = 0;
    b_kin_hadJ2_CvsL = 0;
    b_kin_hadB_CvsL = 0;

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
      b_kin_addJetByPt1_pt = b_kin_addJetByPt2_pt = b_kin_addJetByCSV1_pt = b_kin_addJetByCSV2_pt = 0;
      b_kin_addJetByPt1_CSV = b_kin_addJetByPt2_CSV = b_kin_addJetByCSV1_CSV = b_kin_addJetByCSV2_CSV = -10;
      b_kin_addJetsByPt_m = b_kin_addJetsByCSV_m = 0;
      b_kin_addJetsByPt_dR = b_kin_addJetsByCSV_dR = 0;
    }
    else {
      auto addJetIdxsByPt = addJetIdxs;
      std::sort(addJetIdxsByPt.begin(), addJetIdxsByPt.end(), [&](size_t a, size_t b) { return jets_pt[a] > jets_pt[b]; });
      std::sort(addJetIdxs.begin(), addJetIdxs.end(), [&](size_t a, size_t b) { return jets_bTag[a] > jets_bTag[b]; });
      const size_t jByPt1 = addJetIdxsByPt[0], jByPt2 = addJetIdxsByPt[1];
      const size_t jByCSV1 = addJetIdxs[0], jByCSV2 = addJetIdxs[1];
      TLorentzVector addJetByCSV1, addJetByCSV2, addJetByPt1, addJetByPt2;
      addJetByCSV1.SetPtEtaPhiM(jets_pt[jByPt1], jets_eta[jByPt1], jets_phi[jByPt1], jets_m[jByPt1]);
      addJetByCSV1.SetPtEtaPhiM(jets_pt[jByPt2], jets_eta[jByPt2], jets_phi[jByPt2], jets_m[jByPt2]);
      addJetByPt1.SetPtEtaPhiM(jets_pt[jByCSV1], jets_eta[jByCSV1], jets_phi[jByCSV1], jets_m[jByCSV1]);
      addJetByPt1.SetPtEtaPhiM(jets_pt[jByCSV2], jets_eta[jByCSV2], jets_phi[jByCSV2], jets_m[jByCSV2]);

      b_kin_addJetByPt1_pt = jets_pt[jByPt1];
      b_kin_addJetByPt2_pt = jets_pt[jByPt2];
      b_kin_addJetByCSV1_pt = jets_pt[jByCSV1];
      b_kin_addJetByCSV2_pt = jets_pt[jByCSV2];
      b_kin_addJetByPt1_CSV = jets_bTag[jByPt1];
      b_kin_addJetByPt2_CSV = jets_bTag[jByPt2];
      b_kin_addJetByCSV1_CSV = jets_bTag[jByCSV1];
      b_kin_addJetByCSV2_CSV = jets_bTag[jByCSV2];
      b_kin_addJetsByPt_dR = addJetByPt1.DeltaR(addJetByPt2);
      b_kin_addJetsByCSV_dR = addJetByCSV1.DeltaR(addJetByCSV2);
      b_kin_addJetsByPt_m = (addJetByPt1+addJetByPt2).M();
      b_kin_addJetsByCSV_m = (addJetByCSV1+addJetByCSV2).M();

      hAddJJ_dR->Fill(b_kin_addJetsByCSV_dR);
      hAddJJ_m->Fill(b_kin_addJetsByCSV_m);

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
