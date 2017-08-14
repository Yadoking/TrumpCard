#define AnalyzeM3Delphes_cxx

#include "AnalyzeM3Delphes.h"
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>

void AnalyzeM3Delphes::Loop(const string modeStr, const string outFileName, string algoTypeName)
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
  enum class Mode { TT=0, FCNC };
  Mode mode = (modeStr == "FCNC") ? Mode::FCNC : Mode::TT;
  std::transform(algoTypeName.begin(), algoTypeName.end(), algoTypeName.begin(),
                 [](unsigned char c) -> unsigned char { return std::toupper(c); });
  enum class AlgoType { M3=0, dR };
  AlgoType algoType = (algoTypeName == "M3" ? AlgoType::M3 : AlgoType::dR);

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
  float b_met_pt, b_met_dphi;
  int b_bjetcode;
  float b_lepB_pt, b_lepB_eta, b_lepB_dphi, b_lepB_m;
  float b_lepW_pt, b_lepW_eta, b_lepW_dphi, b_lepW_m;
  float b_lepT_pt, b_lepT_eta, b_lepT_dphi, b_lepT_m;
  float b_hadJ1_pt, b_hadJ1_eta, b_hadJ1_dphi, b_hadJ1_m;
  float b_hadJ2_pt, b_hadJ2_eta, b_hadJ2_dphi, b_hadJ2_m;
  float b_hadB_pt, b_hadB_eta, b_hadB_dphi, b_hadB_m;
  float b_hadW12_pt, b_hadW12_eta, b_hadW12_dphi, b_hadW12_m, b_hadW12_dR;
  float b_hadW23_pt, b_hadW23_eta, b_hadW23_dphi, b_hadW23_m, b_hadW23_dR;
  float b_hadW13_pt, b_hadW13_eta, b_hadW13_dphi, b_hadW13_m, b_hadW13_dR;
  float b_hadT_pt, b_hadT_eta, b_hadT_dphi, b_hadT_m;
  float b_theta1, b_theta2;
  float b_lepB_CSV, b_hadB_CSV, b_hadJ1_CSV, b_hadJ2_CSV;
  float b_lepB_CvsB, b_hadB_CvsB, b_hadJ1_CvsB, b_hadJ2_CvsB;
  float b_lepB_CvsL, b_hadB_CvsL, b_hadJ1_CvsL, b_hadJ2_CvsL;
  float b_addJetByPt1_pt, b_addJetByPt1_CSV;
  float b_addJetByPt2_pt, b_addJetByPt2_CSV;
  float b_addJetByCSV1_pt, b_addJetByCSV1_CSV;
  float b_addJetByCSV2_pt, b_addJetByCSV2_CSV;
  float b_addJetsByPt_m, b_addJetsByPt_dR;
  float b_addJetsByCSV_m, b_addJetsByCSV_dR;
  int b_genMatch; // [lep][nu][lepB][hadJ1][hadJ2][hadB][addJ1][addJ2]

  TH2F* b_hJetImage_ch_n  = new TH2F("hJetImage_ch_n", "Jet image ch n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* b_hJetImage_nh_n  = new TH2F("hJetImage_nh_n", "Jet image nh n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* b_hJetImage_ph_n  = new TH2F("hJetImage_ph_n", "Jet image ph n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* b_hJetImage_ch_pt = new TH2F("hJetImage_ch_pt", "Jet image ch pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* b_hJetImage_nh_pt = new TH2F("hJetImage_nh_pt", "Jet image nh pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* b_hJetImage_ph_pt = new TH2F("hJetImage_ph_pt", "Jet image ph pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  b_hJetImage_ch_n->SetDirectory(0);
  b_hJetImage_nh_n->SetDirectory(0);
  b_hJetImage_ph_n->SetDirectory(0);
  b_hJetImage_ch_pt->SetDirectory(0);
  b_hJetImage_nh_pt->SetDirectory(0);
  b_hJetImage_ph_pt->SetDirectory(0);

  tree->Branch("run", &b_run, "run/I");
  tree->Branch("event", &b_event, "event/I");
  tree->Branch("vertex_n", &b_vertex_n, "vertex_n/I");
  tree->Branch("weight_gen", &b_weight_gen, "weight_gen/F");
  tree->Branch("lepton_pt", &b_lepton_pt, "lepton_pt/F");
  tree->Branch("lepton_eta", &b_lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &b_lepton_phi, "lepton_phi/F");
  tree->Branch("met_pt", &b_met_pt, "met_pt/F");
  tree->Branch("met_dphi", &b_met_dphi, "met_dphi/F");

  tree->Branch("jets_n", &b_jets_n, "jets_n/I");
  tree->Branch("bjets_n", &b_bjets_n, "bjets_n/I");
  tree->Branch("lepW_m", &b_lepW_m, "lepW_m/F");
  tree->Branch("lepT_m", &b_lepT_m, "lepT_m/F");
  tree->Branch("hadW12_m", &b_hadW12_m, "hadW12_m/F");
  tree->Branch("hadW23_m", &b_hadW23_m, "hadW23_m/F");
  tree->Branch("hadW13_m", &b_hadW13_m, "hadW13_m/F");
  tree->Branch("hadT_m", &b_hadT_m, "hadT_m/F");

  tree->Branch("bjetcode", &b_bjetcode, "bjetcode/I"); // b-jet contribution "code". Format=[nbjetInLepT][nbjetInHadT]
  tree->Branch("lepB_pt", &b_lepB_pt, "lepB_pt/F");
  tree->Branch("lepB_eta", &b_lepB_eta, "lepB_eta/F");
  tree->Branch("lepB_dphi", &b_lepB_dphi, "lepB_dphi/F");
  tree->Branch("lepB_m", &b_lepB_m, "lepB_m/F");
  tree->Branch("lepW_pt", &b_lepW_pt, "lepW_pt/F");
  tree->Branch("lepW_eta", &b_lepW_eta, "lepW_eta/F");
  tree->Branch("lepW_dphi", &b_lepW_dphi, "lepW_dphi/F");
  tree->Branch("lepT_pt", &b_lepT_pt, "lepT_pt/F");
  tree->Branch("lepT_eta", &b_lepT_eta, "lepT_eta/F");
  tree->Branch("lepT_dphi", &b_lepT_dphi, "lepT_dphi/F");
  tree->Branch("hadJ1_pt", &b_hadJ1_pt, "hadJ1_pt/F");
  tree->Branch("hadJ1_eta", &b_hadJ1_eta, "hadJ1_eta/F");
  tree->Branch("hadJ1_dphi", &b_hadJ1_dphi, "hadJ1_dphi/F");
  tree->Branch("hadJ1_m", &b_hadJ1_m, "hadJ1_m/F");
  tree->Branch("hadJ2_pt", &b_hadJ2_pt, "hadJ2_pt/F");
  tree->Branch("hadJ2_eta", &b_hadJ2_eta, "hadJ2_eta/F");
  tree->Branch("hadJ2_dphi", &b_hadJ2_dphi, "hadJ2_dphi/F");
  tree->Branch("hadJ2_m", &b_hadJ2_m, "hadJ2_m/F");
  tree->Branch("hadB_pt", &b_hadB_pt, "hadB_pt/F");
  tree->Branch("hadB_eta", &b_hadB_eta, "hadB_eta/F");
  tree->Branch("hadB_dphi", &b_hadB_dphi, "hadB_dphi/F");
  tree->Branch("hadB_m", &b_hadB_m, "hadB_m/F");
  tree->Branch("hadW12_pt", &b_hadW12_pt, "hadW12_pt/F");
  tree->Branch("hadW12_eta", &b_hadW12_eta, "hadW12_eta/F");
  tree->Branch("hadW12_dphi", &b_hadW12_dphi, "hadW12_dphi/F");
  tree->Branch("hadW12_dR", &b_hadW12_dR, "hadW12_dR/F");
  tree->Branch("hadW23_pt", &b_hadW23_pt, "hadW23_pt/F");
  tree->Branch("hadW23_eta", &b_hadW23_eta, "hadW23_eta/F");
  tree->Branch("hadW23_dphi", &b_hadW23_dphi, "hadW23_dphi/F");
  tree->Branch("hadW23_dR", &b_hadW23_dR, "hadW23_dR/F");
  tree->Branch("hadW13_pt", &b_hadW13_pt, "hadW13_pt/F");
  tree->Branch("hadW13_eta", &b_hadW13_eta, "hadW13_eta/F");
  tree->Branch("hadW13_dphi", &b_hadW13_dphi, "hadW13_dphi/F");
  tree->Branch("hadW13_dR", &b_hadW13_dR, "hadW13_dR/F");
  tree->Branch("hadT_pt", &b_hadT_pt, "hadT_pt/F");
  tree->Branch("hadT_eta", &b_hadT_eta, "hadT_eta/F");
  tree->Branch("hadT_dphi", &b_hadT_dphi, "hadT_dphi/F");

  tree->Branch("theta1", &b_theta1, "theta1/F"); // Angle between top and b
  tree->Branch("theta2", &b_theta2, "theta2/F"); // Angle between t-b and w->jj plane

  tree->Branch("lepB_CSV", &b_lepB_CSV, "lepB_CSV/F");
  tree->Branch("hadB_CSV", &b_hadB_CSV, "hadB_CSV/F");
  tree->Branch("hadJ1_CSV", &b_hadJ1_CSV, "hadJ1_CSV/F");
  tree->Branch("hadJ2_CSV", &b_hadJ2_CSV, "hadJ2_CSV/F");

  tree->Branch("lepB_CvsB", &b_lepB_CvsB, "lepB_CvsB/F");
  tree->Branch("hadB_CvsB", &b_hadB_CvsB, "hadB_CvsB/F");
  tree->Branch("hadJ1_CvsB", &b_hadJ1_CvsB, "hadJ1_CvsB/F");
  tree->Branch("hadJ2_CvsB", &b_hadJ2_CvsB, "hadJ2_CvsB/F");

  tree->Branch("lepB_CvsL", &b_lepB_CvsL, "lepB_CvsL/F");
  tree->Branch("hadB_CvsL", &b_hadB_CvsL, "hadB_CvsL/F");
  tree->Branch("hadJ1_CvsL", &b_hadJ1_CvsL, "hadJ1_CvsL/F");
  tree->Branch("hadJ2_CvsL", &b_hadJ2_CvsL, "hadJ2_CvsL/F");

  tree->Branch("addJetsByPt_m", &b_addJetsByPt_m, "addJetsByPt_m/F");
  tree->Branch("addJetsByPt_dR", &b_addJetsByPt_dR, "addJetsByPt_dR/F");
  tree->Branch("addJetsByCSV_m", &b_addJetsByCSV_m, "addJetsByCSV_m/F");
  tree->Branch("addJetsByCSV_dR", &b_addJetsByCSV_dR, "addJetsByCSV_dR/F");

  tree->Branch("addJetByPt1_pt", &b_addJetByPt1_pt, "addJetByPt1_pt/F");
  tree->Branch("addJetByPt2_pt", &b_addJetByPt2_pt, "addJetByPt2_pt/F");
  tree->Branch("addJetByPt1_CSV", &b_addJetByPt1_CSV, "addJetByPt1_CSV/F");
  tree->Branch("addJetByPt2_CSV", &b_addJetByPt2_CSV, "addJetByPt2_CSV/F");

  tree->Branch("addJetByCSV1_pt", &b_addJetByCSV1_pt, "addJetByCSV1_pt/F");
  tree->Branch("addJetByCSV2_pt", &b_addJetByCSV2_pt, "addJetByCSV2_pt/F");
  tree->Branch("addJetByCSV1_CSV", &b_addJetByCSV1_CSV, "addJetByCSV1_CSV/F");
  tree->Branch("addJetByCSV2_CSV", &b_addJetByCSV2_CSV, "addJetByCSV2_CSV/F");

  tree->Branch("genMatch", &b_genMatch, "genMatch/I");

  tree->Branch("hJetImage_ch_n", "TH2F", b_hJetImage_ch_n);
  tree->Branch("hJetImage_nh_n", "TH2F", b_hJetImage_nh_n);
  tree->Branch("hJetImage_ph_n", "TH2F", b_hJetImage_ph_n);
  tree->Branch("hJetImage_ch_pt", "TH2F", b_hJetImage_ch_pt);
  tree->Branch("hJetImage_nh_pt", "TH2F", b_hJetImage_nh_pt);
  tree->Branch("hJetImage_ph_pt", "TH2F", b_hJetImage_ph_pt);

  TH2F* hJetImage_ch_n  = new TH2F("hJetImage_ch_n", "Jet image ch n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* hJetImage_nh_n  = new TH2F("hJetImage_nh_n", "Jet image nh n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* hJetImage_ph_n  = new TH2F("hJetImage_ph_n", "Jet image ph n;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* hJetImage_ch_pt = new TH2F("hJetImage_ch_pt", "Jet image ch pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* hJetImage_nh_pt = new TH2F("hJetImage_nh_pt", "Jet image nh pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);
  TH2F* hJetImage_ph_pt = new TH2F("hJetImage_ph_pt", "Jet image ph pt;#Delta#eta;#Delta#phi", 50, -2, 2, 50, -2, 2);

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
    b_hJetImage_ch_n->Reset();
    b_hJetImage_ph_n->Reset();
    b_hJetImage_nh_n->Reset();
    b_hJetImage_ch_pt->Reset();
    b_hJetImage_ph_pt->Reset();
    b_hJetImage_nh_pt->Reset();

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
    b_met_pt = met_pt;
    b_met_dphi = met_phi;

    std::vector<size_t> jetIdxs;
    b_bjets_n = 0;
    for ( size_t j=0; j<jets_n; ++j ) {
      if ( jets_pt[j] < 30 or std::abs(jets_eta[j]) > 2.5 ) continue;
      if ( jets_bTag[j] > CSVM ) ++b_bjets_n;
      jetIdxs.push_back(j);
    }
    if ( jetIdxs.size() < 4 ) continue;
    b_jets_n = jetIdxs.size();

    std::vector<size_t> bestIdxs;
    TLorentzVector jetP4s[4];
    if ( algoType == AlgoType::M3 ) {
      double maxM3Pt = 0;
      for ( auto ii1 = jetIdxs.begin(); ii1 != jetIdxs.end(); ++ii1 ) {
        jetP4s[1].SetPtEtaPhiM(jets_pt[*ii1], jets_eta[*ii1], jets_phi[*ii1], jets_m[*ii1]);
        for ( auto ii2 = ii1+1; ii2 != jetIdxs.end(); ++ii2 ) {
          jetP4s[2].SetPtEtaPhiM(jets_pt[*ii2], jets_eta[*ii2], jets_phi[*ii2], jets_m[*ii2]);
          for ( auto ii3 = ii2+1; ii3 != jetIdxs.end(); ++ii3 ) {
            int nbjetsInHadT = 0;
            if ( jets_bTag[*ii1] > CSVM ) ++nbjetsInHadT;
            if ( jets_bTag[*ii2] > CSVM ) ++nbjetsInHadT;
            if ( jets_bTag[*ii3] > CSVM ) ++nbjetsInHadT;
            if ( mode == Mode::TT ) {
              // SM ttbar mode: require at least one b jet used
              if ( b_bjets_n >= 2 and nbjetsInHadT < 1 ) continue;
            }
            else if ( mode == Mode::FCNC ) {
              // FCNC mode: require b jets in the trijet system, t->Hc / H->bb
              // can be used for the charged Higgs case as well, t->H+b / H+ -> cb
              if ( b_bjets_n >= 3 and nbjetsInHadT < 2 ) continue; // at least two b jets in hadronic side
              else if ( b_bjets_n == 2 and nbjetsInHadT < 1 ) continue; // at least one b jet in hadronic side
            }

            jetP4s[3].SetPtEtaPhiM(jets_pt[*ii3], jets_eta[*ii3], jets_phi[*ii3], jets_m[*ii3]);

            const double m3Pt = (jetP4s[1]+jetP4s[2]+jetP4s[3]).Pt();
            if ( m3Pt > maxM3Pt ) {
              maxM3Pt = m3Pt;
              bestIdxs = {jets_n, *ii1, *ii2, *ii3};
            }
          }
        }
      }
    }
    else if ( algoType == AlgoType::dR ) {
      double minDR = 1e9;
      for ( auto ii1 = jetIdxs.begin(); ii1 != jetIdxs.end(); ++ii1 ) {
        jetP4s[1].SetPtEtaPhiM(jets_pt[*ii1], jets_eta[*ii1], jets_phi[*ii1], jets_m[*ii1]);
        for ( auto ii2 = ii1+1; ii2 != jetIdxs.end(); ++ii2 ) {
          int nbjetsInHadW = 0;
          if ( jets_bTag[*ii1] > CSVM ) ++nbjetsInHadW;
          if ( jets_bTag[*ii2] > CSVM ) ++nbjetsInHadW;
          if ( mode == Mode::FCNC ) {
            // FCNC mode: require b jets in the trijet system, t->Hc / H->bb
            // can be used for the charged Higgs case as well, t->H+b / H+ -> cb
            if ( b_bjets_n >= 3 and nbjetsInHadW < 2 ) continue; // at least two b jets in hadronic side
            else if ( b_bjets_n == 2 and nbjetsInHadW < 1 ) continue; // at least one b jet in hadronic side
          }

          jetP4s[2].SetPtEtaPhiM(jets_pt[*ii2], jets_eta[*ii2], jets_phi[*ii2], jets_m[*ii2]);
          const double dR = jetP4s[1].DeltaR(jetP4s[2]);
          if ( dR < minDR ) {
            bestIdxs = {jets_n, *ii1, *ii2, jets_n};
            minDR = dR;
          }
        }
      }
      if ( !bestIdxs.empty() ) {
        const auto i1 = bestIdxs[1], i2 = bestIdxs[2];
        jetP4s[1].SetPtEtaPhiM(jets_pt[i1], jets_eta[i1], jets_phi[i1], jets_m[i1]);
        jetP4s[2].SetPtEtaPhiM(jets_pt[i2], jets_eta[i2], jets_phi[i2], jets_m[i2]);
        const auto wP4 = jetP4s[1]+jetP4s[2];
        double minDR2 = 1e9;
        for ( auto i3 : jetIdxs ) {
          if ( i3 == i1 or i3 == i2 ) continue;
          jetP4s[3].SetPtEtaPhiM(jets_pt[i3], jets_eta[i3], jets_phi[i3], jets_m[i3]);

          const double dR = jetP4s[3].DeltaR(wP4);
          if ( dR < minDR2 ) {
            bestIdxs[3] = i3;
            minDR2 = dR;
          }
        }
        if ( minDR2 == 1e9 ) bestIdxs.clear();
      }
    }
    if ( bestIdxs.empty() ) continue;

    // Sort by CSV in increasing order - j1+j2 will prefer 80GeV for SM top, j2+j3 will prefer 125GeV for FCNC top
    // Keep j0 as is which is the jet from leptonically decaying top
    std::sort(std::next(bestIdxs.begin()), bestIdxs.end(),
              [&](size_t a, size_t b){ return jets_pt[a] > jets_pt[b]; });
    if ( mode == Mode::TT ) {
      std::stable_sort(std::next(bestIdxs.begin()), bestIdxs.end(),
                       [&](size_t a, size_t b){ return jets_bTag[a] < jets_bTag[b]; });
    }
    else if ( mode == Mode::FCNC ) {
      std::stable_sort(std::next(bestIdxs.begin()), bestIdxs.end(),
                       [&](size_t a, size_t b){ return jets_bTag[a] > jets_bTag[b]; });
    }

    for ( auto i : jetIdxs ) {
      if ( i == bestIdxs[1] or i == bestIdxs[2] or i == bestIdxs[3] ) continue;
      if ( bestIdxs[0] == jets_n or jets_pt[bestIdxs[0]] < jets_pt[i] ) {
        bestIdxs[0] = i;
      }
    }
    b_bjetcode = 0;
    for ( size_t i=0; i<4; ++i ) {
      const size_t j = bestIdxs[i];
      jetP4s[i].SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
      if ( jets_bTag[j] > CSVM ) {
        if ( i == 0 ) b_bjetcode = 10;
        else b_bjetcode += 1;
      }
    }

    // Generator level objects
    b_genMatch = 0; // [lep][nu][lepB][hadJ1][hadJ2][hadB][addJ1][addJ2]
    TLorentzVector gen_lep, gen_nu, gen_lepB, gen_hadJ1, gen_hadJ2, gen_hadB;
    for ( size_t i=0; i<b_gen_n; ++i ) {
    }

    hLW_m->Fill( (leptonP4+metP4).M() );
    hLT_m->Fill( (leptonP4+metP4+jetP4s[0]).M() );
    hHW_m->Fill( (jetP4s[1]+jetP4s[2]).M() );
    hHW_dR->Fill( jetP4s[1].DeltaR(jetP4s[2]) );
    hHT_m->Fill( (jetP4s[1]+jetP4s[2]+jetP4s[3]).M() );

    b_lepB_pt = jetP4s[0].Pt(); b_lepB_eta = jetP4s[0].Eta(); b_lepB_dphi = jetP4s[0].Phi(); b_lepB_m = jetP4s[0].M();
    b_hadJ1_pt = jetP4s[1].Pt(); b_hadJ1_eta = jetP4s[1].Eta(); b_hadJ1_dphi = jetP4s[1].Phi(); b_hadJ1_m = jetP4s[1].M();
    b_hadJ2_pt = jetP4s[2].Pt(); b_hadJ2_eta = jetP4s[2].Eta(); b_hadJ2_dphi = jetP4s[2].Phi(); b_hadJ2_m = jetP4s[2].M();
    b_hadB_pt = jetP4s[3].Pt(); b_hadB_eta = jetP4s[3].Eta(); b_hadB_dphi = jetP4s[3].Phi(); b_hadB_m = jetP4s[3].M();

    const auto lepW = leptonP4+metP4;
    const auto lepT = lepW+jetP4s[0];
    b_lepW_pt = lepW.Pt(); b_lepW_eta = lepW.Eta(); b_lepW_dphi = lepW.Phi(); b_lepW_m = lepW.M();
    b_lepT_pt = lepT.Pt(); b_lepT_eta = lepT.Eta(); b_lepT_dphi = lepT.Phi(); b_lepT_m = lepT.M();

    const auto hadW12 = jetP4s[1]+jetP4s[2];
    const auto hadW23 = jetP4s[2]+jetP4s[3];
    const auto hadW13 = jetP4s[1]+jetP4s[3];
    const auto hadT = hadW12+jetP4s[3];
    b_hadW12_pt = hadW12.Pt(); b_hadW12_eta = hadW12.Eta(); b_hadW12_dphi = hadW12.Phi(); b_hadW12_m = hadW12.M();
    b_hadW23_pt = hadW23.Pt(); b_hadW23_eta = hadW23.Eta(); b_hadW23_dphi = hadW23.Phi(); b_hadW23_m = hadW23.M();
    b_hadW13_pt = hadW13.Pt(); b_hadW13_eta = hadW13.Eta(); b_hadW13_dphi = hadW13.Phi(); b_hadW13_m = hadW13.M();
    b_hadW12_dR = jetP4s[1].DeltaR(jetP4s[2]);
    b_hadW23_dR = jetP4s[2].DeltaR(jetP4s[3]);
    b_hadW13_dR = jetP4s[1].DeltaR(jetP4s[3]);
    b_hadT_pt = hadT.Pt(); b_hadT_eta = hadT.Eta(); b_hadT_dphi = hadT.Phi(); b_hadT_m = hadT.M();

    TLorentzVector cm_hb = jetP4s[3], cm_hj1 = jetP4s[1], cm_hj2 = jetP4s[2];
    TLorentzVector cm_top = cm_hb+cm_hj1+cm_hj2;
    cm_hb.Boost(-cm_top.BoostVector());
    cm_hj1.Boost(-cm_top.BoostVector());
    cm_hj2.Boost(-cm_top.BoostVector());
    cm_hb *= 1./cm_hb.P();
    cm_hj1 *= 1./cm_hj1.P();
    cm_hj2 *= 1./cm_hj2.P();
    cm_top *= 1./cm_top.P();
    b_theta1 = cm_hb.Vect().Dot(cm_top.Vect());
    b_theta2 = cm_hb.Vect().Cross(cm_hj1.Vect()).Dot(cm_hb.Vect().Cross(cm_top.Vect()));

    b_lepB_CSV = jets_bTag[bestIdxs[0]];
    b_hadJ1_CSV = jets_bTag[bestIdxs[1]];
    b_hadJ2_CSV = jets_bTag[bestIdxs[2]];
    b_hadB_CSV = jets_bTag[bestIdxs[3]];

    b_lepB_CvsB = 0;
    b_hadJ1_CvsB = 0;
    b_hadJ2_CvsB = 0;
    b_hadB_CvsB = 0;

    b_lepB_CvsL = 0;
    b_hadJ1_CvsL = 0;
    b_hadJ2_CvsL = 0;
    b_hadB_CvsL = 0;

    // Fill the jet image
    const double eta1 = jets_eta[bestIdxs[1]]-b_hadT_eta, phi1 = deltaPhi(jets_phi[bestIdxs[1]], b_hadT_dphi);
    const double eta2 = jets_eta[bestIdxs[2]]-b_hadT_eta, phi2 = deltaPhi(jets_phi[bestIdxs[2]], b_hadT_dphi);
    const double eta3 = jets_eta[bestIdxs[3]]-b_hadT_eta, phi3 = deltaPhi(jets_phi[bestIdxs[3]], b_hadT_dphi);
    const double theta3 = atan2(phi3, eta3); // rotational symmetry. rotate trijet along the 3rd jet
    const double y1 = -eta1*sin(theta3) + phi1*cos(theta3);
    const double y2 = -eta2*sin(theta3) + phi2*cos(theta3);
    const int flipSign = (y2 > y1 ? -1 : 1);
    for ( int i=0; i<subjets_n; ++i ) {
      const size_t jetIdx = subjets_jetIdx[i];
      if ( jetIdx != bestIdxs[1] and jetIdx != bestIdxs[2] and jetIdx != bestIdxs[3] ) continue;

      const double eta = subjets_eta[i]-b_hadT_eta; // translate particle to the trijet centre
      const double phi = subjets_phi[i]-b_hadT_dphi; // translate particle to the trijet centre
      const double x =   eta*cos(theta3) + phi*sin(theta3); // Rotate particle by the 3rd jet orientation
      const double y = (-eta*sin(theta3) + phi*cos(theta3))*flipSign; // Rotation ,but also the mirror symmetry

      const double pt = subjets_pt[i];
      const int q = subjets_q[i];
      const int pid = subjets_pdgId[i];

      if ( q != 0 ) {
        hJetImage_ch_n->Fill(x, y);
        hJetImage_ch_pt->Fill(x, y, pt);
        b_hJetImage_ch_n->Fill(x, y);
        b_hJetImage_ch_pt->Fill(x, y, pt);
      }
      else if ( pid == 22 ) {
        hJetImage_ph_n->Fill(x, y);
        hJetImage_ph_pt->Fill(x, y, pt);
        b_hJetImage_ph_n->Fill(x, y);
        b_hJetImage_ph_pt->Fill(x, y, pt);
      }
      else {
        hJetImage_nh_n->Fill(x, y);
        hJetImage_nh_pt->Fill(x, y, pt);
        b_hJetImage_nh_n->Fill(x, y);
        b_hJetImage_nh_pt->Fill(x, y, pt);
      }
    }

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
      b_addJetByPt1_pt = b_addJetByPt2_pt = b_addJetByCSV1_pt = b_addJetByCSV2_pt = 0;
      b_addJetByPt1_CSV = b_addJetByPt2_CSV = b_addJetByCSV1_CSV = b_addJetByCSV2_CSV = -10;
      b_addJetsByPt_m = b_addJetsByCSV_m = 0;
      b_addJetsByPt_dR = b_addJetsByCSV_dR = 0;
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

      b_addJetByPt1_pt = jets_pt[jByPt1];
      b_addJetByPt2_pt = jets_pt[jByPt2];
      b_addJetByCSV1_pt = jets_pt[jByCSV1];
      b_addJetByCSV2_pt = jets_pt[jByCSV2];
      b_addJetByPt1_CSV = jets_bTag[jByPt1];
      b_addJetByPt2_CSV = jets_bTag[jByPt2];
      b_addJetByCSV1_CSV = jets_bTag[jByCSV1];
      b_addJetByCSV2_CSV = jets_bTag[jByCSV2];
      b_addJetsByPt_dR = addJetByPt1.DeltaR(addJetByPt2);
      b_addJetsByCSV_dR = addJetByCSV1.DeltaR(addJetByCSV2);
      b_addJetsByPt_m = (addJetByPt1+addJetByPt2).M();
      b_addJetsByCSV_m = (addJetByCSV1+addJetByCSV2).M();

      hAddJJ_dR->Fill(b_addJetsByCSV_dR);
      hAddJJ_m->Fill(b_addJetsByCSV_m);

    }

    // Rotate by lepton phi
    rotate(b_met_dphi, b_lepton_phi);
    rotate(b_lepB_dphi, b_lepton_phi);
    rotate(b_lepW_dphi, b_lepton_phi);
    rotate(b_lepT_dphi, b_lepton_phi);
    rotate(b_hadJ1_dphi, b_lepton_phi);
    rotate(b_hadJ2_dphi, b_lepton_phi);
    rotate(b_hadB_dphi, b_lepton_phi);
    rotate(b_hadW12_dphi, b_lepton_phi);
    rotate(b_hadW23_dphi, b_lepton_phi);
    rotate(b_hadW13_dphi, b_lepton_phi);
    rotate(b_hadT_dphi, b_lepton_phi);

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
  c = new TCanvas("c_hJetImage_ch_n", "hJetImage_ch_n", 500, 500); hJetImage_ch_n->Draw("COLZ");
  c = new TCanvas("c_hJetImage_nh_n", "hJetImage_nh_n", 500, 500); hJetImage_nh_n->Draw("COLZ");
  c = new TCanvas("c_hJetImage_ph_n", "hJetImage_ph_n", 500, 500); hJetImage_ph_n->Draw("COLZ");
  c = new TCanvas("c_hJetImage_ch_pt", "hJetImage_ch_pt", 500, 500); hJetImage_ch_pt->Draw("COLZ");
  c = new TCanvas("c_hJetImage_nh_pt", "hJetImage_nh_pt", 500, 500); hJetImage_nh_pt->Draw("COLZ");
  c = new TCanvas("c_hJetImage_ph_pt", "hJetImage_ph_pt", 500, 500); hJetImage_ph_pt->Draw("COLZ");
}
