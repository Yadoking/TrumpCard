#define AnalyzerHYTuple_cxx

#include "TTLJKinFit.h"

#include "AnalyzerHYTuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>

void AnalyzerHYTuple::Loop(const string modeStr, const string outFileName)
{
  //   In a ROOT session, you can do:
  //      root> .L AnalyzerHYTuple.C
  //      root> AnalyzerHYTuple t
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

  TFile* fout = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");
  tree->SetDirectory(fout);

  TH1F* hLW_m = new TH1F("hLW_m", "Leptonic W Mass", 100, 0, 300);
  TH1F* hLT_m = new TH1F("hLT_m", "Leptonic Top Mass", 100, 0, 300);
  TH1F* hHW_m = new TH1F("hHW_m", "Hadronic W Mass", 100, 0, 300);
  TH1F* hHW_dR = new TH1F("hHW_dR", "Hadronic W DeltaR", 100, 0, 5);
  TH1F* hHT_m = new TH1F("hHT_m", "Hadronic Top Mass", 100, 0, 300);
  TH1F* hJES = new TH1F("hJES", "Residual JES", 100, 0, 2);
  TH1F* heta = new TH1F("heta", "eta", 100, -5, 5);
  TH1F* hAddJJ_dR = new TH1F("hAddJJ_dR", "hAddJJ_dR", 100, 0, 5);
  TH1F* hAddJJ_m = new TH1F("hAddJJ_m", "hAddJJ_m", 100, 0, 300);

  int b_run, b_event, b_vertex_n;
  int b_jets_n, b_bjets_n;
  float b_weight_gen;
  float b_lepton_pt, b_lepton_eta, b_lepton_phi;
  float b_met_pt, b_met_phi;
  float b_kin_chi2;
  int b_kin_bjetcode;
  float b_kin_lep_pt, b_kin_lep_eta, b_kin_lep_phi;
  float b_kin_nu_pt, b_kin_nu_eta, b_kin_nu_phi;
  float b_kin_lepB_pt, b_kin_lepB_eta, b_kin_lepB_phi, b_kin_lepB_m;
  float b_kin_lepW_pt, b_kin_lepW_eta, b_kin_lepW_phi, b_kin_lepW_m;
  float b_kin_lepT_pt, b_kin_lepT_eta, b_kin_lepT_phi, b_kin_lepT_m;
  float b_kin_hadJ1_pt, b_kin_hadJ1_eta, b_kin_hadJ1_phi, b_kin_hadJ1_m;
  float b_kin_hadJ2_pt, b_kin_hadJ2_eta, b_kin_hadJ2_phi, b_kin_hadJ2_m;
  float b_kin_hadB_pt, b_kin_hadB_eta, b_kin_hadB_phi, b_kin_hadB_m;
  float b_kin_hadW12_pt, b_kin_hadW12_eta, b_kin_hadW12_phi, b_kin_hadW12_m;
  float b_kin_hadW23_pt, b_kin_hadW23_eta, b_kin_hadW23_phi, b_kin_hadW23_m;
  float b_kin_hadT_pt, b_kin_hadT_eta, b_kin_hadT_phi, b_kin_hadT_m;
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
  tree->Branch("kin_hadT_m", &b_kin_hadT_m, "kin_hadT_m/F");

  tree->Branch("kin_chi2", &b_kin_chi2, "kin_chi2/F");
  tree->Branch("kin_bjetcode", &b_kin_bjetcode, "kin_bjetcode/I"); // b-jet contribution "code". Format=[nbjetInLepT][nbjetInHadT]
  tree->Branch("kin_lep_pt", &b_kin_lep_pt, "kin_lep_pt/F");
  tree->Branch("kin_lep_eta", &b_kin_lep_eta, "kin_lep_eta/F");
  tree->Branch("kin_lep_phi", &b_kin_lep_phi, "kin_lep_phi/F");
  tree->Branch("kin_nu_pt", &b_kin_nu_pt, "kin_nu_pt/F");
  tree->Branch("kin_nu_eta", &b_kin_nu_eta, "kin_nu_eta/F");
  tree->Branch("kin_nu_phi", &b_kin_nu_phi, "kin_nu_phi/F");
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
  tree->Branch("kin_hadW23_pt", &b_kin_hadW23_pt, "kin_hadW23_pt/F");
  tree->Branch("kin_hadW23_eta", &b_kin_hadW23_eta, "kin_hadW23_eta/F");
  tree->Branch("kin_hadW23_phi", &b_kin_hadW23_phi, "kin_hadW23_phi/F");
  tree->Branch("kin_hadT_pt", &b_kin_hadT_pt, "kin_hadT_pt/F");
  tree->Branch("kin_hadT_eta", &b_kin_hadT_eta, "kin_hadT_eta/F");
  tree->Branch("kin_hadT_phi", &b_kin_hadT_phi, "kin_hadT_phi/F");

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

  TTLJKinFit fit(-1); // do not apply mass constraint on the dijet among trijets system

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

    double bestChi2 = 1e9;
    std::vector<size_t> bestIdxs;
    TLorentzVector jetP4s[4];
    for ( size_t j1 : jetIdxs ) {
      jetP4s[0].SetPtEtaPhiE(jets_pt[j1], jets_eta[j1], jets_phi[j1], jets_e[j1]);
      const int nbjetsInLepT = (jets_CSV[j1] > CSVM) ? 1 : 0;
      for ( size_t j2 : jetIdxs ) {
        if ( j2 == j1 ) continue;
        jetP4s[1].SetPtEtaPhiE(jets_pt[j2], jets_eta[j2], jets_phi[j2], jets_e[j2]);
        for ( size_t j3 : jetIdxs ) {
          if ( j2 >= j3 ) continue;
          if ( j3 == j1 ) continue;
          jetP4s[2].SetPtEtaPhiE(jets_pt[j3], jets_eta[j3], jets_phi[j3], jets_e[j3]);
          for ( size_t j4 : jetIdxs ) {
            if ( j4 == j3 or j4 == j2 or j4 == j1 ) continue;
            // Apply b-tag requirements
            int nbjetsInHadT = 0;
            if ( jets_CSV[j2] > CSVM ) ++nbjetsInHadT;
            if ( jets_CSV[j3] > CSVM ) ++nbjetsInHadT;
            if ( jets_CSV[j4] > CSVM ) ++nbjetsInHadT;
            if ( mode == Mode::TT ) {
              // SM ttbar mode: require at least one b jet used
              if ( b_bjets_n >= 2 and nbjetsInHadT < 1 ) continue;
              else if ( b_bjets_n == 1 and nbjetsInHadT+nbjetsInLepT < 1 ) continue;
            }
            else if ( mode == Mode::FCNC ) {
              // FCNC mode: require b jets in the trijet system, t->Hc / H->bb
              // can be used for the charged Higgs case as well, t->H+b / H+ -> cb
              if ( b_bjets_n >= 3 and nbjetsInHadT < 2 ) continue; // at least two b jets in hadronic side
              else if ( b_bjets_n == 2 and nbjetsInHadT < 1 ) continue; // at least one b jet in hadronic side
            }

            jetP4s[3].SetPtEtaPhiE(jets_pt[j4], jets_eta[j4], jets_phi[j4], jets_e[j4]);
            const double chi2 = fit.compute(metP4, leptonP4, jetP4s[0], jetP4s[1], jetP4s[2], jetP4s[3]);

            if ( chi2 < bestChi2 ) {
              bestChi2 = chi2;
              bestIdxs = {j1, j2, j3, j4};
              b_kin_bjetcode = nbjetsInLepT*10 + nbjetsInHadT; // b jet contribution "code"
            }
          }
        }
      }
    }
    if ( bestIdxs.empty() ) continue;
    // Sort by CSV in increasing order - j1+j2 will prefer 80GeV for SM top, j2+j3 will prefer 125GeV for FCNC top
    // Keep j0 as is which is the jet from leptonically decaying top
    std::sort(std::next(bestIdxs.begin()), bestIdxs.end(),
              [&](size_t a, size_t b){ return jets_CSV[a] < jets_CSV[b]; });

    for ( size_t i=0; i<4; ++i ) {
      const size_t j = bestIdxs[i];
      jetP4s[i].SetPtEtaPhiE(jets_pt[j], jets_eta[j], jets_phi[j], jets_e[j]);
    }
    b_kin_chi2 = fit.compute(metP4, leptonP4, jetP4s[0], jetP4s[1], jetP4s[2], jetP4s[3]);
    const std::vector<TLorentzVector> solution = fit.getSolution();
    const auto& sol_nuP4 = solution[0], sol_lepP4 = solution[1], sol_ljP4 = solution[2];
    const auto& sol_wj1P4 = solution[3], sol_wj2P4 = solution[4], sol_hbP4 = solution[5];

    hLW_m->Fill( (sol_lepP4+sol_nuP4).M() );
    hLT_m->Fill( (sol_lepP4+sol_nuP4+sol_ljP4).M() );
    if ( mode == Mode::TT ) {
      hHW_m->Fill( (sol_wj1P4+sol_wj2P4).M() );
      hHW_dR->Fill( sol_wj1P4.DeltaR(sol_wj2P4) );
    }
    else if ( mode == Mode::FCNC ) {
      hHW_m->Fill( (sol_wj2P4+sol_hbP4).M() );
      hHW_dR->Fill( sol_wj2P4.DeltaR(sol_hbP4) );
    }
    hHT_m->Fill( (sol_wj1P4+sol_wj2P4+sol_hbP4).M() );

    hJES->Fill(fit.min_->X()[0]);
    heta->Fill(fit.min_->X()[1]);

    b_kin_lep_pt = sol_lepP4.Pt(); b_kin_lep_eta = sol_lepP4.Eta(); b_kin_lep_phi = sol_lepP4.Phi();
    b_kin_nu_pt = sol_nuP4.Pt(); b_kin_nu_eta = sol_nuP4.Eta(); b_kin_nu_phi = sol_nuP4.Phi(); 
    b_kin_lepB_pt = sol_ljP4.Pt(); b_kin_lepB_eta = sol_ljP4.Eta(); b_kin_lepB_phi = sol_ljP4.Phi(); b_kin_lepB_m = sol_ljP4.M(); 
    b_kin_hadJ1_pt = sol_wj1P4.Pt(); b_kin_hadJ1_eta = sol_wj1P4.Eta(); b_kin_hadJ1_phi = sol_wj1P4.Phi(); b_kin_hadJ1_m = sol_wj1P4.M(); 
    b_kin_hadJ2_pt = sol_wj2P4.Pt(); b_kin_hadJ2_eta = sol_wj2P4.Eta(); b_kin_hadJ2_phi = sol_wj2P4.Phi(); b_kin_hadJ2_m = sol_wj2P4.M(); 
    b_kin_hadB_pt = sol_hbP4.Pt(); b_kin_hadB_eta = sol_hbP4.Eta(); b_kin_hadB_phi = sol_hbP4.Phi(); b_kin_hadB_m = sol_hbP4.M(); 

    const auto lepW = sol_lepP4+sol_nuP4;
    const auto lepT = lepW+sol_ljP4;
    b_kin_lepW_pt = lepW.Pt(); b_kin_lepW_eta = lepW.Eta(); b_kin_lepW_phi = lepW.Phi(); b_kin_lepW_m = lepW.M(); 
    b_kin_lepT_pt = lepT.Pt(); b_kin_lepT_eta = lepT.Eta(); b_kin_lepT_phi = lepT.Phi(); b_kin_lepT_m = lepT.M(); 

    const auto hadW12 = sol_wj1P4+sol_wj2P4;
    const auto hadW23 = sol_wj2P4+sol_hbP4;
    const auto hadT = hadW12+sol_hbP4;
    b_kin_hadW12_pt = hadW12.Pt(); b_kin_hadW12_eta = hadW12.Eta(); b_kin_hadW12_phi = hadW12.Phi(); b_kin_hadW12_m = hadW12.M(); 
    b_kin_hadW23_pt = hadW23.Pt(); b_kin_hadW23_eta = hadW23.Eta(); b_kin_hadW23_phi = hadW23.Phi(); b_kin_hadW23_m = hadW23.M(); 
    b_kin_hadT_pt = hadT.Pt(); b_kin_hadT_eta = hadT.Eta(); b_kin_hadT_phi = hadT.Phi(); b_kin_hadT_m = hadT.M(); 

    b_kin_lepB_CSV = jets_CSV[bestIdxs[0]];
    b_kin_hadJ1_CSV = jets_CSV[bestIdxs[1]];
    b_kin_hadJ2_CSV = jets_CSV[bestIdxs[2]];
    b_kin_hadB_CSV = jets_CSV[bestIdxs[3]];

    b_kin_lepB_CvsB = jets_CvsB[bestIdxs[0]];
    b_kin_hadJ1_CvsB = jets_CvsB[bestIdxs[1]];
    b_kin_hadJ2_CvsB = jets_CvsB[bestIdxs[2]];
    b_kin_hadB_CvsB = jets_CvsB[bestIdxs[3]];

    b_kin_lepB_CvsL = jets_CvsL[bestIdxs[0]];
    b_kin_hadJ1_CvsL = jets_CvsL[bestIdxs[1]];
    b_kin_hadJ2_CvsL = jets_CvsL[bestIdxs[2]];
    b_kin_hadB_CvsL = jets_CvsL[bestIdxs[3]];

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
      std::sort(addJetIdxs.begin(), addJetIdxs.end(), [&](size_t a, size_t b) { return jets_CSV[a] > jets_CSV[b]; });
      const size_t jByPt1 = addJetIdxsByPt[0], jByPt2 = addJetIdxsByPt[1];
      const size_t jByCSV1 = addJetIdxs[0], jByCSV2 = addJetIdxs[1];
      TLorentzVector addJetByCSV1, addJetByCSV2, addJetByPt1, addJetByPt2;
      addJetByCSV1.SetPtEtaPhiE(jets_pt[jByPt1], jets_eta[jByPt1], jets_phi[jByPt1], jets_e[jByPt1]);
      addJetByCSV1.SetPtEtaPhiE(jets_pt[jByPt2], jets_eta[jByPt2], jets_phi[jByPt2], jets_e[jByPt2]);
      addJetByPt1.SetPtEtaPhiE(jets_pt[jByCSV1], jets_eta[jByCSV1], jets_phi[jByCSV1], jets_e[jByCSV1]);
      addJetByPt1.SetPtEtaPhiE(jets_pt[jByCSV2], jets_eta[jByCSV2], jets_phi[jByCSV2], jets_e[jByCSV2]);

      b_kin_addJetByPt1_pt = jets_pt[jByPt1];
      b_kin_addJetByPt2_pt = jets_pt[jByPt2];
      b_kin_addJetByCSV1_pt = jets_pt[jByCSV1];
      b_kin_addJetByCSV2_pt = jets_pt[jByCSV2];
      b_kin_addJetByPt1_CSV = jets_CSV[jByPt1];
      b_kin_addJetByPt2_CSV = jets_CSV[jByPt2];
      b_kin_addJetByCSV1_CSV = jets_CSV[jByCSV1];
      b_kin_addJetByCSV2_CSV = jets_CSV[jByCSV2];
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
  c = new TCanvas("cJES", "JES", 500, 500); hJES->Draw();
  c = new TCanvas("ceta", "eta", 500, 500); heta->Draw();
}
