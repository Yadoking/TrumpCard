#define AnalyzeBasic_cxx

#include "AnalyzeBasic.h"
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>

const double CSVM = 0.5; // this is just to distinguish 0 or 1 for Delphes

AnalyzeBasic::TTLJSolution AnalyzeBasic::solveByDeltaR(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const
{
  // Find solution with minimum deltaR
  if ( jetIdxs.size() < 4 ) return TTLJSolution();

  // Sort by pT
  std::sort(jetIdxs.begin(), jetIdxs.end(),
            [&](size_t a, size_t b) {return jets_pt[a] > jets_pt[b];});
  // Then move b-jets to front, keeping the pT ordering
  std::stable_sort(jetIdxs.begin(), jetIdxs.end(),
                   [&](size_t a, size_t b) {return (jets_bTag[a]>CSVM) > (jets_bTag[b]>CSVM);});

  TTLJSolution sol;
  // Find a dijet with minDR
  double minDR = 1e9;
  TLorentzVector j1P4, j2P4;
  for ( auto i = jetIdxs.begin(); i != jetIdxs.end(); ++i ) {
    const size_t ii = *i;
    if ( jets_bTag[ii] < CSVM ) continue;
    j1P4.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    for ( auto j = std::next(i); j != jetIdxs.end(); ++j ) { 
      const size_t jj = *j;
      if ( jets_bTag[jj] < CSVM ) continue;
      j2P4.SetPtEtaPhiM(jets_pt[jj], jets_eta[jj], jets_phi[jj], jets_m[jj]);
      const double dR = j1P4.DeltaR(j2P4);
      if ( dR < minDR ) {
        minDR = dR;
        sol.hadJ1 = j1P4;
        sol.hadJ2 = j2P4;
        sol.hadJ1Idx = ii;
        sol.hadJ2Idx = jj;
      }
    }
  }
  if ( minDR >= 1e8 ) return TTLJSolution(); // Can happen if there's no b-jet pair
  sol.quality = minDR;

  // Continue to match trijet with minDR
  minDR = 1e9;
  TLorentzVector wP4 = j1P4+j2P4, j3P4;
  for ( auto ii : jetIdxs ) {
    if ( ii == sol.hadJ1Idx or ii == sol.hadJ2Idx ) continue;
    j3P4.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    const double dR = wP4.DeltaR(j3P4);
    if ( dR < minDR ) {
      minDR = dR;
      sol.hadJ3 = j3P4;
      sol.hadJ3Idx = ii;
    }
  }
  if ( minDR >= 1e8 ) return TTLJSolution();

  // Then move to the leptonic part, just pick the first non-overalpping one
  for ( auto ii : jetIdxs ) {
    if ( ii == sol.hadJ1Idx or ii == sol.hadJ2Idx or ii == sol.hadJ3Idx ) continue;
    sol.lepB.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    sol.lepBIdx = ii;
    break;
  }
  sol.lep = lepP4;
  sol.met = metP4;

  sol.isValid = true;
  return sol;
}

AnalyzeBasic::TTLJSolution AnalyzeBasic::solveByM3(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const
{
  // Find solution by M3 - trijet with maximum pT
  if ( jetIdxs.size() < 4 ) return TTLJSolution();

  // Sort by pT
  std::sort(jetIdxs.begin(), jetIdxs.end(),
            [&](size_t a, size_t b) {return jets_pt[a] > jets_pt[b];});
  // Then move b-jets to front, keeping the pT ordering
  std::stable_sort(jetIdxs.begin(), jetIdxs.end(),
                   [&](size_t a, size_t b) {return (jets_bTag[a]>CSVM) > (jets_bTag[b]>CSVM);});

  TTLJSolution sol;
  // Find a dijet with minDR
  double maxPt = 0;
  TLorentzVector j1P4, j2P4, j3P4;
  for ( auto i = jetIdxs.begin(); i != jetIdxs.end(); ++i ) {
    const size_t ii = *i;
    if ( jets_bTag[ii] < CSVM ) continue;
    j1P4.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);

    for ( auto j = std::next(i); j != jetIdxs.end(); ++j ) { 
      const size_t jj = *j;
      if ( jets_bTag[jj] < CSVM ) continue;
      j2P4.SetPtEtaPhiM(jets_pt[jj], jets_eta[jj], jets_phi[jj], jets_m[jj]);

      for ( auto kk : jetIdxs ) {
        if ( ii == kk or jj == kk ) continue;
        j3P4.SetPtEtaPhiM(jets_pt[kk], jets_eta[kk], jets_phi[kk], jets_m[kk]);

        const double pt = (j1P4+j2P4+j3P4).Pt();
        if ( pt > maxPt ) {
          maxPt = pt;
          sol.hadJ1 = j1P4;
          sol.hadJ2 = j2P4;
          sol.hadJ3 = j3P4;
          sol.hadJ1Idx = ii;
          sol.hadJ2Idx = jj;
          sol.hadJ3Idx = kk;
        }
      }
    }
  }
  if ( maxPt <= 0 ) return TTLJSolution(); // Can happen if there's no b-jet pair
  sol.quality = maxPt;

  // Then move to the leptonic part, just pick the first non-overalpping one
  for ( auto ii : jetIdxs ) {
    if ( ii == sol.hadJ1Idx or ii == sol.hadJ2Idx or ii == sol.hadJ3Idx ) continue;
    sol.lepB.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    sol.lepBIdx = ii;
    break;
  }
  sol.lep = lepP4;
  sol.met = metP4;

  sol.isValid = true;
  return sol;
}

AnalyzeBasic::TTLJSolution AnalyzeBasic::solveByMTop(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const
{
  // Find solution by closest mass to the PDG average
  if ( jetIdxs.size() < 4 ) return TTLJSolution();

  // Sort by pT
  std::sort(jetIdxs.begin(), jetIdxs.end(),
            [&](size_t a, size_t b) {return jets_pt[a] > jets_pt[b];});
  // Then move b-jets to front, keeping the pT ordering
  std::stable_sort(jetIdxs.begin(), jetIdxs.end(),
                   [&](size_t a, size_t b) {return (jets_bTag[a]>CSVM) > (jets_bTag[b]>CSVM);});

  TTLJSolution sol;

  double minChi2 = 1e9;
  TLorentzVector j1P4, j2P4, j3P4;
  for ( auto i = jetIdxs.begin(); i != jetIdxs.end(); ++i ) {
    const size_t ii = *i;
    if ( jets_bTag[ii] < CSVM ) continue;
    j1P4.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    for ( auto j = std::next(i); j != jetIdxs.end(); ++j ) {
      const size_t jj = *j;
      if ( jets_bTag[jj] < CSVM ) continue;
      j2P4.SetPtEtaPhiM(jets_pt[jj], jets_eta[jj], jets_phi[jj], jets_m[jj]);
      const auto wP4 = j1P4+j2P4;
      const double dmW = 0; // Just a dummy
      //const double dmW = std::abs(wP4.M()-125);

      for ( auto k = jetIdxs.begin(); k != jetIdxs.end(); ++k ) {
        const size_t kk = *k;
        if ( kk == ii or kk == jj ) continue;
        j3P4.SetPtEtaPhiM(jets_pt[kk], jets_eta[kk], jets_phi[kk], jets_m[kk]);
        const auto tP4 = wP4+j3P4;
        const double dmT = std::abs(tP4.M()-172.5);
        if ( dmW+dmT < minChi2 ) {
          minChi2 = dmW+dmT;
          sol.hadJ1 = j1P4;
          sol.hadJ2 = j2P4;
          sol.hadJ3 = j3P4;
          sol.hadJ1Idx = ii;
          sol.hadJ2Idx = jj;
          sol.hadJ3Idx = kk;
        }
      }
    }
  }

  if ( minChi2 >= 1e8 ) return TTLJSolution(); // Can happen if there's no b-jet pair
  sol.quality = minChi2;

  // Then move to the leptonic part, just pick the first non-overalpping one
  for ( auto ii : jetIdxs ) {
    if ( ii == sol.hadJ1Idx or ii == sol.hadJ2Idx or ii == sol.hadJ3Idx ) continue;
    sol.lepB.SetPtEtaPhiM(jets_pt[ii], jets_eta[ii], jets_phi[ii], jets_m[ii]);
    sol.lepBIdx = ii;
    break;
  }
  sol.lep = lepP4;
  sol.met = metP4;

  sol.isValid = true;
  return sol;
}

void AnalyzeBasic::Loop(const string outFileName)
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

  TFile* fout = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");
  tree->SetDirectory(fout);

  // Book output tree branches
  int b_run, b_event;
  int b_vertex_n;

  float b_weight_gen;
  float b_genLepT_pt, b_genLepT_eta, b_genLepT_phi, b_genLepT_m;
  float b_genHadT_pt, b_genHadT_eta, b_genHadT_phi, b_genHadT_m;
  float b_genHadJJ_pt, b_genHadJJ_eta, b_genHadJJ_phi, b_genHadJJ_m, b_genHadJJ_dR;
  int b_genMatch;

  int b_jets_n, b_bjets_n;
  
  float b_lepton_pt, b_lepton_eta, b_lepton_phi;
  float b_met_pt, b_met_phi;

  //int b_bjetcode;
  float b_lepB_pt, b_lepB_eta, b_lepB_phi, b_lepB_m;
  float b_lepW_pt, b_lepW_eta, b_lepW_phi, b_lepW_m;
  float b_lepT_pt, b_lepT_eta, b_lepT_phi, b_lepT_m;
  float b_hadJ1_pt, b_hadJ1_eta, b_hadJ1_phi, b_hadJ1_m;
  float b_hadJ2_pt, b_hadJ2_eta, b_hadJ2_phi, b_hadJ2_m;
  float b_hadJ3_pt, b_hadJ3_eta, b_hadJ3_phi, b_hadJ3_m;
  float b_hadW12_pt, b_hadW12_eta, b_hadW12_phi, b_hadW12_m, b_hadW12_dR;
  float b_hadW23_pt, b_hadW23_eta, b_hadW23_phi, b_hadW23_m, b_hadW23_dR;
  float b_hadW13_pt, b_hadW13_eta, b_hadW13_phi, b_hadW13_m, b_hadW13_dR;
  float b_hadT_pt, b_hadT_eta, b_hadT_phi, b_hadT_m;
  float b_theta1, b_theta2;
  float b_lepB_bTag, b_hadJ1_bTag, b_hadJ2_bTag, b_hadJ3_bTag;

  tree->Branch("run", &b_run, "run/I");
  tree->Branch("event", &b_event, "event/I");
  tree->Branch("vertex_n", &b_vertex_n, "vertex_n/I");

  tree->Branch("weight_gen", &b_weight_gen, "weight_gen/F");

  tree->Branch("genLepT_pt" , &b_genLepT_pt , "genLepT_pt/F" );
  tree->Branch("genLepT_eta", &b_genLepT_eta, "genLepT_eta/F");
  tree->Branch("genLepT_phi", &b_genLepT_phi, "genLepT_phi/F");
  tree->Branch("genLepT_m"  , &b_genLepT_m  , "genLepT_m/F"  );

  tree->Branch("genHadT_pt" , &b_genHadT_pt , "genHadT_pt/F" );
  tree->Branch("genHadT_eta", &b_genHadT_eta, "genHadT_eta/F");
  tree->Branch("genHadT_phi", &b_genHadT_phi, "genHadT_phi/F");
  tree->Branch("genHadT_m"  , &b_genHadT_m  , "genHadT_m/F"  );

  tree->Branch("genHadJJ_pt" , &b_genHadJJ_pt , "genHadJJ_pt/F" );
  tree->Branch("genHadJJ_eta", &b_genHadJJ_eta, "genHadJJ_eta/F");
  tree->Branch("genHadJJ_phi", &b_genHadJJ_phi, "genHadJJ_phi/F");
  tree->Branch("genHadJJ_m"  , &b_genHadJJ_m  , "genHadJJ_m/F"  );
  tree->Branch("genHadJJ_dR" , &b_genHadJJ_dR , "genHadJJ_dR/F" );

  tree->Branch("lepton_pt", &b_lepton_pt, "lepton_pt/F");
  tree->Branch("lepton_eta", &b_lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &b_lepton_phi, "lepton_phi/F");
  tree->Branch("met_pt", &b_met_pt, "met_pt/F");
  tree->Branch("met_phi", &b_met_phi, "met_phi/F");

  tree->Branch("jets_n", &b_jets_n, "jets_n/I");
  tree->Branch("bjets_n", &b_bjets_n, "bjets_n/I");

  tree->Branch("lepW_m", &b_lepW_m, "lepW_m/F");
  tree->Branch("lepT_m", &b_lepT_m, "lepT_m/F");
  tree->Branch("hadW12_m", &b_hadW12_m, "hadW12_m/F");
  tree->Branch("hadW23_m", &b_hadW23_m, "hadW23_m/F");
  tree->Branch("hadW13_m", &b_hadW13_m, "hadW13_m/F");
  tree->Branch("hadT_m", &b_hadT_m, "hadT_m/F");

  //tree->Branch("bjetcode", &b_bjetcode, "bjetcode/I"); // b-jet contribution "code". Format=[nbjetInLepT][nbjetInHadT]
  tree->Branch("lepB_pt", &b_lepB_pt, "lepB_pt/F");
  tree->Branch("lepB_eta", &b_lepB_eta, "lepB_eta/F");
  tree->Branch("lepB_phi", &b_lepB_phi, "lepB_phi/F");
  tree->Branch("lepB_m", &b_lepB_m, "lepB_m/F");
  tree->Branch("lepW_pt", &b_lepW_pt, "lepW_pt/F");
  tree->Branch("lepW_eta", &b_lepW_eta, "lepW_eta/F");
  tree->Branch("lepW_phi", &b_lepW_phi, "lepW_phi/F");
  tree->Branch("lepT_pt", &b_lepT_pt, "lepT_pt/F");
  tree->Branch("lepT_eta", &b_lepT_eta, "lepT_eta/F");
  tree->Branch("lepT_phi", &b_lepT_phi, "lepT_phi/F");
  tree->Branch("hadJ1_pt", &b_hadJ1_pt, "hadJ1_pt/F");
  tree->Branch("hadJ1_eta", &b_hadJ1_eta, "hadJ1_eta/F");
  tree->Branch("hadJ1_phi", &b_hadJ1_phi, "hadJ1_phi/F");
  tree->Branch("hadJ1_m", &b_hadJ1_m, "hadJ1_m/F");
  tree->Branch("hadJ2_pt", &b_hadJ2_pt, "hadJ2_pt/F");
  tree->Branch("hadJ2_eta", &b_hadJ2_eta, "hadJ2_eta/F");
  tree->Branch("hadJ2_phi", &b_hadJ2_phi, "hadJ2_phi/F");
  tree->Branch("hadJ2_m", &b_hadJ2_m, "hadJ2_m/F");
  tree->Branch("hadJ3_pt", &b_hadJ3_pt, "hadJ3_pt/F");
  tree->Branch("hadJ3_eta", &b_hadJ3_eta, "hadJ3_eta/F");
  tree->Branch("hadJ3_phi", &b_hadJ3_phi, "hadJ3_phi/F");
  tree->Branch("hadJ3_m", &b_hadJ3_m, "hadJ3_m/F");
  tree->Branch("hadW12_pt", &b_hadW12_pt, "hadW12_pt/F");
  tree->Branch("hadW12_eta", &b_hadW12_eta, "hadW12_eta/F");
  tree->Branch("hadW12_phi", &b_hadW12_phi, "hadW12_phi/F");
  tree->Branch("hadW12_dR", &b_hadW12_dR, "hadW12_dR/F");
  tree->Branch("hadW23_pt", &b_hadW23_pt, "hadW23_pt/F");
  tree->Branch("hadW23_eta", &b_hadW23_eta, "hadW23_eta/F");
  tree->Branch("hadW23_phi", &b_hadW23_phi, "hadW23_phi/F");
  tree->Branch("hadW23_dR", &b_hadW23_dR, "hadW23_dR/F");
  tree->Branch("hadW13_pt", &b_hadW13_pt, "hadW13_pt/F");
  tree->Branch("hadW13_eta", &b_hadW13_eta, "hadW13_eta/F");
  tree->Branch("hadW13_phi", &b_hadW13_phi, "hadW13_phi/F");
  tree->Branch("hadW13_dR", &b_hadW13_dR, "hadW13_dR/F");
  tree->Branch("hadT_pt", &b_hadT_pt, "hadT_pt/F");
  tree->Branch("hadT_eta", &b_hadT_eta, "hadT_eta/F");
  tree->Branch("hadT_phi", &b_hadT_phi, "hadT_phi/F");

  tree->Branch("theta1", &b_theta1, "theta1/F"); // Angle between top and b
  tree->Branch("theta2", &b_theta2, "theta2/F"); // Angle between t-b and w->jj plane
  tree->Branch("lepB_bTag", &b_lepB_bTag, "lepB_bTag/F");
  tree->Branch("hadJ1_bTag", &b_hadJ1_bTag, "hadJ1_bTag/F");
  tree->Branch("hadJ2_bTag", &b_hadJ2_bTag, "hadJ2_bTag/F");
  tree->Branch("hadJ3_bTag", &b_hadJ3_bTag, "hadJ3_bTag/F");

  tree->Branch("genMatch", &b_genMatch, "genMatch/I");

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
    b_weight_gen = weight;

    // Start from the generaror level objects
    // Construct decay tree
    b_genLepT_pt = b_genLepT_eta = b_genLepT_phi = b_genLepT_m = 0;
    b_genHadT_pt = b_genHadT_eta = b_genHadT_phi = b_genHadT_m = 0;
    b_genHadJJ_pt = b_genHadJJ_eta = b_genHadJJ_phi = b_genHadJJ_m = b_genHadJJ_dR = 0;
    int gen_lepId = 0, gen_nuId = 0, gen_hadJ1Id = 0, gen_hadJ2Id = 0, gen_hadJ3Id = 0;
    TLorentzVector gen_lepT, gen_hadT;
    TLorentzVector gen_lep, gen_nu, gen_lepB, gen_hadJ1, gen_hadJ2, gen_hadJ3;
    for ( size_t i=0; i<gen_n; ++i ) {
      const int absId = abs(gen_pdgId[i]);
      if ( absId != 24 and absId != 25 ) continue;

      const int mother = gen_mother[i];
      if ( mother < 0 or abs(gen_pdgId[mother]) != 6 ) continue; // Should be from top quark decay

      const int sibling1 = gen_dau1[mother], sibling2 = gen_dau2[mother];
      if ( sibling2-sibling1 != 1 or sibling2 < 0 or sibling1 < 0 ) continue; // two siblingings (including itself), t->Wb or t->Hc
      const int sibling = (sibling1 == int(i)) ? sibling2 : sibling1;

      int dau1 = gen_dau1[i], dau2 = gen_dau2[i];
      if ( dau2-dau1 != 1 or dau1 < 0 or dau2 < 0 ) continue; // should have two daughters only
      int dau1Aid = abs(gen_pdgId[dau1]), dau2Aid = abs(gen_pdgId[dau2]);
      if ( dau1Aid > dau2Aid ) {
        swap(dau1Aid, dau2Aid);
        swap(dau1, dau2);
      }

      if ( dau1Aid <= 5 and dau2Aid <= 5 ) { // W->jj or H->bb
        gen_hadT.SetPtEtaPhiM(gen_pt[mother], gen_eta[mother], gen_phi[mother], gen_m[mother]);
        gen_hadJ1.SetPtEtaPhiM(gen_pt[dau1], gen_eta[dau1], gen_phi[dau1], gen_m[dau1]);
        gen_hadJ2.SetPtEtaPhiM(gen_pt[dau2], gen_eta[dau2], gen_phi[dau2], gen_m[dau2]);
        gen_hadJ3.SetPtEtaPhiM(gen_pt[sibling], gen_eta[sibling], gen_phi[sibling], gen_m[sibling]);
        gen_hadJ1Id = gen_pdgId[dau1];
        gen_hadJ2Id = gen_pdgId[dau2];
        gen_hadJ3Id = gen_pdgId[sibling];
      }
      else if ( dau1Aid >= 11 and dau1Aid <= 15 ) {
        if ( dau1Aid == 15 ) {
          for ( int j=gen_dau1[dau1]; j<=gen_dau2[dau1]; ++j ) {
            if ( j < 0 ) continue;
            const int aid = abs(gen_pdgId[j]);
            if ( aid == 11 or aid == 13 ) { dau1 = j; break; }
          }
        }

        gen_lepT.SetPtEtaPhiM(gen_pt[mother], gen_eta[mother], gen_phi[mother], gen_m[mother]);
        gen_lep.SetPtEtaPhiM(gen_pt[dau1], gen_eta[dau1], gen_phi[dau1], gen_m[dau1]);
        gen_nu.SetPtEtaPhiM(gen_pt[dau2], gen_eta[dau2], gen_phi[dau2], gen_m[dau2]);
        gen_lepB.SetPtEtaPhiM(gen_pt[sibling], gen_eta[sibling], gen_phi[sibling], gen_m[sibling]);
        gen_lepId = gen_pdgId[dau1];
        gen_nuId = gen_pdgId[dau2];
      }
    }
    if ( gen_lepId != 0 and gen_hadJ3Id != 0 ) { // If we could build full decay tree
      b_genLepT_pt = gen_lepT.Pt(); b_genLepT_eta = gen_lepT.Eta(); b_genLepT_phi = gen_lepT.Phi(); b_genLepT_m = gen_lepT.M();
      b_genHadT_pt = gen_hadT.Pt(); b_genHadT_eta = gen_hadT.Eta(); b_genHadT_phi = gen_hadT.Phi(); b_genHadT_m = gen_hadT.M();
      TLorentzVector gen_hadJJ = gen_hadJ1+gen_hadJ2;
      b_genHadJJ_pt = gen_hadJJ.Pt(); b_genHadJJ_eta = gen_hadJJ.Eta(); b_genHadJJ_phi = gen_hadJJ.Phi(); b_genHadJJ_m = gen_hadJJ.M();
      b_genHadJJ_dR = gen_hadJ1.DeltaR(gen_hadJ2);
    }

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
    b_met_phi = met_phi;

    std::vector<size_t> jetIdxs;
    b_bjets_n = 0;
    for ( size_t j=0; j<jets_n; ++j ) {
      if ( jets_pt[j] < 30 or std::abs(jets_eta[j]) > 2.5 ) continue;
      TLorentzVector jetP4;
      jetP4.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
      if ( jetP4.DeltaR(leptonP4) < 0.3 ) continue;
      if ( jets_bTag[j] > CSVM ) ++b_bjets_n;
      jetIdxs.push_back(j);
    }
    if ( jetIdxs.size() < 4 ) continue;
    b_jets_n = jetIdxs.size();

    auto sol = solveByDeltaR(leptonP4, metP4, jetIdxs);
    //auto sol = solveByM3(leptonP4, metP4, jetIdxs);
    //auto sol = solveByMTop(leptonP4, metP4, jetIdxs);
    if ( !sol.isValid ) continue;

    // Do the deltaR matching to the reconstructed objects
    b_genMatch = 0; // [lep][nu][lepB][hadJ3][hadJ2][hadJ1]
    if ( gen_lep.Pt()   > 0 and gen_lep.DeltaR(leptonP4)    < 0.1 ) b_genMatch += 100000;
    if ( gen_nu.Pt()    > 0 and gen_nu.DeltaPhi(metP4)      < 0.1 ) b_genMatch += 10000;
    if ( gen_lepB.Pt()  > 0 and gen_lepB.DeltaR(sol.lepB)   < 0.4 ) b_genMatch += 1000;
    if ( gen_hadJ3.Pt()  > 0 and gen_hadJ3.DeltaR(sol.hadJ3)  < 0.4 ) b_genMatch += 100;
    if ( gen_hadJ2.Pt() > 0 and (gen_hadJ2.DeltaR(sol.hadJ1) < 0.4 or gen_hadJ2.DeltaR(sol.hadJ2) < 0.4) ) b_genMatch += 10;
    if ( gen_hadJ1.Pt() > 0 and (gen_hadJ1.DeltaR(sol.hadJ1) < 0.4 or gen_hadJ1.DeltaR(sol.hadJ2) < 0.4) ) b_genMatch += 1;

    b_lepB_pt = sol.lepB.Pt(); b_lepB_eta = sol.lepB.Eta(); b_lepB_phi = sol.lepB.Phi(); b_lepB_m = sol.lepB.M();
    b_hadJ1_pt = sol.hadJ1.Pt(); b_hadJ1_eta = sol.hadJ1.Eta(); b_hadJ1_phi = sol.hadJ1.Phi(); b_hadJ1_m = sol.hadJ1.M();
    b_hadJ2_pt = sol.hadJ2.Pt(); b_hadJ2_eta = sol.hadJ2.Eta(); b_hadJ2_phi = sol.hadJ2.Phi(); b_hadJ2_m = sol.hadJ2.M();
    b_hadJ3_pt = sol.hadJ3.Pt(); b_hadJ3_eta = sol.hadJ3.Eta(); b_hadJ3_phi = sol.hadJ3.Phi(); b_hadJ3_m = sol.hadJ3.M();

    const auto lepW = leptonP4+metP4;
    const auto lepT = lepW+sol.lepB;
    b_lepW_pt = lepW.Pt(); b_lepW_eta = lepW.Eta(); b_lepW_phi = lepW.Phi(); b_lepW_m = lepW.M();
    b_lepT_pt = lepT.Pt(); b_lepT_eta = lepT.Eta(); b_lepT_phi = lepT.Phi(); b_lepT_m = lepT.M();

    const auto hadW12 = sol.hadJ1+sol.hadJ2;
    const auto hadW23 = sol.hadJ2+sol.hadJ3;
    const auto hadW13 = sol.hadJ1+sol.hadJ3;
    const auto hadT = hadW12+sol.hadJ3;
    b_hadW12_pt = hadW12.Pt(); b_hadW12_eta = hadW12.Eta(); b_hadW12_phi = hadW12.Phi(); b_hadW12_m = hadW12.M();
    b_hadW23_pt = hadW23.Pt(); b_hadW23_eta = hadW23.Eta(); b_hadW23_phi = hadW23.Phi(); b_hadW23_m = hadW23.M();
    b_hadW13_pt = hadW13.Pt(); b_hadW13_eta = hadW13.Eta(); b_hadW13_phi = hadW13.Phi(); b_hadW13_m = hadW13.M();
    b_hadW12_dR = sol.hadJ1.DeltaR(sol.hadJ2);
    b_hadW23_dR = sol.hadJ2.DeltaR(sol.hadJ3);
    b_hadW13_dR = sol.hadJ1.DeltaR(sol.hadJ3);
    b_hadT_pt = hadT.Pt(); b_hadT_eta = hadT.Eta(); b_hadT_phi = hadT.Phi(); b_hadT_m = hadT.M();

    b_lepB_bTag = jets_bTag[sol.lepBIdx];
    b_hadJ1_bTag = jets_bTag[sol.hadJ1Idx];
    b_hadJ2_bTag = jets_bTag[sol.hadJ2Idx];
    b_hadJ3_bTag = jets_bTag[sol.hadJ3Idx];

    TLorentzVector cm_hb = sol.hadJ3, cm_hj1 = sol.hadJ1, cm_hj2 = sol.hadJ2;
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

/*
    // Fill the jet image
    const double eta1 = jets_eta[bestIdxs[1]]-b_hadT_eta, phi1 = deltaPhi(jets_phi[bestIdxs[1]], b_hadT_phi);
    const double eta2 = jets_eta[bestIdxs[2]]-b_hadT_eta, phi2 = deltaPhi(jets_phi[bestIdxs[2]], b_hadT_phi);
    const double eta3 = jets_eta[bestIdxs[3]]-b_hadT_eta, phi3 = deltaPhi(jets_phi[bestIdxs[3]], b_hadT_phi);
    const double theta3 = atan2(phi3, eta3); // rotational symmetry. rotate trijet along the 3rd jet
    const double y1 = -eta1*sin(theta3) + phi1*cos(theta3);
    const double y2 = -eta2*sin(theta3) + phi2*cos(theta3);
    const int flipSign = (y2 > y1 ? -1 : 1);
    for ( int i=0; i<subjets_n; ++i ) {
      const size_t jetIdx = subjets_jetIdx[i];
      if ( jetIdx != bestIdxs[1] and jetIdx != bestIdxs[2] and jetIdx != bestIdxs[3] ) continue;

      const double eta = subjets_eta[i]-b_hadT_eta; // translate particle to the trijet centre
      const double phi = subjets_phi[i]-b_hadT_phi; // translate particle to the trijet centre
      const double x =   eta*cos(theta3) + phi*sin(theta3); // Rotate particle by the 3rd jet orientation
      const double y = (-eta*sin(theta3) + phi*cos(theta3))*flipSign; // Rotation ,but also the mirror symmetry

      const double pt = subjets_pt[i];
      const int q = subjets_q[i];
      const int pid = subjets_pdgId[i];
    }
*/

    // Rotate by lepton phi
    //rotate(b_met_phi, b_lepton_phi);

    tree->Fill();
  }

  fout->Write();

}
