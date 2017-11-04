//#ifdef __CLING__
#ifdef __CINT__
R__LOAD_LIBRARY(libDelphes)
#endif
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include <iostream>

//------------------------------------------------------------------------------
int getLast(TClonesArray* branch, const int iGen)
{
  const GenParticle* p = (const GenParticle*)branch->At(iGen);
  if ( p->D1 == -1 or p->D2 == -1 ) return iGen;

  for ( int i=p->D1; i<=p->D2; ++i ) {
    const GenParticle* dau = (const GenParticle*)branch->At(i);
    if ( p->PID == dau->PID ) return getLast(branch, i);
  }
  return iGen;
}

void makeFlatTuple(const std::string finName, const std::string foutName)
{
  // Prepare output tree
  TFile* fout = TFile::Open(foutName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");

  unsigned short b_run = 1;
  unsigned int b_event = 0;
  float b_weight = 0;

  float b_met_pt, b_met_phi;

  const unsigned short muons_N = 100;
  unsigned short b_muons_n;
  float b_muons_pt[muons_N], b_muons_eta[muons_N], b_muons_phi[muons_N], b_muons_m[muons_N];
  short b_muons_q[muons_N];
  float b_muons_relIso[muons_N];

  const unsigned short electrons_N = 100;
  unsigned short b_electrons_n;
  float b_electrons_pt[electrons_N], b_electrons_eta[electrons_N], b_electrons_phi[electrons_N], b_electrons_m[electrons_N];
  short b_electrons_q[electrons_N];
  float b_electrons_relIso[electrons_N];

  const unsigned short jets_N = 100;
  unsigned short b_jets_n;
  float b_jets_pt[jets_N], b_jets_eta[jets_N], b_jets_phi[jets_N], b_jets_m[jets_N];
  short b_jets_flav[jets_N];
  float b_jets_bTag[jets_N];

  const unsigned short gen_N = 1000;
  unsigned short b_gen_n;
  float b_gen_pt[gen_N], b_gen_eta[gen_N], b_gen_phi[gen_N], b_gen_m[gen_N];
  short b_gen_pdgId[gen_N], b_gen_q3[gen_N];
  short b_gen_mother[gen_N], b_gen_dau1[gen_N], b_gen_dau2[gen_N];

  const unsigned short subjets_N = 10000;
  unsigned short b_subjets_n;
  float b_subjets_pt[subjets_N], b_subjets_eta[subjets_N], b_subjets_phi[subjets_N];
  short b_subjets_q[subjets_N], b_subjets_pdgId[subjets_N];
  unsigned short b_subjets_jetIdx[subjets_N];

  tree->Branch("run", &b_run, "run/s");
  tree->Branch("event", &b_event, "event/i");
  tree->Branch("weight", &b_weight, "weight/F");

  tree->Branch("met_pt", &b_met_pt, "met_pt/F");
  tree->Branch("met_phi", &b_met_phi, "met_phi/F");

  tree->Branch("muons_n", &b_muons_n, "muons_n/s");
  tree->Branch("muons_pt", b_muons_pt, "muons_pt[muons_n]/F");
  tree->Branch("muons_eta", b_muons_eta, "muons_eta[muons_n]/F");
  tree->Branch("muons_phi", b_muons_phi, "muons_phi[muons_n]/F");
  tree->Branch("muons_m", b_muons_m, "muons_m[muons_n]/F");
  tree->Branch("muons_q", b_muons_q, "muons_q[muons_n]/S");
  tree->Branch("muons_relIso", b_muons_relIso, "muons_relIso[muons_n]/F");

  tree->Branch("electrons_n", &b_electrons_n, "electrons_n/s");
  tree->Branch("electrons_pt", b_electrons_pt, "electrons_pt[electrons_n]/F");
  tree->Branch("electrons_eta", b_electrons_eta, "electrons_eta[electrons_n]/F");
  tree->Branch("electrons_phi", b_electrons_phi, "electrons_phi[electrons_n]/F");
  tree->Branch("electrons_m", b_electrons_m, "electrons_m[electrons_n]/F");
  tree->Branch("electrons_q", b_electrons_q, "electrons_q[electrons_n]/S");
  tree->Branch("electrons_relIso", b_electrons_relIso, "electrons_relIso[electrons_n]/F");

  tree->Branch("jets_n", &b_jets_n, "jets_n/s");
  tree->Branch("jets_pt", b_jets_pt, "jets_pt[jets_n]/F");
  tree->Branch("jets_eta", b_jets_eta, "jets_eta[jets_n]/F");
  tree->Branch("jets_phi", b_jets_phi, "jets_phi[jets_n]/F");
  tree->Branch("jets_m", b_jets_m, "jets_m[jets_n]/F");
  tree->Branch("jets_flav", b_jets_flav, "jets_flav[jets_n]/S");
  tree->Branch("jets_bTag", b_jets_bTag, "jets_bTag[jets_n]/F");

  tree->Branch("gen_n", &b_gen_n, "gen_n/s");
  tree->Branch("gen_pt", b_gen_pt, "gen_pt[gen_n]/F");
  tree->Branch("gen_eta", b_gen_eta, "gen_eta[gen_n]/F");
  tree->Branch("gen_phi", b_gen_phi, "gen_phi[gen_n]/F");
  tree->Branch("gen_m", b_gen_m, "gen_m[gen_n]/F");
  tree->Branch("gen_pdgId", b_gen_pdgId, "gen_pdgId[gen_n]/S");
  tree->Branch("gen_q3", b_gen_q3, "gen_q3[gen_n]/S");
  tree->Branch("gen_mother", b_gen_mother, "gen_mother[gen_n]/S");
  tree->Branch("gen_dau1", b_gen_dau1, "gen_dau1[gen_n]/S");
  tree->Branch("gen_dau2", b_gen_dau2, "gen_dau2[gen_n]/S");

  tree->Branch("subjets_n", &b_subjets_n, "subjets_n/s");
  tree->Branch("subjets_pt", b_subjets_pt, "subjets_pt[subjets_n]/F");
  tree->Branch("subjets_eta", b_subjets_eta, "subjets_eta[subjets_n]/F");
  tree->Branch("subjets_phi", b_subjets_phi, "subjets_phi[subjets_n]/F");
  tree->Branch("subjets_q", b_subjets_q, "subjets_q[subjets_n]/S");
  tree->Branch("subjets_pdgId", b_subjets_pdgId, "subjets_pdgId[subjets_n]/S");
  tree->Branch("subjets_jetIdx", b_subjets_jetIdx, "subjets_jetIdx[subjets_n]/S");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(finName.c_str());

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout << entry << '/' << numberOfEntries << '\r';

    const HepMCEvent* event = (const HepMCEvent*)branchEvent->At(0);
    b_event = event->Number;
    b_weight = event->Weight;

    // Build gen particle collection, for the ttbar decays
    std::vector<std::vector<int> > gen_topDaus;
    b_gen_n = 0;
    for ( int i=0; i<branchGen->GetEntries(); ++i ) {
      const GenParticle* p = (const GenParticle*)branchGen->At(i);
      const int pid = p->PID, absId = abs(p->PID);
      if ( absId != 6 ) continue; // top only.
      if ( p->D1 == -1 or p->D2 == -1 ) continue; // should have valid daughters

      bool isDupl = false;
      for ( int j=p->D1; j<=p->D2; ++j ) {
        const GenParticle* dau  = (const GenParticle*)branchGen->At(j);
        if ( pid == dau->PID ) { isDupl = true; break; }
      }
      if ( isDupl ) continue;

      gen_topDaus.emplace_back();
      for ( int j=p->D1; j<=p->D2; ++j ) gen_topDaus.back().push_back(j);

      // Fill top quarks
      b_gen_pt[b_gen_n] = p->PT;
      b_gen_eta[b_gen_n] = p->Eta;
      b_gen_phi[b_gen_n] = p->Phi;
      b_gen_m[b_gen_n] = p->Mass;
      b_gen_pdgId[b_gen_n] = p->PID;
      b_gen_q3[b_gen_n] = p->Charge*3;
      b_gen_dau1[b_gen_n] = b_gen_dau2[b_gen_n] = -1;
      b_gen_mother[b_gen_n] = -1;

      ++b_gen_n;
      if ( b_gen_n >= gen_N ) break;
    }
    std::vector<int> dauIdx;
    for ( int i=0, n=gen_topDaus.size(); i<n; ++i ) {
      const auto& dauIdxs = gen_topDaus.at(i);
      if ( dauIdxs.empty() ) continue; // no daughters

      b_gen_dau1[i] = b_gen_n;
      b_gen_dau2[i] = b_gen_n+dauIdxs.size()-1;

      for ( auto j : dauIdxs ) {
        const GenParticle* dau = (const GenParticle*)branchGen->At(j);

        // Fill top quark daughters
        b_gen_pt[b_gen_n] = dau->PT;
        b_gen_eta[b_gen_n] = dau->Eta;
        b_gen_phi[b_gen_n] = dau->Phi;
        b_gen_m[b_gen_n] = dau->Mass;
        b_gen_pdgId[b_gen_n] = dau->PID;
        b_gen_q3[b_gen_n] = dau->Charge*3;
        b_gen_dau1[b_gen_n] = b_gen_dau2[b_gen_n] = -1;
        b_gen_mother[b_gen_n] = i;

        const int iDau = b_gen_n;
        ++b_gen_n;
        if ( b_gen_n >= gen_N ) break;

        if ( abs(dau->PID) < 23 or abs(dau->PID) > 25 ) continue;
        if ( dau->D1 == -1 or dau->D2 == -1 ) continue;
        dau = (const GenParticle*)branchGen->At(getLast(branchGen, j));

        int ngdau = 0;
        for ( int k=dau->D1; k<=dau->D2; ++k ) {
          const GenParticle* gdau = (const GenParticle*)branchGen->At(k);

          // Fill W/Z/H daughters
          b_gen_pt[b_gen_n] = gdau->PT;
          b_gen_eta[b_gen_n] = gdau->Eta;
          b_gen_phi[b_gen_n] = gdau->Phi;
          b_gen_m[b_gen_n] = gdau->Mass;
          b_gen_pdgId[b_gen_n] = gdau->PID;
          b_gen_q3[b_gen_n] = gdau->Charge*3;
          b_gen_dau1[b_gen_n] = b_gen_dau2[b_gen_n] = -1;
          b_gen_mother[b_gen_n] = iDau;

          ++ngdau;
          ++b_gen_n;
          if ( b_gen_n >= gen_N ) break;

          const int iGdau = b_gen_n;
          // For the tau decays
          if ( abs(gdau->PID) != 15 ) continue;
          if ( gdau->D1 == -1 or gdau->D2 == -1 ) continue;
          gdau = (const GenParticle*)branchGen->At(getLast(branchGen, k));

          int nggdau = 0;
          for ( int l=gdau->D1; l<=gdau->D2; ++l ) {
            const GenParticle* ggdau = (const GenParticle*)branchGen->At(l);
            const int absId = abs(ggdau->PID);
            if ( absId < 11 or absId >= 15 ) continue;

            // Fill W/Z/H daughters
            b_gen_pt[b_gen_n] = ggdau->PT;
            b_gen_eta[b_gen_n] = ggdau->Eta;
            b_gen_phi[b_gen_n] = ggdau->Phi;
            b_gen_m[b_gen_n] = ggdau->Mass;
            b_gen_pdgId[b_gen_n] = ggdau->PID;
            b_gen_q3[b_gen_n] = ggdau->Charge*3;
            b_gen_dau1[b_gen_n] = b_gen_dau2[b_gen_n] = -1;
            b_gen_mother[b_gen_n] = iGdau;

            ++nggdau;
            ++b_gen_n;
            if ( b_gen_n >= gen_N ) break;
          }
          b_gen_dau1[iGdau] = iGdau+1;
          b_gen_dau2[iGdau] = iGdau+nggdau;
          if ( b_gen_n >= gen_N ) break;
        }
        b_gen_dau1[iDau] = iDau+1;
        b_gen_dau2[iDau] = iDau+ngdau;
        
        if ( b_gen_n >= gen_N ) break;
      }
      if ( b_gen_n >= gen_N ) break;
    }

    if ( branchMET->GetEntries() > 0 ) {
      const MissingET* met = (const MissingET*)branchMET->At(0);
      b_met_pt = met->MET;
      b_met_phi = met->Phi;
    }
    else {
      b_met_pt = b_met_phi = 0;
    }

    b_muons_n = 0;
    for ( int i=0; i<branchMuon->GetEntries(); ++i ) {
      const Muon* muon = (const Muon*) branchMuon->At(i);
      const TLorentzVector p4 = muon->P4();

      b_muons_pt[b_muons_n] = muon->PT;
      b_muons_eta[b_muons_n] = muon->Eta;
      b_muons_phi[b_muons_n] = muon->Phi;
      b_muons_m[b_muons_n] = p4.M();
      b_muons_q[b_muons_n] = muon->Charge;

      b_muons_relIso[b_muons_n] = muon->IsolationVar;///muon->PT;

      ++b_muons_n;
      if ( b_muons_n >= muons_N ) break;
    }

    b_electrons_n = 0;
    for ( int i=0; i<branchElectron->GetEntries(); ++i ) {
      const Electron* electron = (const Electron*) branchElectron->At(i);
      const TLorentzVector p4 = electron->P4();

      b_electrons_pt[b_electrons_n] = electron->PT;
      b_electrons_eta[b_electrons_n] = electron->Eta;
      b_electrons_phi[b_electrons_n] = electron->Phi;
      b_electrons_m[b_electrons_n] = p4.M();
      b_electrons_q[b_electrons_n] = electron->Charge;

      b_electrons_relIso[b_electrons_n] = electron->IsolationVarRhoCorr;///electron->PT;

      ++b_electrons_n;
      if ( b_electrons_n >= electrons_N ) break;
    }

    b_jets_n = b_subjets_n = 0;
    for ( int i=0; i<branchJet->GetEntries(); ++i ) {
      const Jet* jet = (const Jet*) branchJet->At(i);
      //const TLorentzVector p4 = jet->P4();

      b_jets_pt[b_jets_n] = jet->PT;
      b_jets_eta[b_jets_n] = jet->Eta;
      b_jets_phi[b_jets_n] = jet->Phi;
      b_jets_m[b_jets_n] = jet->Mass;
      b_jets_flav[b_jets_n] = jet->Flavor;
      b_jets_bTag[b_jets_n] = jet->BTag;

      // Keep the subjet particles
      TRefArray cons = jet->Constituents;
      for ( int j=0; j<cons.GetEntriesFast(); ++j ) {
        if ( b_subjets_n > subjets_N ) break;

        const TObject* obj = cons.At(j);
        if ( !obj ) continue;

        //const GenParticle* p = dynamic_cast<const GenParticle*>(obj);
        const Track* track = dynamic_cast<const Track*>(obj);
        const Tower* tower = dynamic_cast<const Tower*>(obj);
        if ( track ) {
          b_subjets_pt[b_subjets_n] = track->PT;
          b_subjets_eta[b_subjets_n] = track->Eta;
          b_subjets_phi[b_subjets_n] = track->Phi;
          b_subjets_q[b_subjets_n] = track->Charge;
          b_subjets_pdgId[b_subjets_n] = track->Charge*211;
        }
        else if ( tower ) {
          b_subjets_pt[b_subjets_n] = tower->ET;
          b_subjets_eta[b_subjets_n] = tower->Eta;
          b_subjets_phi[b_subjets_n] = tower->Phi;
          b_subjets_q[b_subjets_n] = 0;
          const bool isPhoton = ( tower->Eem > tower->Ehad ); // Crude estimation
          if ( isPhoton ) b_subjets_pdgId[b_subjets_n] = 22; // photons
          else b_subjets_pdgId[b_subjets_n] = 2112; // set as neutron
        }
        else {
          std::cout << obj->IsA()->GetName() << endl;
          continue;
        }
        b_subjets_jetIdx[b_subjets_n] = b_jets_n;
        ++b_subjets_n;
      }

      ++b_jets_n;
      if ( b_jets_n >= jets_N ) break;
    }

    tree->Fill();
  }

  tree->Write();
  fout->Close();

}

