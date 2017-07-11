#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TTree.h"
#include "TFile.h"
#endif

//------------------------------------------------------------------------------
const std::vector<std::string> inputFiles = {
  //"/home/minerva1993/public/delphes_analysis/tch_run01.root",
  "/home/minerva1993/public/delphes_analysis/ttbb_run02.root",
};
//const std::string outputFile = "ntuple_tch.root";
const std::string outputFile = "ntuple_ttbb.root";

void makeFlatTuple()
{
  gSystem->Load("libDelphes");

  // Prepare output tree
  TFile* fout = TFile::Open(outputFile.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");

  unsigned short b_run = 1;
  unsigned long b_event = 0;
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

  tree->Branch("run", &b_run, "run/s");
  tree->Branch("event", &b_event, "event/i");

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

  // Create chain of root trees
  TChain chain("Delphes");
  for ( auto& inputFile : inputFiles ) chain.Add(inputFile.c_str());

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    cout << entry << '/' << numberOfEntries << '\r';

    const HepMCEvent* event = (const HepMCEvent*)branchEvent->At(0);
    b_event = event->Number;
    b_weight = event->Weight;

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

    b_jets_n = 0;
    for ( int i=0; i<branchJet->GetEntries(); ++i ) {
      const Jet* jet = (const Jet*) branchJet->At(i);
      //const TLorentzVector p4 = jet->P4();

      b_jets_pt[b_jets_n] = jet->PT;
      b_jets_eta[b_jets_n] = jet->Eta;
      b_jets_phi[b_jets_n] = jet->Phi;
      b_jets_m[b_jets_n] = jet->Mass;
      b_jets_flav[b_jets_n] = jet->Flavor;
      b_jets_bTag[b_jets_n] = jet->BTag;

      ++b_jets_n;
      if ( b_jets_n >= jets_N ) break;
    }

    tree->Fill();
  }

  tree->Write();
  fout->Close();
}

