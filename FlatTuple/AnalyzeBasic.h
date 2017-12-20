//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  5 19:33:50 2017 by ROOT version 6.08/02
// from TTree tree/tree
// found on file: ntuple_tch.root
//////////////////////////////////////////////////////////

#ifndef AnalyzeBasic_h
#define AnalyzeBasic_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include "TLorentzVector.h"

// Header file for the classes stored in the TTree if any.

class AnalyzeBasic {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UShort_t        run;
   UInt_t          event;
   Float_t         weight;
   Float_t         met_pt;
   Float_t         met_phi;
   UShort_t        muons_n;
   Float_t         muons_pt[100];   //[muons_n]
   Float_t         muons_eta[100];   //[muons_n]
   Float_t         muons_phi[100];   //[muons_n]
   Float_t         muons_m[100];   //[muons_n]
   Short_t         muons_q[100];   //[muons_n]
   Float_t         muons_relIso[100];   //[muons_n]
   UShort_t        electrons_n;
   Float_t         electrons_pt[100];   //[electrons_n]
   Float_t         electrons_eta[100];   //[electrons_n]
   Float_t         electrons_phi[100];   //[electrons_n]
   Float_t         electrons_m[100];   //[electrons_n]
   Short_t         electrons_q[100];   //[electrons_n]
   Float_t         electrons_relIso[100];   //[electrons_n]
   UShort_t        jets_n;
   Float_t         jets_pt[100];   //[jets_n]
   Float_t         jets_eta[100];   //[jets_n]
   Float_t         jets_phi[100];   //[jets_n]
   Float_t         jets_m[100];   //[jets_n]
   Short_t         jets_flav[100];   //[jets_n]
   Float_t         jets_bTag[100];   //[jets_n]
   UShort_t        gen_n;
   Float_t         gen_pt[100];   //[gen_n]
   Float_t         gen_eta[100];   //[gen_n]
   Float_t         gen_phi[100];   //[gen_n]
   Float_t         gen_m[100];   //[gen_n]
   Short_t         gen_pdgId[100];   //[gen_n]
   Short_t         gen_q3[100];   //[gen_n]
   Short_t         gen_mother[100];   //[gen_n]
   Short_t         gen_dau1[100];   //[gen_n]
   Short_t         gen_dau2[100];   //[gen_n]
   UShort_t        subjets_n;
   Float_t         subjets_pt[1000];
   Float_t         subjets_eta[1000];
   Float_t         subjets_phi[1000];
   Short_t         subjets_q[1000];
   Short_t         subjets_pdgId[1000];
   UShort_t        subjets_jetIdx[1000];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_muons_n;   //!
   TBranch        *b_muons_pt;   //!
   TBranch        *b_muons_eta;   //!
   TBranch        *b_muons_phi;   //!
   TBranch        *b_muons_m;   //!
   TBranch        *b_muons_q;   //!
   TBranch        *b_muons_relIso;   //!
   TBranch        *b_electrons_n;   //!
   TBranch        *b_electrons_pt;   //!
   TBranch        *b_electrons_eta;   //!
   TBranch        *b_electrons_phi;   //!
   TBranch        *b_electrons_m;   //!
   TBranch        *b_electrons_q;   //!
   TBranch        *b_electrons_relIso;   //!
   TBranch        *b_jets_n;   //!
   TBranch        *b_jets_pt;   //!
   TBranch        *b_jets_eta;   //!
   TBranch        *b_jets_phi;   //!
   TBranch        *b_jets_m;   //!
   TBranch        *b_jets_flav;   //!
   TBranch        *b_jets_bTag;   //!
   TBranch        *b_gen_n;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_m;   //!
   TBranch        *b_gen_pdgId;   //!
   TBranch        *b_gen_q3;   //!
   TBranch        *b_gen_mother;   //!
   TBranch        *b_gen_dau1;   //!
   TBranch        *b_gen_dau2;   //!
   TBranch        *b_subjets_n;
   TBranch        *b_subjets_pt;
   TBranch        *b_subjets_eta;
   TBranch        *b_subjets_phi;
   TBranch        *b_subjets_q;
   TBranch        *b_subjets_pdgId;
   TBranch        *b_subjets_jetIdx;

   AnalyzeBasic(TTree *tree=0);
   virtual ~AnalyzeBasic();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const std::string outFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   struct TTLJSolution {
     TTLJSolution() { isValid = false; }
     bool isValid;
     TLorentzVector lep, met, lepB;
     TLorentzVector hadJ1, hadJ2, hadJ3;
     size_t lepBIdx, hadJ1Idx, hadJ2Idx, hadJ3Idx;
     double quality;
   };

   TTLJSolution solveByDeltaR(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const;
   TTLJSolution solveByM3(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const;
   TTLJSolution solveByMTop(TLorentzVector lepP4, TLorentzVector metP4, std::vector<size_t> jetIdxs) const;

   void rotate(float& phi, const float refPhi) {
     phi -= refPhi;
     const int n2Pi = floor(phi/2/TMath::Pi());
     phi -= n2Pi*(2*TMath::Pi());
     if ( phi > TMath::Pi() ) phi -= 2*TMath::Pi();
   }

   double deltaPhi(const double phi1, const double phi2)
   {
     double dphi = phi1 - phi2;
     while ( dphi < -TMath::Pi() ) dphi += 2*TMath::Pi();
     while ( dphi >  TMath::Pi() ) dphi -= 2*TMath::Pi();
     return dphi;
   }
};

#endif

#ifdef AnalyzeBasic_cxx
AnalyzeBasic::AnalyzeBasic(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple_tch.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ntuple_tch.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

AnalyzeBasic::~AnalyzeBasic()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeBasic::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeBasic::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalyzeBasic::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("MET_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("MET_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("nMuon", &muons_n, &b_muons_n);
   fChain->SetBranchAddress("Muon_pt", muons_pt, &b_muons_pt);
   fChain->SetBranchAddress("Muon_eta", muons_eta, &b_muons_eta);
   fChain->SetBranchAddress("Muon_phi", muons_phi, &b_muons_phi);
   fChain->SetBranchAddress("Muon_m", muons_m, &b_muons_m);
   fChain->SetBranchAddress("Muon_q", muons_q, &b_muons_q);
   fChain->SetBranchAddress("Muon_relIso", muons_relIso, &b_muons_relIso);
   fChain->SetBranchAddress("nElectron", &electrons_n, &b_electrons_n);
   fChain->SetBranchAddress("Electron_pt", electrons_pt, &b_electrons_pt);
   fChain->SetBranchAddress("Electron_eta", electrons_eta, &b_electrons_eta);
   fChain->SetBranchAddress("Electron_phi", electrons_phi, &b_electrons_phi);
   fChain->SetBranchAddress("Electron_m", electrons_m, &b_electrons_m);
   fChain->SetBranchAddress("Electron_q", electrons_q, &b_electrons_q);
   fChain->SetBranchAddress("Electron_relIso", electrons_relIso, &b_electrons_relIso);
   fChain->SetBranchAddress("nJet", &jets_n, &b_jets_n);
   fChain->SetBranchAddress("Jet_pt", jets_pt, &b_jets_pt);
   fChain->SetBranchAddress("Jet_eta", jets_eta, &b_jets_eta);
   fChain->SetBranchAddress("Jet_phi", jets_phi, &b_jets_phi);
   fChain->SetBranchAddress("Jet_m", jets_m, &b_jets_m);
   fChain->SetBranchAddress("Jet_flav", jets_flav, &b_jets_flav);
   fChain->SetBranchAddress("Jet_bTag", jets_bTag, &b_jets_bTag);
   fChain->SetBranchAddress("nGenParticle", &gen_n, &b_gen_n);
   fChain->SetBranchAddress("GenParticle_pt", gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("GenParticle_eta", gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("GenParticle_phi", gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("GenParticle_m", gen_m, &b_gen_m);
   fChain->SetBranchAddress("GenParticle_pdgId", gen_pdgId, &b_gen_pdgId);
   fChain->SetBranchAddress("GenParticle_q3", gen_q3, &b_gen_q3);
   fChain->SetBranchAddress("GenParticle_mother", gen_mother, &b_gen_mother);
   fChain->SetBranchAddress("GenParticle_dau1", gen_dau1, &b_gen_dau1);
   fChain->SetBranchAddress("GenParticle_dau2", gen_dau2, &b_gen_dau2);
   fChain->SetBranchAddress("nSubjet", &subjets_n, &b_subjets_n);
   fChain->SetBranchAddress("Subjet_pt", &subjets_pt, &b_subjets_pt);
   fChain->SetBranchAddress("Subjet_eta", &subjets_eta, &b_subjets_eta);
   fChain->SetBranchAddress("Subjet_phi", &subjets_phi, &b_subjets_phi);
   fChain->SetBranchAddress("Subjet_q", &subjets_q, &b_subjets_q);
   fChain->SetBranchAddress("Subjet_pdgId", &subjets_pdgId, &b_subjets_pdgId);
   fChain->SetBranchAddress("Subjet_jetIdx", &subjets_jetIdx, &b_subjets_jetIdx);
   Notify();
}

Bool_t AnalyzeBasic::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzeBasic::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeBasic::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzeBasic_cxx
