//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  5 19:33:50 2017 by ROOT version 6.08/02
// from TTree tree/tree
// found on file: ntuple_tch.root
//////////////////////////////////////////////////////////

#ifndef AnalyzeM3Delphes_h
#define AnalyzeM3Delphes_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"

// Header file for the classes stored in the TTree if any.

class AnalyzeM3Delphes {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UShort_t        run;
   UInt_t          event;
   Float_t         met_pt;
   Float_t         met_phi;
   UShort_t        muons_n;
   Float_t         muons_pt[3];   //[muons_n]
   Float_t         muons_eta[3];   //[muons_n]
   Float_t         muons_phi[3];   //[muons_n]
   Float_t         muons_m[3];   //[muons_n]
   Short_t         muons_q[3];   //[muons_n]
   Float_t         muons_relIso[3];   //[muons_n]
   UShort_t        electrons_n;
   Float_t         electrons_pt[3];   //[electrons_n]
   Float_t         electrons_eta[3];   //[electrons_n]
   Float_t         electrons_phi[3];   //[electrons_n]
   Float_t         electrons_m[3];   //[electrons_n]
   Short_t         electrons_q[3];   //[electrons_n]
   Float_t         electrons_relIso[3];   //[electrons_n]
   UShort_t        jets_n;
   Float_t         jets_pt[17];   //[jets_n]
   Float_t         jets_eta[17];   //[jets_n]
   Float_t         jets_phi[17];   //[jets_n]
   Float_t         jets_m[17];   //[jets_n]
   Short_t         jets_flav[17];   //[jets_n]
   Float_t         jets_bTag[17];   //[jets_n]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
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

   AnalyzeM3Delphes(TTree *tree=0);
   virtual ~AnalyzeM3Delphes();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const std::string modeStr, const std::string outFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void rotate(float& phi, const float refPhi) {
     phi -= refPhi;
     const int n2Pi = floor(phi/2/TMath::Pi());
     phi -= n2Pi*(2*TMath::Pi());
     if ( phi > TMath::Pi() ) phi -= 2*TMath::Pi();
   }
};

#endif

#ifdef AnalyzeM3Delphes_cxx
AnalyzeM3Delphes::AnalyzeM3Delphes(TTree *tree) : fChain(0)
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

AnalyzeM3Delphes::~AnalyzeM3Delphes()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeM3Delphes::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeM3Delphes::LoadTree(Long64_t entry)
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

void AnalyzeM3Delphes::Init(TTree *tree)
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
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("muons_n", &muons_n, &b_muons_n);
   fChain->SetBranchAddress("muons_pt", muons_pt, &b_muons_pt);
   fChain->SetBranchAddress("muons_eta", muons_eta, &b_muons_eta);
   fChain->SetBranchAddress("muons_phi", muons_phi, &b_muons_phi);
   fChain->SetBranchAddress("muons_m", muons_m, &b_muons_m);
   fChain->SetBranchAddress("muons_q", muons_q, &b_muons_q);
   fChain->SetBranchAddress("muons_relIso", muons_relIso, &b_muons_relIso);
   fChain->SetBranchAddress("electrons_n", &electrons_n, &b_electrons_n);
   fChain->SetBranchAddress("electrons_pt", electrons_pt, &b_electrons_pt);
   fChain->SetBranchAddress("electrons_eta", electrons_eta, &b_electrons_eta);
   fChain->SetBranchAddress("electrons_phi", electrons_phi, &b_electrons_phi);
   fChain->SetBranchAddress("electrons_m", electrons_m, &b_electrons_m);
   fChain->SetBranchAddress("electrons_q", electrons_q, &b_electrons_q);
   fChain->SetBranchAddress("electrons_relIso", electrons_relIso, &b_electrons_relIso);
   fChain->SetBranchAddress("jets_n", &jets_n, &b_jets_n);
   fChain->SetBranchAddress("jets_pt", jets_pt, &b_jets_pt);
   fChain->SetBranchAddress("jets_eta", jets_eta, &b_jets_eta);
   fChain->SetBranchAddress("jets_phi", jets_phi, &b_jets_phi);
   fChain->SetBranchAddress("jets_m", jets_m, &b_jets_m);
   fChain->SetBranchAddress("jets_flav", jets_flav, &b_jets_flav);
   fChain->SetBranchAddress("jets_bTag", jets_bTag, &b_jets_bTag);
   Notify();
}

Bool_t AnalyzeM3Delphes::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzeM3Delphes::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeM3Delphes::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzeM3Delphes_cxx
