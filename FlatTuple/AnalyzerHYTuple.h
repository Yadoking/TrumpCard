//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 10 22:12:13 2017 by ROOT version 6.08/02
// from TTree tree/TopTree
// found on file: /home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root
//////////////////////////////////////////////////////////

#ifndef AnalyzerHYTuple_h
#define AnalyzerHYTuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class AnalyzerHYTuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           luminumber;
   Float_t         genweight;
   Int_t           GoodPV;
   Int_t           channel;
   vector<float>   *PUWeight;
   Float_t         MET;
   Float_t         MET_phi;
   Float_t         lepton_pT;
   Float_t         lepton_eta;
   Float_t         lepton_phi;
   Float_t         lepton_E;
   Float_t         lepton_LES;
   vector<float>   *lepton_SF;
   Float_t         lepton_relIso;
   Bool_t          lepton_isIso;
   vector<float>   *jet_pT;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_E;
   vector<int>     *jet_index;
   vector<int>     *jet_gencone_mom;
   vector<float>   *jet_CSV;
   vector<float>   *jet_SF_CSV_25;
   vector<float>   *jet_SF_CSV_30;
   vector<float>   *jet_SF_CSV_35;
   vector<float>   *jet_SF_CSV_40;
   vector<float>   *jet_SF_CSV;
   vector<float>   *jet_CvsL;
   vector<float>   *jet_CvsB;
   Int_t           jet_number;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_JES_Up;
   vector<float>   *jet_JES_Down;
   vector<float>   *jet_JER_Up;
   vector<float>   *jet_JER_Nom;
   vector<float>   *jet_JER_Down;
   Float_t         kin_chi2;
   Float_t         kinnu_pT;
   Float_t         kinnu_eta;
   Float_t         kinnu_phi;
   Float_t         kinnu_E;
   vector<float>   *kinjet_pT;
   vector<float>   *kinjet_eta;
   vector<float>   *kinjet_phi;
   vector<float>   *kinjet_E;
   vector<int>     *kinjet_index;
   Float_t         fcnhkin_chi2;
   Float_t         fcnhkinnu_pT;
   Float_t         fcnhkinnu_eta;
   Float_t         fcnhkinnu_phi;
   Float_t         fcnhkinnu_E;
   vector<float>   *fcnhkinjet_pT;
   vector<float>   *fcnhkinjet_eta;
   vector<float>   *fcnhkinjet_phi;
   vector<float>   *fcnhkinjet_E;
   vector<int>     *fcnhkinjet_index;
   Float_t         fcnhkinnu_M;
   Float_t         kinnu_M;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_luminumber;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_GoodPV;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_lepton_pT;   //!
   TBranch        *b_lepton_eta;   //!
   TBranch        *b_lepton_phi;   //!
   TBranch        *b_lepton_E;   //!
   TBranch        *b_lepton_LES;   //!
   TBranch        *b_lepton_SF;   //!
   TBranch        *b_lepton_relIso;   //!
   TBranch        *b_lepton_isIso;   //!
   TBranch        *b_jet_pT;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_index;   //!
   TBranch        *b_jet_gencone_mom;   //!
   TBranch        *b_jet_CSV;   //!
   TBranch        *b_jet_SF_CSV_25;   //!
   TBranch        *b_jet_SF_CSV_30;   //!
   TBranch        *b_jet_SF_CSV_35;   //!
   TBranch        *b_jet_SF_CSV_40;   //!
   TBranch        *b_jet_SF_CSV;   //!
   TBranch        *b_jet_CvsL;   //!
   TBranch        *b_jet_CvsB;   //!
   TBranch        *b_jet_number;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_JES_Up;   //!
   TBranch        *b_jet_JES_Down;   //!
   TBranch        *b_jet_JER_Up;   //!
   TBranch        *b_jet_JER_Nom;   //!
   TBranch        *b_jet_JER_Down;   //!
   TBranch        *b_kin_chi2;   //!
   TBranch        *b_kinnu_pT;   //!
   TBranch        *b_kinnu_eta;   //!
   TBranch        *b_kinnu_phi;   //!
   TBranch        *b_kinnu_E;   //!
   TBranch        *b_kinjet_pT;   //!
   TBranch        *b_kinjet_eta;   //!
   TBranch        *b_kinjet_phi;   //!
   TBranch        *b_kinjet_E;   //!
   TBranch        *b_kinjet_index;   //!
   TBranch        *b_fcnhkin_chi2;   //!
   TBranch        *b_fcnhkinnu_pT;   //!
   TBranch        *b_fcnhkinnu_eta;   //!
   TBranch        *b_fcnhkinnu_phi;   //!
   TBranch        *b_fcnhkinnu_E;   //!
   TBranch        *b_fcnhkinjet_pT;   //!
   TBranch        *b_fcnhkinjet_eta;   //!
   TBranch        *b_fcnhkinjet_phi;   //!
   TBranch        *b_fcnhkinjet_E;   //!
   TBranch        *b_fcnhkinjet_index;   //!
   TBranch        *b_fcnhkinnu_M;   //!
   TBranch        *b_kinnu_M;   //!

   AnalyzerHYTuple(TTree *tree=0);
   virtual ~AnalyzerHYTuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const std::string modeStr, const std::string outFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyzerHYTuple_cxx
AnalyzerHYTuple::AnalyzerHYTuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root:/ttbbLepJets");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

AnalyzerHYTuple::~AnalyzerHYTuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzerHYTuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzerHYTuple::LoadTree(Long64_t entry)
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

void AnalyzerHYTuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PUWeight = 0;
   lepton_SF = 0;
   jet_pT = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_E = 0;
   jet_index = 0;
   jet_gencone_mom = 0;
   jet_CSV = 0;
   jet_SF_CSV_25 = 0;
   jet_SF_CSV_30 = 0;
   jet_SF_CSV_35 = 0;
   jet_SF_CSV_40 = 0;
   jet_SF_CSV = 0;
   jet_CvsL = 0;
   jet_CvsB = 0;
   jet_partonFlavour = 0;
   jet_hadronFlavour = 0;
   jet_JES_Up = 0;
   jet_JES_Down = 0;
   jet_JER_Up = 0;
   jet_JER_Nom = 0;
   jet_JER_Down = 0;
   kinjet_pT = 0;
   kinjet_eta = 0;
   kinjet_phi = 0;
   kinjet_E = 0;
   kinjet_index = 0;
   fcnhkinjet_pT = 0;
   fcnhkinjet_eta = 0;
   fcnhkinjet_phi = 0;
   fcnhkinjet_E = 0;
   fcnhkinjet_index = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminumber", &luminumber, &b_luminumber);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("GoodPV", &GoodPV, &b_GoodPV);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("lepton_pT", &lepton_pT, &b_lepton_pT);
   fChain->SetBranchAddress("lepton_eta", &lepton_eta, &b_lepton_eta);
   fChain->SetBranchAddress("lepton_phi", &lepton_phi, &b_lepton_phi);
   fChain->SetBranchAddress("lepton_E", &lepton_E, &b_lepton_E);
   fChain->SetBranchAddress("lepton_LES", &lepton_LES, &b_lepton_LES);
   fChain->SetBranchAddress("lepton_SF", &lepton_SF, &b_lepton_SF);
   fChain->SetBranchAddress("lepton_relIso", &lepton_relIso, &b_lepton_relIso);
   fChain->SetBranchAddress("lepton_isIso", &lepton_isIso, &b_lepton_isIso);
   fChain->SetBranchAddress("jet_pT", &jet_pT, &b_jet_pT);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_index", &jet_index, &b_jet_index);
   fChain->SetBranchAddress("jet_gencone_mom", &jet_gencone_mom, &b_jet_gencone_mom);
   fChain->SetBranchAddress("jet_CSV", &jet_CSV, &b_jet_CSV);
   fChain->SetBranchAddress("jet_SF_CSV_25", &jet_SF_CSV_25, &b_jet_SF_CSV_25);
   fChain->SetBranchAddress("jet_SF_CSV_30", &jet_SF_CSV_30, &b_jet_SF_CSV_30);
   fChain->SetBranchAddress("jet_SF_CSV_35", &jet_SF_CSV_35, &b_jet_SF_CSV_35);
   fChain->SetBranchAddress("jet_SF_CSV_40", &jet_SF_CSV_40, &b_jet_SF_CSV_40);
   fChain->SetBranchAddress("jet_SF_CSV", &jet_SF_CSV, &b_jet_SF_CSV);
   fChain->SetBranchAddress("jet_CvsL", &jet_CvsL, &b_jet_CvsL);
   fChain->SetBranchAddress("jet_CvsB", &jet_CvsB, &b_jet_CvsB);
   fChain->SetBranchAddress("jet_number", &jet_number, &b_jet_number);
   fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
   fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
   fChain->SetBranchAddress("jet_JES_Up", &jet_JES_Up, &b_jet_JES_Up);
   fChain->SetBranchAddress("jet_JES_Down", &jet_JES_Down, &b_jet_JES_Down);
   fChain->SetBranchAddress("jet_JER_Up", &jet_JER_Up, &b_jet_JER_Up);
   fChain->SetBranchAddress("jet_JER_Nom", &jet_JER_Nom, &b_jet_JER_Nom);
   fChain->SetBranchAddress("jet_JER_Down", &jet_JER_Down, &b_jet_JER_Down);
   fChain->SetBranchAddress("kin_chi2", &kin_chi2, &b_kin_chi2);
   fChain->SetBranchAddress("kinnu_pT", &kinnu_pT, &b_kinnu_pT);
   fChain->SetBranchAddress("kinnu_eta", &kinnu_eta, &b_kinnu_eta);
   fChain->SetBranchAddress("kinnu_phi", &kinnu_phi, &b_kinnu_phi);
   fChain->SetBranchAddress("kinnu_E", &kinnu_E, &b_kinnu_E);
   fChain->SetBranchAddress("kinjet_pT", &kinjet_pT, &b_kinjet_pT);
   fChain->SetBranchAddress("kinjet_eta", &kinjet_eta, &b_kinjet_eta);
   fChain->SetBranchAddress("kinjet_phi", &kinjet_phi, &b_kinjet_phi);
   fChain->SetBranchAddress("kinjet_E", &kinjet_E, &b_kinjet_E);
   fChain->SetBranchAddress("kinjet_index", &kinjet_index, &b_kinjet_index);
   fChain->SetBranchAddress("fcnhkin_chi2", &fcnhkin_chi2, &b_fcnhkin_chi2);
   fChain->SetBranchAddress("fcnhkinnu_pT", &fcnhkinnu_pT, &b_fcnhkinnu_pT);
   fChain->SetBranchAddress("fcnhkinnu_eta", &fcnhkinnu_eta, &b_fcnhkinnu_eta);
   fChain->SetBranchAddress("fcnhkinnu_phi", &fcnhkinnu_phi, &b_fcnhkinnu_phi);
   fChain->SetBranchAddress("fcnhkinnu_E", &fcnhkinnu_E, &b_fcnhkinnu_E);
   fChain->SetBranchAddress("fcnhkinjet_pT", &fcnhkinjet_pT, &b_fcnhkinjet_pT);
   fChain->SetBranchAddress("fcnhkinjet_eta", &fcnhkinjet_eta, &b_fcnhkinjet_eta);
   fChain->SetBranchAddress("fcnhkinjet_phi", &fcnhkinjet_phi, &b_fcnhkinjet_phi);
   fChain->SetBranchAddress("fcnhkinjet_E", &fcnhkinjet_E, &b_fcnhkinjet_E);
   fChain->SetBranchAddress("fcnhkinjet_index", &fcnhkinjet_index, &b_fcnhkinjet_index);
   fChain->SetBranchAddress("fcnhkinnu_M", &fcnhkinnu_M, &b_fcnhkinnu_M);
   fChain->SetBranchAddress("kinnu_M", &kinnu_M, &b_kinnu_M);
   Notify();
}

Bool_t AnalyzerHYTuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzerHYTuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzerHYTuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzerHYTuple_cxx
