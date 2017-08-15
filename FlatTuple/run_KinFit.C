#include "AnalyzeKinFitHYTuple.h"
#include "AnalyzeKinFitDelphes.h"
R__LOAD_LIBRARY(TTLJKinFit.C+)
R__LOAD_LIBRARY(AnalyzeKinFitHYTuple.C+)
R__LOAD_LIBRARY(AnalyzeKinFitDelphes.C+)

const string mode = "FCNC"; // mode set to FCNC, require 3 b jets
//const string mode = "tt";

void run_CMSKinFit();
void run_DelphesKinFit();

void run_KinFit()
{
  //run_CMSKinFit();
  run_DelphesKinFit();
}

void run_CMSKinFit()
{
  TChain chainFCNC("ttbbLepJets/tree");
  chainFCNC.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct.root");
  chainFCNC.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
  ///home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut.root
  ///home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hut.root

  TChain chainTTBB("ttbbLepJets/tree");
  chainTTBB.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_powheg_ttbb.root");

  AnalyzeKinFitHYTuple tFCNC(&chainFCNC);
  AnalyzeKinFitHYTuple tTTBB(&chainTTBB);

  tFCNC.Loop(mode, "kin/cmsTuple_FCNC.root");
  tTTBB.Loop(mode, "kin/cmsTuple_ttbb.root");
}

void run_DelphesKinFit()
{
  TChain chainFCNC("tree");
  chainFCNC.Add("../Delphes2Flat/ntuple_tch.root");

  TChain chainTTBB("tree");
  chainTTBB.Add("../Delphes2Flat/ntuple_ttbb.root");

  TChain chainTTJJ("tree");
  chainTTJJ.Add("../Delphes2Flat/ntuple_tt.root");

  AnalyzeKinFitDelphes tFCNC(&chainFCNC);
  AnalyzeKinFitDelphes tTTBB(&chainTTBB);
  AnalyzeKinFitDelphes tTTJJ(&chainTTJJ);

  tFCNC.Loop(mode, "kin/delphes_FCNC.root");
  tTTBB.Loop(mode, "kin/delphes_ttbb.root");
  tTTJJ.Loop(mode, "kin/delphes_tt.root");
}
