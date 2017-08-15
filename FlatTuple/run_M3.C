#include "AnalyzeM3HYTuple.h"
#include "AnalyzeM3Delphes.h"
R__LOAD_LIBRARY(AnalyzeM3HYTuple.C+)
R__LOAD_LIBRARY(AnalyzeM3Delphes.C+)

const char* mode = "FCNC";
//const string mode = "tt";
//const char* algo = "m3";
const char* algo = "deltaR";

void run_CMSM3();
void run_DelphesM3();

void run_M3()
{
  //run_CMSM3();
  run_DelphesM3();
}

void run_CMSM3()
{
  TChain chainFCNC("ttbbLepJets/tree");
  chainFCNC.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct.root");
  chainFCNC.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
  ///home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut.root
  ///home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hut.root

  TChain chainTTBB("ttbbLepJets/tree");
  chainTTBB.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_powheg_ttbb.root");

  AnalyzeM3HYTuple tFCNC(&chainFCNC);
  AnalyzeM3HYTuple tTTBB(&chainTTBB);

  tFCNC.Loop(mode, Form("%s/cmsTuple_FCNC.root", algo), algo);
  tTTBB.Loop(mode, Form("%s/cmsTuple_ttbb.root", algo), algo);
}

void run_DelphesM3()
{
  TChain chainFCNC("tree");
  chainFCNC.Add("../Delphes2Flat/ntuple_tch.root");

  TChain chainTTBB("tree");
  chainTTBB.Add("../Delphes2Flat/ntuple_ttbb.root");

  TChain chainTTJJ("tree");
  chainTTJJ.Add("../Delphes2Flat/ntuple_ttbb.root");

  AnalyzeM3Delphes tFCNC(&chainFCNC);
  AnalyzeM3Delphes tTTBB(&chainTTBB);
  AnalyzeM3Delphes tTTJJ(&chainTTJJ);

  tFCNC.Loop(mode, Form("%s/delphes_FCNC.root", algo), algo);
  tTTBB.Loop(mode, Form("%s/delphes_ttbb.root", algo), algo);
  tTTJJ.Loop(mode, Form("%s/delphes_tt.root", algo), algo);
}
