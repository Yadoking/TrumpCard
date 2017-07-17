#include "AnalyzeM3HYTuple.h"
#include "AnalyzeM3Delphes.h"
R__LOAD_LIBRARY(AnalyzeM3HYTuple.C+)
R__LOAD_LIBRARY(AnalyzeM3Delphes.C+)

const string sample = "FCNC";
//const string sample = "ttbb";

void run_CMSM3();
void run_DelphesM3();

void run_M3()
{
  run_CMSM3();
  //run_DelphesM3();
}

void run_CMSM3()
{
  TChain chain("ttbbLepJets/tree");
  if ( sample == "FCNC" ) {
    chain.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hct.root");
    chain.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hct.root");
    ///home/minerva1993/fcnc/ntuple_jw/v2/TT_AntitopLeptonicDecay_TH_1L3B_Eta_Hut.root
    ///home/minerva1993/fcnc/ntuple_jw/v2/TT_TopLeptonicDecay_TH_1L3B_Eta_Hut.root
  }
  else if ( sample == "ttbb" ) {
    chain.Add("/home/minerva1993/fcnc/ntuple_jw/v2/TT_powheg_ttbb.root");
  }
  AnalyzeM3HYTuple t(&chain);

  t.Loop(Form("cmsTuple_%s.root", sample.c_str()));
}

void run_DelphesM3()
{
  TChain chain("tree");
  if ( sample == "FCNC" ) {
    chain.Add("../Delphes2Flat/ntuple_tch.root");
  }
  else if ( sample == "ttbb" ) {
    chain.Add("../Delphes2Flat/ntuple_ttbb.root");
  }
  AnalyzeM3Delphes t(&chain);

  t.Loop(Form("delphes_%s.root", sample.c_str()));
}
