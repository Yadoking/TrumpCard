#define INPUTTYPE 1

#if ( INPUTTYPE == 1 )
#include "AnalyzerHYTuple.h"
R__LOAD_LIBRARY(TTLJKinFit.C+)
R__LOAD_LIBRARY(AnalyzerHYTuple.C+)
#else
#include "AnalyzerDelphes.h"
R__LOAD_LIBRARY(TTLJKinFit.C+)
R__LOAD_LIBRARY(AnalyzerDelphes.C+)
#endif

void run_HYTuple()
{
  const string mode = "FCNC";
  //const string mode = "tt";
  const string sample = "FCNC";
  //const string sample = "ttbb";
  
#if ( INPUTTYPE == 1 )
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
  AnalyzerHYTuple t(&chain);
#else
  TChain chain("tree");
  if ( sample == "FCNC" ) {
    chain.Add("../ntuple_tch.root");
  }
  else if ( sample == "ttbb" ) {
    chain.Add("../ntuple_ttbb.root");
  }
  AnalyzerDelphes t(&chain);
#endif

  t.Loop(mode, Form("ntuple_%s.root", sample.c_str())); // mode set to FCNC, require 3 b jets
}
