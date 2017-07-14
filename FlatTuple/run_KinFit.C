#include "AnalyzerHYTuple.h"
#include "AnalyzerDelphes.h"
R__LOAD_LIBRARY(TTLJKinFit.C+)
R__LOAD_LIBRARY(AnalyzerHYTuple.C+)
R__LOAD_LIBRARY(AnalyzerDelphes.C+)

const string mode = "FCNC";
//const string mode = "tt";
const string sample = "FCNC";
//const string sample = "ttbb";

void run_CMSAnalyzer();
void run_DelphesAnalyzer();

void run_KinFit()
{
  //run_CMSAnalyzer();
  run_DelphesAnalyzer();
}

void run_CMSAnalyzer()
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
  AnalyzerHYTuple t(&chain);

  t.Loop(mode, Form("cmsTuple_%s.root", sample.c_str())); // mode set to FCNC, require 3 b jets
}

void run_DelphesAnalyzer()
{
  TChain chain("tree");
  if ( sample == "FCNC" ) {
    chain.Add("../ntuple_tch.root");
  }
  else if ( sample == "ttbb" ) {
    chain.Add("../ntuple_ttbb.root");
  }
  AnalyzerDelphes t(&chain);

  t.Loop(mode, Form("delphes_%s.root", sample.c_str())); // mode set to FCNC, require 3 b jets
}
