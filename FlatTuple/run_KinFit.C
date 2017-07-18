#include "AnalyzeKinFitHYTuple.h"
#include "AnalyzeKinFitDelphes.h"
//R__LOAD_LIBRARY(TTLJKinFit.C+)
R__LOAD_LIBRARY(AnalyzeKinFitHYTuple.C+)
R__LOAD_LIBRARY(AnalyzeKinFitDelphes.C+)

const string mode = "FCNC";
//const string mode = "tt";
const string sample = "FCNC";
//const string sample = "ttbb";

void run_CMSKinFit();
void run_DelphesKinFit();

void run_KinFit()
{
  run_CMSKinFit();
  //run_DelphesKinFit();
}

void run_CMSKinFit()
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
  AnalyzeKinFitHYTuple t(&chain);

  t.Loop(mode, Form("kin/cmsTuple_%s.root", sample.c_str())); // mode set to FCNC, require 3 b jets
}

void run_DelphesKinFit()
{
  TChain chain("tree");
  if ( sample == "FCNC" ) {
    chain.Add("../Delphes2Flat/ntuple_tch.root");
  }
  else if ( sample == "ttbb" ) {
    chain.Add("../Delphes2Flat/ntuple_ttbb.root");
  }
  AnalyzeKinFitDelphes t(&chain);

  t.Loop(mode, Form("kin/delphes_%s.root", sample.c_str())); // mode set to FCNC, require 3 b jets
}
