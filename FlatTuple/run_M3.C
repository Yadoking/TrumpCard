#include "AnalyzeBasic.h"
R__LOAD_LIBRARY(AnalyzeBasic.C+)

void run_M3()
{
  TChain chainFCNC("tree");
  chainFCNC.Add("../Delphes2Flat/ntuple_tch.root");

  TChain chainTTBB("tree");
  chainTTBB.Add("../Delphes2Flat/ntuple_ttbb.root");

  AnalyzeBasic tFCNC(&chainFCNC);
  AnalyzeBasic tTTBB(&chainTTBB);

  tFCNC.Loop("delphes_FCNC.root");
  tTTBB.Loop("delphes_ttbb.root");
}
