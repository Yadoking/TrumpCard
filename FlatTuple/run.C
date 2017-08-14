#include "Yosi.h"
R__LOAD_LIBRARY(Yosi.C+)

void run()
{
   gSystem->CompileMacro("Yosi.C");

   TChain chain("tree");
   //chain.Add("out.root");
   chain.Add("../Delphes2Flat/ntuple_tch.root");

   Yosi t(&chain);
   t.Loop();
}
