#include "Yosi.h"
R__LOAD_LIBRARY(Yosi.C+)

void run()
{
   gSystem->CompileMacro("Yosi.C");

   TChain chain("tree");
   //chain.Add("out.root");
   chain.Add("out_fcnc.root");

   Yosi t(&chain);
   t.Loop();
}
