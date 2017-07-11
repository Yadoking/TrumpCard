#define Yosi_cxx
#include "Yosi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Yosi::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Yosi.C
//      root> Yosi t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   TH1F* hdR = new TH1F("dR", "dR", 100, 0, 5);
   TH1F* hdijetM = new TH1F("dijet_Mass", "dijet Mass", 300, 0, 300);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      //
      // Yosi!Programming Season!
      int njets = 0;
      int nbjets = 0;

      vector<TLorentzVector> p4s;

      for (int i = 0; i < jets_n; ++i)
      {
         if (jets_pt[i] < 30 or std::abs(jets_eta[i]) > 2.5) continue;
         njets++;
         if (jets_bTag[i] == 1)
         {
            ++nbjets;
            TLorentzVector p4;
            p4.SetPtEtaPhiM(jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
            p4s.push_back(p4);
         }
      }

/*
      double dr = 0;
      double diM = 0;
      vector<double> min;
      vector<double> dijetM;
      TLorentzVector Bjet1, Bjet2;
*/

      if (nbjets > 1)
      {
         double minDR = 999, minM = -1;
         for (int j = 1; j < nbjets; ++j)
         {
            const auto& Bjet1 = p4s.at(j);
            for (int k = 0; k < j; ++k)
            {
               const auto& Bjet2 = p4s.at(k);
               const double dR = Bjet1.DeltaR(Bjet2);
               if ( dR < minDR ) {
                 minDR = dR;
                 minM = (Bjet1+Bjet2).M();
               }
            }
         }
/*
vector<pair<double, double> > mins;
std::sort(mins.begin(), mins.end(), 
   [&](  pair<double, double>a, pair<double, double> b)
   {
       return a.first < b.first; 
   }
);
*/

	 hdR->Fill(minDR);
	 hdijetM->Fill(minM);
      }
   }
   
   TCanvas *c = new TCanvas("cDr", "delta R", 500, 500);
   hdR->Draw();
   c = new TCanvas("cdiM", "dijet Mass", 500, 500);
   hdijetM->Draw();
}
