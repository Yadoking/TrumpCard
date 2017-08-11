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

	TFile* fout = new TFile("temp.root", "recreate");
	TTree* tree = new TTree("namu", "namu");

	int b_njet;
	int b_nbjet;

	float b_lepton_pt;
	float b_lepton_eta;
	float b_lepton_phi;
	float b_lepton_m;

	float b_jet_pt;
	float b_jet_eta;
	float b_jet_phi;
	float b_jet_m;

	float b_m3_m;
	float b_m3_m12;
	float b_m3_m23;
	float b_m3_m31;
	float b_m3_dR12;
	float b_m3_dR23;
	float b_m3_dR31;

	float b_dr_m;
	float b_dr_m12;
	float b_dr_m23;
	float b_dr_m31;
	float b_dr_dR12;
	float b_dr_dR23;
	float b_dr_dR31;

	tree->Branch("njet", &b_njet, "njet/I");
	tree->Branch("nbjet", &b_nbjet, "nbjet/I");

	tree->Branch("lepton_pt", &b_lepton_pt, "lepton_pt/F");
	tree->Branch("lepton_eta", &b_lepton_eta, "lepton_eta/F");
	tree->Branch("lepton_phi", &b_lepton_phi, "lepton_phi/F");
	tree->Branch("lepton_m", &b_lepton_m, "lepton_m/F");

	tree->Branch("jet_pt", &b_jet_pt, "jet_pt/F");
	tree->Branch("jet_eta", &b_jet_eta, "jet_eta/F");
	tree->Branch("jet_phi", &b_jet_phi, "jet_phi/F");
	tree->Branch("jet_m", &b_jet_m, "jet_m/F");

	tree->Branch("M3_Massj1j2j3", &b_m3_m, "M3_Massj1j2j3/F");
	tree->Branch("M3_Massj1j2", &b_m3_m12, "M3_Massj1j2/F");
	tree->Branch("M3_Massj2j3", &b_m3_m23, "M3_Massj2j3/F");
	tree->Branch("M3_Massj3j1", &b_m3_m31, "M3_Massj3j1/F");
	tree->Branch("M3_dRj1j2", &b_m3_dR12, "M3_dRj1j2/F");
	tree->Branch("M3_dRj2j3", &b_m3_dR23, "M3_dRj2j3/F");
	tree->Branch("M3_dRj3j1", &b_m3_dR31, "M3_dRj3j1/F");

	tree->Branch("dR_Massj1j2j3", &b_dr_m, "dR_Massj1j2j3/F");
	tree->Branch("dR_Massj1j2", &b_dr_m12, "dR_Massj1j2/F");
	tree->Branch("dR_Massj2j3", &b_dr_m23, "dR_Massj2j3/F");
	tree->Branch("dR_Massj3j1", &b_dr_m31, "dR_Massj3j1/F");
	tree->Branch("dR_dRj1j2", &b_dr_dR12, "dR_dRj1j2/F");
	tree->Branch("dR_dRj2j3", &b_dr_dR23, "dR_dRj2j3/F");
	tree->Branch("dR_dRj3j1", &b_dr_dR31, "dR_dRj3j1/F");


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

		// Information about lepton (muon or electron)
		TLorentzVector leptonP4;
		if ( muons_n >= 1 and muons_pt[0] > 30 and std::abs(muons_eta[0]) < 2.5)   leptonP4.SetPtEtaPhiM(muons_pt[0], muons_eta[0], muons_phi[0], muons_m[0]);
		else if ( electrons_n >= 1 and electrons_pt[0] and std::abs(electrons_eta[0]) < 2.5)   leptonP4.SetPtEtaPhiM(electrons_pt[0], electrons_eta[0], electrons_phi[0], electrons_m[0]);
		else continue;

		b_lepton_pt = leptonP4.Pt();
		b_lepton_eta = leptonP4.Eta();
		b_lepton_phi = leptonP4.Phi();
		b_lepton_m = leptonP4.M();


		// Information about jet
		int njets = 0;
		int nbjets = 0;
		vector<TLorentzVector> p4s; // b jets
		set<int> bjets; // not b jets

		for (int i = 0; i < jets_n; ++i)
		{
			if (jets_pt[i] < 30 or std::abs(jets_eta[i]) > 2.5) continue;
			njets++;
      TLorentzVector p4;
      p4.SetPtEtaPhiM(jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
      p4s.push_back(p4);
			if (jets_bTag[i] == 1)
			{
        bjets.insert(i);
				++nbjets;
				b_jet_pt = p4.Pt();
				b_jet_eta = p4.Eta();
				b_jet_phi = p4.Phi();
				b_jet_m = p4.M();
			}
		}
		b_njet = njets;
		b_nbjet = nbjets;


		// m3 variable 
		/*
		vector<TLorentzVector> keepm3;
		if (njets < 2)   continue;
		double m3pt = 0;
		for (int i = 2; i < njets; ++i)
		{
			for (int j = 1; j < i; ++j)
			{
				for (int k = 0; k < j; ++k)
				{
					TLorentzVector j1, j2, j3;
					j1.SetPtEtaPhiM(jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
					j2.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]); 
					j3.SetPtEtaPhiM(jets_pt[k], jets_eta[k], jets_phi[k], jets_m[k]);
					const double M3pt = (j1+j2+j3).Pt();
					if (M3pt >= m3pt)
					{
						m3pt = M3pt;
						keepm3 = { j1, j2, j3 };
					}
				}
			}
		}
		if ( keepm3.size() < 3 ) continue;
		b_m3_m = (keepm3[0]+keepm3[1]+keepm3[2]).M();
		// sorting (j1 j2 : Higgs candidate, j3 : c candidate)
		TLorentzVector temp;
		if (keepm3[0].M() > keepm3[1].M())
		{
			temp = keepm3[0];
			keepm3[0] = keepm3[1];
			keepm3[1] = temp;
		}
		if (keepm3[1].M() > keepm3[2].M())
		{
			temp = keepm3[1];
			keepm3[1] = keepm3[2];
			keepm3[2] = temp;
		}
		if (keepm3[0].M() > keepm3[1].M())
		{
			temp = keepm3[0];
			keepm3[0] = keepm3[1];
			keepm3[1] = temp;
		}
		b_m3_m12 = (keepm3[0]+keepm3[1]).M();
		b_m3_m23 = (keepm3[1]+keepm3[2]).M();
		b_m3_m31 = (keepm3[2]+keepm3[0]).M();
		b_m3_dR12 = keepm3[0].DeltaR(keepm3[1]);
		b_m3_dR23 = keepm3[1].DeltaR(keepm3[2]);
		b_m3_dR31 = keepm3[2].DeltaR(keepm3[0]);
*/

		// minimum delta R
		vector<TLorentzVector> keepdR;
		vector<TLorentzVector> keepdR3;
		int nnotb = njets-nbjets;
		if (nbjets < 2 or nnotb < 0)   continue;
		// find H > b ~b candidate
		double minDR = 999, minM = -1;
		double maxPt = 0;
		for (int i = 0; i < njets; ++i)
		{
			const auto& Bjet1 = p4s.at(i);
			for (int j = i+1; j < njets; ++j)
			{
        if ( bjets.count(i) == 0 or bjets.count(j) == 0 ) continue; // Require BOTH b jets
        //if ( bjets.count(i) == 0 and bjets.count(j) == 0 ) continue; // Require AT LEAST ONE b jet
				const auto& Bjet2 = p4s.at(j);
				const double dR = Bjet1.DeltaR(Bjet2);
				const double PT = Bjet1.Pt() >= Bjet2.Pt() ? Bjet1.Pt() : Bjet2.Pt();
				if ( dR < minDR ) 
				{
					// using Lambda
					/*
					   vector<pair<double, double> > mins;
					   std::sort(mins.begin(), mins.end(),
					   [&]( pair<double, double>a, pair<double, double> b)
					   {
					   return a.first < b.first;
					   }
					   );
					   */
					minDR = dR;
					minM = (Bjet1+Bjet2).M();
					maxPt = PT;
					keepdR = { Bjet1, Bjet2 };
				}
			}
		}
		// find c candidate
		if ( keepdR.size() < 2 ) continue;
		minDR = 999;
		const auto sumH = keepdR[0]+keepdR[1];
		for (int i = 0; i < njets; ++i)
		{
      //if ( bjets.count(i) != 0 ) continue; // Require non-b tagged jet
			const auto& Cjet = p4s.at(i);
      if ( Cjet.DeltaR(keepdR[0]) < 0.5 or Cjet.DeltaR(keepdR[1]) < 0.5 ) continue;
			const double dR = Cjet.DeltaR(sumH);
			if ( dR < minDR ) {
        keepdR3 = {keepdR[0], keepdR[1], Cjet};
        minDR = dR;
      }
		}
		if ( keepdR3.size() < 3 ) continue;
		b_dr_m = (keepdR3[0]+keepdR3[1]+keepdR3[2]).M();
		b_dr_m12 = (keepdR3[0]+keepdR3[1]).M();
		b_dr_m23 = (keepdR3[1]+keepdR3[2]).M();
		b_dr_m31 = (keepdR3[2]+keepdR3[0]).M();
		b_dr_dR12 = keepdR3[0].DeltaR(keepdR3[1]);
		b_dr_dR23 = keepdR3[1].DeltaR(keepdR3[2]);
		b_dr_dR31 = keepdR3[2].DeltaR(keepdR3[0]);

		tree->Fill();
	}
	fout->Write();
}
