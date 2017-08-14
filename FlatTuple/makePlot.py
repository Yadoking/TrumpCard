#!/usr/bin/env python

from ROOT import *
from PlotMaker import *

#basedir = "kin"
#basedir = "m3"
basedir = "deltaR"

h = HistInfo()
h.setLumi(1.0)

xsec_ttbar  = 356.4/3
xsec_ttbb = 13.93
br_tcH = 0.01
xsec_tcH= xsec_ttbar*br_tcH

## Set data
#h.addData("ntuple/DataSingleEG.root")

## Set MC samples
h.addSig("FCNC", "FCNC Br=%g%%" % (br_tcH*100), "%s/delphes_FCNC.root" % basedir, kBlue, xsec_tcH, 500000)
h.addSig("ttbb", "t#bar{t}+b#bar{b}", "%s/delphes_ttbb.root" % basedir, kRed+3, xsec_ttbb, 500000)
h.addBkg("tt", "t#bar{t}", "%s/delphes_tt.root" % basedir, kRed, xsec_ttbar, 1000000)

## Define plots
h.add1D("met_pt", "met_pt", "MET p_{T} (GeV);Events / 5GeV", 40, 0, 200)
h.add1D("jets_n", "jets_n", "Jet multiplicity;Events", 10, 0, 10)
h.add1D("bjets_n", "bjets_n", "B-jet multiplicity;Events", 10, 0, 10)
h.add1D("hadT_m", "hadT_m", "Hadronic top mass (GeV);Events / 10GeV", 50, 0, 500)
h.add1D("hadW12_m", "hadW12_m", "Hadronic W_{12} mass (GeV);Events / 10GeV", 50, 0, 500)
h.add1D("theta1", "theta1", "#theta_{1};Events / 0.1", 64, -3.2/2, 3.2/2)
h.add1D("theta2", "theta2", "#theta_{2};Events / 0.1", 64, -3.2/2, 3.2/2)

## Define cut flows
h.addCutStep("step0", "bjetcode%10>=2", "", "weight_gen")
h.addCutStep("step1", "lepton_pt>30 && abs(lepton_eta)<2.1", "lepton_pt,lepton_eta,jets_n,met_pt", "weight_gen")
h.addCutStep("step2", "met_pt >= 30", "jets_n,bjets_n,met_pt", "weight_gen")
h.addCutStep("step3", "jets_n >= 4", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m,theta1,theta2", "weight_gen")
h.addCutStep("step4", "bjets_n >= 1", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m,theta1,theta2", "weight_gen")
h.addCutStep("step5", "bjets_n >= 2", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m,theta1,theta2", "weight_gen")
h.addCutStep("step6", "bjets_n >= 3", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m,theta1,theta2", "weight_gen")

## Produce histograms
if not os.path.exists("hist.root"):
    hm = HistMaker(h, "tree", True)
    hm.applyCutSteps("hist.root")

## Produce plots
pm = PlotMaker(h, "", "hist.root", {'width':400, 'aspectRatio':1., 'doRatio':False})

