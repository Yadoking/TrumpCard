#!/usr/bin/env python

from ROOT import *
from PlotMaker import *

basedir = "../FlatTuple/full_shift_rotation_flip/m3"

h = HistInfo()
h.setLumi(1.0)

xsec_ttbar = 356.4
xsec_ttbb  = 356.4
br_ttLJ = 0.40
br_tcH = 1.00

## Set data
#h.addData("ntuple/DataSingleEG.root")

## Set MC samples
h.addSig("tcH", "t#bar{t}, t#rightarrow{}Hc (BR=%f%%)" % br_tcH*100, "%s/delphes_FCNC.root" % basedir, kBlue, xsec_ttbar*br_tcH*br_ttLJ, 500000)

h.addBkg("ttbb", "t#bar{t}+b#bar{b}", "%s/delphes_ttbb.root" % basedir, kRed, xsec_ttbb*br_ttLJ, 500000)

## Define plots
h.add1D("met_pt", "met_pt", "MET p_{T} (GeV);Events / 2GeV", 100, 0, 200)
h.add1D("jets_n", "jets_n", "Jet multiplicity;Events", 10, 0, 10)
h.add1D("bjets_n", "bjets_n", "B-jet multiplicity;Events", 10, 0, 10)
h.add1D("hadT_m", "kin_hadT_m", "Hadronic top mass (GeV);Events / 5GeV", 100, 0, 500)
h.add1D("hadW12_m", "kin_hadW12_m", "Hadronic W_{12} mass (GeV);Events / 3GeV", 100, 0, 300)

## Define cut flows
h.addCutStep("step1", "lepton_pt>30 && abs(lepton_eta)<2.1", "lepton_pt,lepton_eta,jets_n,met_pt", "weight_gen")
h.addCutStep("step2", "met_pt >= 30", "jets_n,bjets_n,met_pt", "weight_gen")
h.addCutStep("step3", "jets_n >= 4", "jets_n,bjets_n,met_pt", "weight_gen")
h.addCutStep("step4", "bjets_n >= 1", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m", "weight_gen")
h.addCutStep("step4", "bjets_n >= 2", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m", "weight_gen")
h.addCutStep("step4", "bjets_n >= 3", "jets_n,bjets_n,met_pt,hadT_m,hadW12_m", "weight_gen")

## Produce histograms
if not os.path.exists("hist.root"):
    hm = HistMaker("tree", h)
    hm.applyCutSteps("hist.root")

## Produce plots
pm = PlotMaker("", "hist.root", h, {'width':400, 'aspectRatio':1., 'doRatio':False})

