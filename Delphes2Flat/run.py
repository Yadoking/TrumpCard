#!/usr/bin/env python

from ROOT import *
import os

gSystem.Load("libDelphes")
gROOT.ProcessLine(".L makeFlatTuple.C++")

if not os.path.exists("ntuple_tch.root"):
    print "@@ Processing tch..."
    gROOT.ProcessLine('makeFlatTuple("/home/minerva1993/phenostudy/delphes_analysis/rootfiles/tch_v2.root", "ntuple_tch.root");')

if not os.path.exists("ntuple_ttbb.root"):
    print "@@ Processing ttbb..."
    gROOT.ProcessLine('makeFlatTuple("/home/minerva1993/phenostudy/delphes_analysis/rootfiles/ttbb_v2.root", "ntuple_ttbb.root");')
