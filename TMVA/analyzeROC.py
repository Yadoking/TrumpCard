#!/usr/bin/env python

from ROOT import *
from collections import OrderedDict
import sys, os

if len(sys.argv) < 2:
    print "Usage: %s INPUTFILE.root" % sys.argv[0]
    sys.exit(1)

fName = sys.argv[1]
if not os.path.exists(fName):
    print "Cannot find input root file %s" % fName
    sys.exit(2)

gROOT.ProcessLine(".L ../FlatTuple/tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

grps = OrderedDict()
colors = [kBlack, kBlue, kRed, kMagenta, kGreen+1]

try:
    f = TFile(fName)
except(e):
    print "Input file %f is not valid root file" % fName
    sys.exit(2)

for dName in [key.GetName() for key in f.GetListOfKeys()]:
    dBase = f.Get(dName)
    for dName in [key.GetName() for key in dBase.GetListOfKeys()]:
        if dName.split('_', 1)[0] != 'Method': continue
        d = dBase.GetDirectory(dName)
        if d == None: continue

        methodName = dName.split('_',1)[1]
        hSigEff = d.Get("%s/MVA_%s_effS" % (methodName, methodName))
        hBkgEff = d.Get("%s/MVA_%s_effB" % (methodName, methodName))

        grp = TGraph()
        grp.SetLineColor(colors[len(grps)-1])
        for b in range(hSigEff.GetNbinsX()):
            mvaVal = hSigEff.GetXaxis().GetBinLowEdge(b+1)
            x = hSigEff.GetBinContent(b+1)
            y = 1-hBkgEff.GetBinContent(b+1)
            np = grp.GetN()
            if np > 0:
                prevX = grp.GetX()[np-1]
                prevY = grp.GetY()[np-1]
                if x == prevX or y == prevY: continue
            grp.SetPoint(np, x, y)
        auc = 0.0
        for p in range(grp.GetN()-1):
            x1, x2 = grp.GetX()[p], grp.GetX()[p+1]
            y1, y2 = grp.GetY()[p], grp.GetY()[p+1]
            dx = abs(x2-x1)
            y = (y1+y2)/2
            if dx == 0: continue
            auc += y*dx

        grps[methodName] = (grp, auc)
grps = OrderedDict(reversed(sorted(grps.iteritems(), key=lambda x: x[1][1])))

cROC = TCanvas("cROC", "cROC", 500, 500)
hFrameROC = TH2F("hFrameROC", "hFrame;Signal efficiency;Background rejection", 100, 0, 1, 100, 0, 1)
hFrameROC.Draw()
leg = TLegend(0.2, 0.2, 0.5, 0.5)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
cAUC = TCanvas("cAUC", "cAUC", 500, 500)
hAUC = TH1D("hAUC", "AUC;;Area under curve", len(grps), 0, len(grps))
hAUC.SetMinimum(0.5)
for i, (name, (grp, auc)) in enumerate(grps.iteritems()):
    hAUC.GetXaxis().SetBinLabel(i+1, name)
    hAUC.SetBinContent(i+1, auc)

    cROC.cd()
    grp.Draw("l")

    leg.AddEntry(grp, name, "l")
cROC.cd()
leg.Draw()

cAUC.cd()
hAUC.Draw()
