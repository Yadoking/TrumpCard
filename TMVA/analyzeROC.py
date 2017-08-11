#!/usr/bin/env python

from ROOT import *
from collections import OrderedDict
import sys, os

if len(sys.argv) < 2:
    print "Usage: %s INPUTFILE.root" % sys.argv[0]
    sys.exit(1)

gROOT.ProcessLine(".L ../FlatTuple/tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

grps = OrderedDict()
colors = [kBlack, kBlue, kRed, kMagenta, kGreen+1, kYellow+2, kAzure+2, kRed+4]

files = []
for fName in sys.argv[1:]:
    f = TFile(fName)
    if f == None: continue
    files.append(f)

    for dBaseName in [key.GetName() for key in f.GetListOfKeys()]:
        dBase = f.Get(dBaseName)
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

            grps["%s %s" % (dBaseName, methodName)] = (grp, auc)
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
cBkgAtWP = TCanvas("cBkgAtWP", "cBkgAtWP", 500, 500)
hBkgAt90 = TH1D("hBkgAt90", "90% WP;;Background rejection", len(grps), 0, len(grps))
hBkgAt90.SetMinimum(0)
hBkgAt90.SetMaximum(1)
hBkgAt80 = TH1D("hBkgAt80", "80% WP;;Background rejection", len(grps), 0, len(grps))
hBkgAt90.SetLineColor(kBlue)
hBkgAt80.SetLineColor(kRed)
legBkgAtWP = TLegend(0.65, 0.75, 0.85, 0.85)
legBkgAtWP.SetBorderSize(0)
legBkgAtWP.SetFillStyle(0)
legBkgAtWP.AddEntry(hBkgAt80, "80% WP", "l")
legBkgAtWP.AddEntry(hBkgAt90, "90% WP", "l")
for i, (name, (grp, auc)) in enumerate(grps.iteritems()):
    hAUC.GetXaxis().SetBinLabel(i+1, name)
    hAUC.SetBinContent(i+1, auc)

    hBkgAt90.GetXaxis().SetBinLabel(i+1, name)
    hBkgAt90.SetBinContent(i+1, grp.Eval(0.9))
    hBkgAt80.SetBinContent(i+1, grp.Eval(0.8))

    cROC.cd()
    grp.Draw("l")

    leg.AddEntry(grp, name, "l")
cROC.cd()
leg.Draw()

cAUC.cd()
hAUC.Draw()

cBkgAtWP.cd()
hBkgAt90.Draw()
hBkgAt80.Draw("same")
legBkgAtWP.Draw()
