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
colors = [kRed, kBlue+1, kGreen+3, kOrange-3, kMagenta+2]
colors.extend([kRed-1, kBlue-1, kGreen-5, kOrange+3, kMagenta-1])

files = []
for fName in sys.argv[1:]:
    f = TFile(fName)
    if f == None: continue
    files.append(f)

    for dBaseName0 in [key.GetName() for key in f.GetListOfKeys()]:
        dBase = f.Get(dBaseName0)
        for dName in [key.GetName() for key in dBase.GetListOfKeys()]:
            if dName.split('_', 1)[0] != 'Method': continue
            d = dBase.GetDirectory(dName)
            if d == None: continue

            methodName = dName.split('_',1)[1]
            hSigEff = d.Get("%s/MVA_%s_effS" % (methodName, methodName))
            hBkgEff = d.Get("%s/MVA_%s_effB" % (methodName, methodName))

            grp = TGraph()
            grp.SetLineWidth(2)
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

            dBaseName = dBaseName0.replace("mva_delphes", "").replace("mva_cmsTuple", "").strip("_")
            grps["%s %s" % (dBaseName, methodName)] = (grp, auc)
grps = OrderedDict(reversed(sorted(grps.iteritems(), key=lambda x: x[1][1])))

cROC = TCanvas("cROC", "cROC", 500, 500)
hFrameROC = TH2F("hFrameROC", "hFrame;Signal efficiency;Background rejection", 100, 0, 1, 100, 0, 1)
hFrameROC.Draw()
leg = TLegend(0.2, 0.15, 0.6, 0.6)
leg.SetBorderSize(0)
leg.SetFillStyle(1111)
cAUC = TCanvas("cAUC", "cAUC", 500, 500)
hAUC = TH1D("hAUC", "AUC;;Area under curve", len(grps), 0, len(grps))
hAUC.SetMinimum(0.5)
cAUC2D_DNN = TCanvas("cAUC2D_DNN", "cAUC2D_DNN", 500, 500)
cAUC2D_Keras = TCanvas("cAUC2D_Keras", "cAUC2D_Keras", 500, 500)
dnnFtns, dnnNodes = ["TANH", "ReLU"], [16,32,64,128,256,512]
hAUC2D_DNN= TH2D("hAUC2D_DNN", "hAUC2D_DNN;Number of nodes;Number of Layers", len(dnnFtns)*len(dnnNodes), 0, len(dnnFtns)*len(dnnNodes), 11, 1, 11)
hAUC2D_Keras= TH2D("hAUC2D_Keras", "hAUC2D_Keras;Number of nodes;Number of Layers", len(dnnFtns)*len(dnnNodes), 0, len(dnnFtns)*len(dnnNodes), 11, 1, 11)
for iw, w in enumerate(dnnFtns):
    for jx, x in enumerate(dnnNodes):
        hAUC2D_DNN.GetXaxis().SetBinLabel(iw*len(dnnNodes)+jx+1, "%s_%d" % (w, x))
        hAUC2D_Keras.GetXaxis().SetBinLabel(iw*len(dnnNodes)+jx+1, "%s_%d" % (w, x))
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

    if 'DNN_' in name and '_X' in name:
        w, x, y = name.split('DNN_')[-1].split('_')
        x, y = int(x[1:]), int(y[1:])
        w = dnnFtns.index(w)
        xbin = w*len(dnnNodes)+dnnNodes.index(x)+1
        ybin = range(1,11).index(y)+1

        hAUC2D_DNN.SetBinContent(xbin, ybin, auc)

    if 'Keras_' in name and '_X' in name:
        w, x, y = name.split('Keras_')[-1].split('_')
        x, y = int(x[1:]), int(y[1:])
        w = dnnFtns.index(w)
        xbin = w*len(dnnNodes)+dnnNodes.index(x)+1
        ybin = range(1,11).index(y)+1

        hAUC2D_Keras.SetBinContent(xbin, ybin, auc)

    hBkgAt90.GetXaxis().SetBinLabel(i+1, name)
    hBkgAt90.SetBinContent(i+1, grp.Eval(0.9))
    hBkgAt80.SetBinContent(i+1, grp.Eval(0.8))

if len(grps) > 5:
    leg.AddEntry(0, "Top 5", "p").SetMarkerColor(kWhite)
grpsToDraw = []
for i, (name, (grp, auc)) in enumerate(grps.iteritems()):
    grp.SetLineColor(colors[i])
    grpsToDraw.append(grp)
    leg.AddEntry(grp, name, "l")
    if len(grps) > 10 and i >= 4: break

if len(grps) > 10:
    leg.AddEntry(0, "", "p").SetMarkerColor(kWhite)
    leg.AddEntry(0, "Worst 5", "p").SetMarkerColor(kWhite)

    for i, (name, (grp, auc)) in enumerate(reversed([x for x in grps.iteritems()])):
        grp.SetLineColor(colors[i+5])
        grpsToDraw.append(grp)
        leg.AddEntry(grp, name, "l")
        if i >= 4: break
cROC.cd()
for grp in grpsToDraw: grp.Draw("l")
cROC.SetGridx()
cROC.SetGridy()
leg.Draw()

cAUC.cd()
hAUC.Draw()

cAUC2D_DNN.cd()
hAUC2D_DNN.Draw("COLZTEXT")

cAUC2D_Keras.cd()
hAUC2D_Keras.Draw("COLZTEXT")

cBkgAtWP.cd()
hBkgAt90.Draw()
hBkgAt80.Draw("same")
legBkgAtWP.Draw()

for c in [cROC, cAUC, cBkgAtWP]:
    c.Modified()
    c.Update()
