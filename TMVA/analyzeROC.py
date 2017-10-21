#!/usr/bin/env python

from ROOT import *
from collections import OrderedDict
import sys, os
import re

if len(sys.argv) < 2:
    print "Usage: %s INPUTFILE.root" % sys.argv[0]
    sys.exit(1)

def buildROC(fName, objPath):
    f = TFile(fName)
    h = f.Get(objPath)

    grp = TGraph()
    for i in range(h.GetNbinsX()):
      x = h.GetXaxis().GetBinCenter(i+1)
      y = h.GetBinContent(i+1)
      grp.SetPoint(i, x, y)
    grp.SetLineWidth(2)
    return grp

def buildLegend(x1=0.6, y1=0.75, x2=0.85, y2=0.87):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1111)
    return leg

gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat(".3f")
#gStyle.SetPalette(kDarkBodyRadiator)
gStyle.SetPalette(kRainBow)

grps = OrderedDict()

## Load ROC histograms and compute AUC
for fName in sys.argv[1:]:
    f = TFile(fName)
    if f == None: continue

    for dBaseName0 in [key.GetName() for key in f.GetListOfKeys()]:
        dBase = f.Get(dBaseName0)
        for dName in [key.GetName() for key in dBase.GetListOfKeys()]:
            if dName.split('_', 1)[0] != 'Method': continue
            d = dBase.GetDirectory(dName)
            if d == None: continue

            methodName = dName.split('_', 1)[-1]
            d = d.GetDirectory(methodName)

            hROC = d.Get("MVA_%s_rejBvsS" % methodName)
            objPath = d.GetPath().split(':',1)[-1]+'/'+hROC.GetName()
            auc = hROC.Integral()/hROC.GetNbinsX()

            grps[methodName] = [(fName, objPath), auc]
grps = OrderedDict(reversed(sorted(grps.iteritems(), key=lambda x: x[1][1])))

## Fill histograms
colors = [kRed, kBlue+1, kGreen+3, kOrange-3, kMagenta+2]
colors.extend([kRed-1, kBlue-1, kGreen-5, kOrange+3, kMagenta-1])

hAUC = TH1D("hAUC", "AUC;;Area under curve", len(grps), 1, len(grps)+1)
hAUCRefs = OrderedDict()
#hAUCRefs["BDT_GiniIndex_.*minNode2p5_nTree850"] = ["BDT >=2.5node 850tree", kAzure+1, hAUC.Clone()]
#hAUCRefs["BDT_GiniIndex_.*nCuts80.*"] = ["BDT >=2.5node 850tree", kAzure+1, hAUC.Clone()]
hAUCRefs["DNN_ReLU_.*"] = ["DNN (ReLU)", kRed+2, hAUC.Clone()]
hAUCRefs["BaggedBDT_GiniIndex_nCuts20_maxDepth3_minNode2p5_nTree850"] = ["BDT default", kBlue+1, hAUC.Clone()]

dnnNodes = [16,32,64,128,256,512]
nX, nY = len(dnnNodes), 20
hAUC2D_DNN   = TH2D("hAUC2D_DNN"  , "hAUC2D_DNN;Number of nodes;Number of Layers", nX, 0, nX, nY, 1, nY+1)
hAUC2D_Keras = TH2D("hAUC2D_Keras", "hAUC2D_Keras;Number of nodes;Number of Layers", nX, 0, nX, nY, 1, nY+1)
for i, x in enumerate(dnnNodes):
    hAUC2D_DNN.GetXaxis().SetBinLabel(i+1, "%d" % x)
    hAUC2D_Keras.GetXaxis().SetBinLabel(i+1, "%d" % x)

for i, (name, [(fName, objPath), auc]) in enumerate(grps.iteritems()):
    hAUC.SetBinContent(i+1, auc)
    for pattern, (title, color, h) in hAUCRefs.iteritems():
        if not re.match('^'+pattern+'$', name): continue
        h.SetBinContent(i+1, auc)
    if 'DNN_' in name and '_X' in name:
        w, x, y = name.split('DNN_')[-1].split('_')
        x, y = int(x[1:]), int(y[1:])
        xbin = hAUC2D_DNN.GetXaxis().FindBin("%d" % x)
        ybin = hAUC2D_DNN.GetYaxis().FindBin(y)
        if xbin < 1 or xbin > hAUC2D_DNN.GetNbinsX() or \
           ybin < 1 or ybin > hAUC2D_DNN.GetNbinsY(): continue

        hAUC2D_DNN.SetBinContent(xbin, ybin, auc)
    if 'Keras_' in name and '_X' in name:
        w, x, y = name.split('Keras_')[-1].split('_')
        x, y = int(x[1:]), int(y[1:])
        xbin = hAUC2D_Keras.GetXaxis().FindBin("%d" % x)
        ybin = hAUC2D_Keras.GetYaxis().FindBin(y)
        if xbin < 1 or xbin > hAUC2D_Keras.GetNbinsX() or \
           ybin < 1 or ybin > hAUC2D_Keras.GetNbinsY(): continue

        hAUC2D_Keras.SetBinContent(xbin, ybin, auc)

grpsToDraw = []
legROC = buildLegend(0.2, 0.15, 0.8, 0.6)
if len(grps) > 5:
    legROC.AddEntry(0, "Top 5", "p").SetMarkerColor(kWhite)
for i, (name, [(fName, objPath), auc]) in enumerate(grps.iteritems()):
    if len(grps) > 10 and i >= 5: break
    grp = buildROC(fName, objPath)
    grp.SetLineColor(colors[i])
    grpsToDraw.append(grp)
    legROC.AddEntry(grp, name, "l")
if len(grps) > 10:
    legROC.AddEntry(0, "", "p").SetMarkerColor(kWhite)
    legROC.AddEntry(0, "Worst 5", "p").SetMarkerColor(kWhite)
    for i, (name, [(fName, objPath), auc]) in enumerate(reversed([x for x in grps.iteritems()])):
        if i >= 5: break
        grp = buildROC(fName, objPath)
        grp.SetLineColor(colors[i+5])
        grpsToDraw.append(grp)
        legROC.AddEntry(grp, name, "l")

## Start drawing
cROC = TCanvas("cROC", "cROC", 500, 500)
hFrameROC = TH2F("hFrameROC", "hFrame;Signal efficiency;Background rejection", 100, 0, 1, 100, 0, 1)
hFrameROC.Draw()
legROC.Draw()
for grp in grpsToDraw: grp.Draw("l")
cROC.SetGridx()
cROC.SetGridy()

cAUC = TCanvas("cAUC", "cAUC", 500, 500)
#hAUC.SetMinimum(hAUC.GetMinimum()*0.9)
#hAUC.SetMaximum(hAUC.GetMaximum()*1.1)
hAUC.Draw()
legAUC = buildLegend()
for name, (title, color, h) in hAUCRefs.iteritems():
    legAUC.AddEntry(h, "%s" % title, "f")
    h.SetFillColor(color)
    h.SetLineColor(color)
    h.Draw("same")
legAUC.Draw()

cAUC2D_DNN = TCanvas("cAUC2D_DNN", "cAUC2D_DNN", 500, 500)
cAUC2D_DNN.SetRightMargin(0.16)
cAUC2D_DNN.SetLeftMargin(0.14)
hAUC2D_DNN.GetYaxis().SetTitleOffset(1)
hAUC2D_DNN.SetMinimum(0.60)
hAUC2D_DNN.SetMaximum(0.75)
hAUC2D_DNN.Draw("COLZTEXT")

cAUC2D_Keras = TCanvas("cAUC2D_Keras", "cAUC2D_Keras", 500, 500)
cAUC2D_Keras.SetRightMargin(0.16)
cAUC2D_Keras.SetLeftMargin(0.14)
hAUC2D_Keras.GetYaxis().SetTitleOffset(1)
hAUC2D_Keras.SetMinimum(0.60)
hAUC2D_Keras.SetMaximum(0.75)
hAUC2D_Keras.Draw("COLZTEXT")

for x in grps.keys()[:5]:
    print grps[x][0][0], grps[x][0][1].split('/')[-2], grps[x][1]
