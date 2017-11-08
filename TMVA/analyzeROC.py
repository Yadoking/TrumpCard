#!/usr/bin/env python

from ROOT import *
from collections import OrderedDict
import sys, os
import re

minAUC, maxAUC = 0.58, 0.68
minChi2, maxChi2 = 0.4, 60

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
    leg.SetFillStyle(0)#1111)
    return leg

gROOT.ProcessLine(".L ../plotting/tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat(".3f")
#gStyle.SetPalette(kDarkBodyRadiator)
gStyle.SetPalette(kRainBow)
gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

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

            hMVA_TestS = d.Get("MVA_%s_S" % methodName)
            hMVA_TestB = d.Get("MVA_%s_B" % methodName)
            hMVA_TrainS = d.Get("MVA_%s_Train_S" % methodName)
            hMVA_TrainB = d.Get("MVA_%s_Train_B" % methodName)

            chiS2 = hMVA_TestS.Chi2Test(hMVA_TrainS, "WW CHI2/NDF")
            chiB2 = hMVA_TestB.Chi2Test(hMVA_TrainB, "WW CHI2/NDF")

            grps[methodName] = [(fName, objPath), auc, chiS2, chiB2]

grps = OrderedDict(reversed(sorted(grps.iteritems(), key=lambda x: x[1][1])))

## Fill histograms
colors = [kRed, kBlue+1, kGreen+3, kOrange-3, kMagenta+2]
colors.extend([kRed-1, kBlue-1, kGreen-5, kOrange+3, kMagenta-1])

hAUC = TH1D("hAUC", "AUC;;Area under curve", len(grps), 0, len(grps)) ## Just for frame
grpAUC1, grpAUC2, grpAUC3 = TGraph(), TGraph(), TGraph()
grpAUCRefs = OrderedDict()
#grpAUCRefs["BDT_.*minNode2p5_nTree850"] = ["BDT >=2.5node 850tree", kAzure+1, grpAUC.Clone()]
#grpAUCRefs["DNN_ReLU_.*"] = ["DNN (ReLU)", kRed+2, TGraph()]
#grpAUCRefs["BDTG_.*"] = ["Gradient BDT", kMagenta+1, TGraph()]
grpAUCRefs["BDT_nCuts40_maxDepth3_minNode2p5_nTree850"] = ["BDT default", kBlue+1, TGraph()]
grpAUCRefs["PyKeras"] = ["PyKeras", kGreen+2, TGraph()]
#grpAUCRefs["BDTD_.*"] = ["Decorrelated BDT", kCyan+1, grpAUC.Clone()]
for name, (title, color, grp) in grpAUCRefs.iteritems():
    grp.SetMarkerStyle(kFullCircle)
    grp.SetMarkerSize(0.5)
    grp.SetMarkerColor(color)
    grp.SetLineColor(color)

grpAUCRefs["BDT_nCuts40_maxDepth3_minNode2p5_nTree850"][-1].SetMarkerSize(2)
grpAUCRefs["BDT_nCuts40_maxDepth3_minNode2p5_nTree850"][-1].SetMarkerStyle(kFullTriangleUp)
grpAUCRefs["PyKeras"][-1].SetMarkerSize(2.5)
grpAUCRefs["PyKeras"][-1].SetMarkerStyle(kFullStar)

grpChi2S1, grpChi2B1 = TGraph(), TGraph()
grpChi2S1.SetTitle("AdaBoost sig.")
grpChi2B1.SetTitle("AdaBoost bkg.")
for grp in (grpChi2S1, grpChi2B1):
    grp.SetMarkerStyle(kFullTriangleUp)
    grp.SetMarkerSize(0.7)
grpChi2S1.SetMarkerColor(kAzure+2)
grpChi2B1.SetMarkerColor(kGray+2)
grpChi2S1.SetFillColor(kAzure+2)
grpChi2B1.SetFillColor(kGray+2)

grpChi2S2, grpChi2B2 = TGraph(), TGraph()
grpChi2S2.SetTitle("BDT sig.")
grpChi2B2.SetTitle("BDT bkg.")
for grp in (grpChi2S2, grpChi2B2):
    grp.SetMarkerStyle(kFullTriangleUp)
    grp.SetMarkerSize(0.7)
grpChi2S2.SetMarkerColor(kAzure+1)
grpChi2B2.SetMarkerColor(kGray+1)
grpChi2S2.SetFillColor(kAzure+1)
grpChi2B2.SetFillColor(kGray+1)

grpChi2S3, grpChi2B3 = TGraph(), TGraph()
grpChi2S3.SetTitle("DNN sig.")
grpChi2B3.SetTitle("DNN bkg.")
grpChi2S3.SetMarkerColor(kGreen-5)
grpChi2B3.SetMarkerColor(kMagenta-8)
grpChi2S3.SetFillColor(kGreen-5)
grpChi2B3.SetFillColor(kMagenta-8)

grpChi2RefS1, grpChi2RefB1 = TGraph(), TGraph()
grpChi2RefS1.SetMarkerStyle(kFullStar)
grpChi2RefB1.SetMarkerStyle(kFullStar)
grpChi2RefS1.SetTitle("PyKeras sig.")
grpChi2RefB1.SetTitle("PyKeras bkg.")
grpChi2RefS1.SetMarkerSize(2.5)
grpChi2RefB1.SetMarkerSize(2.5)
grpChi2RefS1.SetMarkerColor(kGreen+2)
grpChi2RefB1.SetMarkerColor(kMagenta)

grpChi2RefS2, grpChi2RefB2 = TGraph(), TGraph()
grpChi2RefS2.SetTitle("BDT sig. (opt.)")
grpChi2RefB2.SetTitle("BDT bkg. (opt.)")
grpChi2RefS2.SetMarkerSize(2.0)
grpChi2RefB2.SetMarkerSize(2.0)
grpChi2RefS2.SetMarkerStyle(kFullTriangleUp)
grpChi2RefB2.SetMarkerStyle(kFullTriangleUp)
grpChi2RefS2.SetMarkerColor(kBlue+1)
grpChi2RefB2.SetMarkerColor(kRed)

dnnNodes = [16,32,64,128,256,512]
nX, nY = len(dnnNodes), 20
hAUC2D_DNN   = TH2D("hAUC2D_DNN"  , "hAUC2D_DNN;Number of nodes;Number of Layers", nX, 0, nX, nY, 1, nY+1)
hAUC2D_Keras = TH2D("hAUC2D_Keras", "hAUC2D_Keras;Number of nodes;Number of Layers", nX, 0, nX, nY, 1, nY+1)
for i, x in enumerate(dnnNodes):
    hAUC2D_DNN.GetXaxis().SetBinLabel(i+1, "%d" % x)
    hAUC2D_Keras.GetXaxis().SetBinLabel(i+1, "%d" % x)

for i, (name, [(fName, objPath), auc, chiS2, chiB2]) in enumerate(grps.iteritems()):
    if name.startswith("BDT_"):
        grpAUC1.SetPoint(grpAUC1.GetN(), i, auc)
        grpChi2S1.SetPoint(i, auc, chiS2)
        grpChi2B1.SetPoint(i, auc, chiB2)
    elif 'BDT' in name:
        grpAUC2.SetPoint(grpAUC2.GetN(), i, auc)
        grpChi2S2.SetPoint(i, auc, chiS2)
        grpChi2B2.SetPoint(i, auc, chiB2)
    else:
        grpAUC3.SetPoint(grpAUC3.GetN(), i, auc)
        grpChi2S3.SetPoint(i, auc, chiS2)
        grpChi2B3.SetPoint(i, auc, chiB2)
    for pattern0, (title, color, grp) in grpAUCRefs.iteritems():
        pattern = pattern0.replace("+", "\\+")
        if not re.match('^'+pattern+'$', name): continue
        grp.SetPoint(grp.GetN(), i, auc)
    if name == "PyKeras":
        grpChi2RefS1.SetPoint(grpChi2RefS1.GetN(), auc, chiS2)
        grpChi2RefB1.SetPoint(grpChi2RefB1.GetN(), auc, chiB2)
    if name == "BDT_nCuts40_maxDepth3_minNode2p5_nTree850":
        grpChi2RefS2.SetPoint(grpChi2RefS2.GetN(), auc, chiS2)
        grpChi2RefB2.SetPoint(grpChi2RefB2.GetN(), auc, chiB2)
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
for i, (name, [(fName, objPath), auc, chiS2, chiB2]) in enumerate(grps.iteritems()):
    if len(grps) > 10 and i >= 5: break
    grp = buildROC(fName, objPath)
    grp.SetLineColor(colors[i])
    grpsToDraw.append(grp)
    legROC.AddEntry(grp, name, "l")
if len(grps) > 10:
    legROC.AddEntry(0, "", "p").SetMarkerColor(kWhite)
    legROC.AddEntry(0, "Worst 5", "p").SetMarkerColor(kWhite)
    for i, (name, [(fName, objPath), auc, chiS2, chiB2]) in enumerate(reversed([x for x in grps.iteritems()])):
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
hAUC.SetMinimum(minAUC)
hAUC.SetMaximum(maxAUC)
hAUC.GetXaxis().SetNdivisions(505)
hAUC.Draw()
legAUC = buildLegend(0.2, 0.2, 0.6, 0.4)
grpAUC1.SetMarkerColor(kGray+2)
grpAUC2.SetMarkerColor(kGray+1)
grpAUC3.SetMarkerColor(kMagenta-8)
for grp in (grpAUC1, grpAUC2, grpAUC3):
    grp.SetMarkerStyle(kFullCircle)
    grp.SetMarkerSize(0.5)
    grp.Draw("p")
    grp.SetEditable(False)
legAUC.AddEntry(grpAUC1, "AdaBoost", "p")
legAUC.AddEntry(grpAUC2, "BDTG/BaggedBDT", "p")
legAUC.AddEntry(grpAUC3, "All DNN", "p")
for name, (title, color, grp) in grpAUCRefs.iteritems():
    legAUC.AddEntry(grp, "%s" % title, "p")
    grp.Draw("p")
    grp.SetEditable(False)
legAUC.Draw()

cChi2 = TCanvas("cChi2", "cChi2", 500, 500)
cChi2.SetLogy()
hFrameChi2 = TH1F("hChi2Frame", "hChi2Frame;Area under curve;#chi^{2}/dof", 100, minAUC, maxAUC)
hFrameChi2.SetMinimum(minChi2)
hFrameChi2.SetMaximum(maxChi2)
hFrameChi2.GetXaxis().SetNdivisions(505)
hFrameChi2.Draw()
legChi2 = buildLegend(0.45, 0.70, 0.87, 0.87)
legChi2.SetNColumns(2)
for grp in (grpChi2B2, grpChi2S2, grpChi2B1, grpChi2S1, grpChi2B3, grpChi2S3):
    grp.Draw("p")
    grp.SetEditable(False)
    legChi2.AddEntry(grp, grp.GetTitle(), "f")
for grp in (grpChi2RefB1, grpChi2RefS1, grpChi2RefB2, grpChi2RefS2):
    grp.Draw("p")
    grp.SetEditable(False)
    legChi2.AddEntry(grp, grp.GetTitle(), "p")
legChi2.Draw()

cAUC2D_DNN = TCanvas("cAUC2D_DNN", "cAUC2D_DNN", 500, 500)
cAUC2D_DNN.SetRightMargin(0.16)
cAUC2D_DNN.SetLeftMargin(0.14)
hAUC2D_DNN.GetYaxis().SetTitleOffset(1)
hAUC2D_DNN.SetMinimum(minAUC)
hAUC2D_DNN.SetMaximum(maxAUC)
hAUC2D_DNN.Draw("COLZTEXT")

cAUC2D_Keras = TCanvas("cAUC2D_Keras", "cAUC2D_Keras", 500, 500)
cAUC2D_Keras.SetRightMargin(0.16)
cAUC2D_Keras.SetLeftMargin(0.14)
hAUC2D_Keras.GetYaxis().SetTitleOffset(1)
hAUC2D_Keras.SetMinimum(minAUC)
hAUC2D_Keras.SetMaximum(maxAUC)
hAUC2D_Keras.Draw("COLZTEXT")

for x in grps.keys()[:5]:
    print grps[x][0][0], grps[x][0][1].split('/')[-2], grps[x][1]
