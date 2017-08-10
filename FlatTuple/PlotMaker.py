#!/usr/bin/env python

import sys, os
from collections import OrderedDict
from ROOT import *
from array import array

class HistInfo:
    def __init__(self, *args):
        self.lumi = 1.0
        self.cutsteps = OrderedDict()
        self.plots = {}

        self.samples_RD = {'fNames':[]}
        self.samples_sig = OrderedDict()
        self.samples_bkg = OrderedDict()

        if len(args) > 0 and type(args[0]) == type(self):
            src = args[0]
            self.lumi = src.lumi
            self.cutsteps = src.cutsteps
            self.plots = src.plots

            self.samples_RD = src.samples_RD
            self.samples_sig = src.samples_sig
            self.samples_bkg = src.samples_bkg

    def setLumi(self, lumi):
        self.lumi = lumi

    def addCutStep(self, name, cut, plots, weight):
        if type(plots) == str: plots = [x.strip() for x in plots.split(',') if x.strip() != '']
        self.cutsteps[name] = (cut, plots, weight)

    def add1D(self, name, varExpr, title, *arg):
        if len(arg) < 1:
            print "!!! Missing arguments in add1D(%s, %s). Please give binning information." % (name, varExpr)
            return
        if type(arg[0]) == list: hArgs = (len(arg[0])-1, array('d', arg[0]))
        elif type(arg[0]) == array: hArgs = (len(arg[0]-1, arg[0]))
        elif len(arg) == 3 and type(arg[0]) == int: hArgs = arg[:]
        else:
            print "!!! Wrong arguments in add1D(%s, %s, args)" % (name, varExpr), "with args=", arg
            return
        
        self.plots[name] = (varExpr, title, hArgs)

    def addData(self, fNames, doReset=False):
        if type(fNames) == str: fNames = fNames.split(',')
        if doReset: self.samples_RD['fNames'] = fNames
        else: self.samples_RD['fNames'].extend(fNames)

    def addSig(self, name, title, fNames, color, xsec, evt, doReset=False):
        if type(fNames) == str: fNames = fNames.split(',')
        nEvent = 0.0
        if type(evt) == int or type(evt) == float:
            nEvent = evt
        elif type(evt) == str:
            for fName in fNames:
                f = TFile(fName)
                nEvent += f.Get(evt).GetBinContent(2)

        if title not in self.samples_sig:
            self.samples_sig[title] = {'color':color, 'subsamples':OrderedDict()}
        if name not in self.samples_sig[title]['subsamples']:
            self.samples_sig[title]['subsamples'][name] = {'fNames':[], 'xsec':xsec, 'evt':0.0}

        if doReset:
            self.samples_sig[title]['subsamples'][name]['fNames'] = fNames
            self.samples_sig[title]['subsamples'][name]['evt'] = nEvent
        else:
            self.samples_sig[title]['subsamples'][name]['fNames'].extend(fNames)
            self.samples_sig[title]['subsamples'][name]['evt'] += nEvent

    def addBkg(self, name, title, fNames, color, xsec, evt, doReset=False):
        if type(fNames) == str: fNames = fNames.split(',')
        nEvent = 0.0
        if type(evt) == int or type(evt) == float:
            nEvent = evt
        elif type(evt) == str:
            for fName in fNames:
                f = TFile(fName)
                nEvent += f.Get(evt).GetBinContent(2)

        if title not in self.samples_bkg:
            self.samples_bkg[title] = {'color':color, 'subsamples':OrderedDict()}
        if name not in self.samples_bkg[title]['subsamples']:
            self.samples_bkg[title]['subsamples'][name] = {'fNames':[], 'xsec':xsec, 'evt':0.0}

        if doReset:
            self.samples_bkg[title]['subsamples'][name]['fNames'] = fNames
            self.samples_bkg[title]['subsamples'][name]['evt'] = nEvent
        else:
            self.samples_bkg[title]['subsamples'][name]['fNames'].extend(fNames)
            self.samples_bkg[title]['subsamples'][name]['evt'] += nEvent

class HistMaker(HistInfo):
    def __init__(self, treeName, doProof, *args):
        HistInfo.__init__(self, *args)
        self.treeName = treeName
        self.prf = None
        if doProof: self.prf = TProof.Open("workers=16")

    def applyCutSteps(self, foutName):
        isBatch = gROOT.IsBatch()
        gROOT.SetBatch(True)

        fout = TFile(foutName, "recreate")
        fout.cd()

        hCutFlows = {}
        hCutFlowsRaw = {}
        if len(self.samples_RD['fNames']) > 0:
            self.samples_RD['chain'] = TChain(self.treeName)
            if self.prf != None: self.samples_RD['chain'].SetProof()
            for fName in self.samples_RD['fNames']:
                self.samples_RD['chain'].Add(fName)
            hCutFlows["RD"] = TH1D("hCutFlows_%s" % "RD", "Cut flow %s;;Events" % "RD", len(self.cutsteps), 0., len(self.cutsteps))
            hCutFlowsRaw["RD"] = TH1D("hCutFlowsRaw_%s" % "RD", "Raw Cut flow %s;;Events" % "RD", len(self.cutsteps), 0., len(self.cutsteps))
        for sInfo in self.samples_sig.values():
            for ssName, ssInfo in sInfo['subsamples'].iteritems():
                hCutFlows[ssName] = TH1D("hCutFlows_%s" % ssName, "Cut flow %s;;Events" % ssName, len(self.cutsteps), 0., len(self.cutsteps))
                hCutFlowsRaw[ssName] = TH1D("hCutFlowsRaw_%s" % ssName, "Raw Cut flow %s;;Events" % ssName, len(self.cutsteps), 0., len(self.cutsteps))
                ssInfo['chain'] = TChain(self.treeName)
                if self.prf != None: ssInfo['chain'].SetProof()
                for fName in ssInfo['fNames']: ssInfo['chain'].Add(fName)
        for sInfo in self.samples_bkg.values():
            for ssName, ssInfo in sInfo['subsamples'].iteritems():
                hCutFlows[ssName] = TH1D("hCutFlows_%s" % ssName, "Cut flow %s;;Events" % ssName, len(self.cutsteps), 0., len(self.cutsteps))
                hCutFlowsRaw[ssName] = TH1D("hCutFlowsRaw_%s" % ssName, "Raw Cut flow %s;;Events" % ssName, len(self.cutsteps), 0., len(self.cutsteps))
                ssInfo['chain'] = TChain(self.treeName)
                if self.prf != None: ssInfo['chain'].SetProof()
                for fName in ssInfo['fNames']: ssInfo['chain'].Add(fName)

        ## Data
        if len(self.samples_RD['fNames']) > 0:
            gROOT.cd()
            eventList1 = TEventList("eventList1", "eventList1")

            chain = self.samples_RD['chain']
            chain.SetEventList(0)
            print "@@@ Processing real data @@@"
            cuts = []
            for i, (cutName, (cut, plots, weight)) in enumerate(self.cutsteps.iteritems()):
                print "@@ Processing %s " % cutName
                cuts.append("(%s)" % cut)
                cut = "&&".join(cuts)
                for h in hCutFlows.values()+hCutFlowsRaw.values():
                    h.GetXaxis().SetBinLabel(i+1, cutName)

                dout = fout.GetDirectory(cutName)
                if dout == None: dout = fout.mkdir(cutName)
                dout.cd()

                chain.SetProof(False)
                eventList2 = TEventList("eventList2", "eventList2")
                nEvent = chain.Draw('>>eventList2', cut)
                eventList1 = eventList2.Clone("eventList1")
                chain.SetEventList(eventList1)
                hCutFlowsRaw["RD"].AddBinContent(i, nEvent)

                if self.prf != None: chain.SetProof(True)
                for plotName in plots:
                    if plotName not in self.plots: continue
                    varExpr, axisTitles, hArgs = self.plots[plotName]
                    h = TH1D("h%s_%s" % (plotName, "RD"), "%s;%s" % ("Data", axisTitles), *hArgs)
                    print '@@@@ Projecting %s' % plotName
                    chain.Project(h.GetName(), varExpr)
                    #chain.Project(h.GetName(), varExpr, cut)
                    h.Write()

            eventList1, eventList2 = None, None

        ## Signal
        print "@@@ Processing signal MC"
        for sInfo in self.samples_sig.values():
            for ssName, ssInfo in sInfo['subsamples'].iteritems():
                gROOT.cd()
                eventList1 = TEventList("eventList1", "eventList1")

                chain = ssInfo['chain']
                chain.SetEventList(0)
                print "@@@ Processing %s @@@" % ssName
                cuts = []
                for i, (cutName, (cut, plots, weight)) in enumerate(self.cutsteps.iteritems()):
                    print "@@ Processing %s " % cutName
                    cuts.append("(%s)" % cut)
                    cut = "&&".join(cuts)
                    for h in hCutFlows.values()+hCutFlowsRaw.values():
                        h.GetXaxis().SetBinLabel(i+1, cutName)

                    chain.SetProof(False)
                    eventList2 = TEventList("eventList2", "eventList2")
                    nEvent = chain.Draw('>>eventList2', cut)
                    eventList1 = eventList2.Clone("eventList1")
                    chain.SetEventList(eventList1)
                    hCutFlowsRaw[ssName].AddBinContent(i, nEvent)

                    dout = fout.GetDirectory(cutName)
                    if dout == None: dout = fout.mkdir(cutName)

                    if self.prf != None: chain.SetProof(True)
                    for plotName in plots:
                        if plotName not in self.plots: continue
                        varExpr, axisTitles, hArgs = self.plots[plotName]
                        dout.cd()
                        h = TH1D("h%s_%s" % (plotName, ssName), "%s;%s" % (ssName, axisTitles), *hArgs)
                        print '@@@@ Projecting %s' % plotName
                        chain.Project(h.GetName(), varExpr, weight)
                        h.Write()

                eventList1, eventList2 = None, None

        ## Background
        print "@@@ Processing background MC"
        for sInfo in self.samples_bkg.values():
            for ssName, ssInfo in sInfo['subsamples'].iteritems():
                gROOT.cd()
                eventList1 = TEventList("eventList1", "eventList1")

                chain = ssInfo['chain']
                chain.SetEventList(0)
                print "@@@ Processing %s @@@" % ssName
                cuts = []
                for i, (cutName, (cut, plots, weight)) in enumerate(self.cutsteps.iteritems()):
                    print "@@ Processing %s" % cutName
                    cuts.append("(%s)" % cut)
                    cut = "&&".join(cuts)
                    for h in hCutFlows.values()+hCutFlowsRaw.values():
                        h.GetXaxis().SetBinLabel(i+1, cutName)

                    chain.SetProof(False)
                    eventList2 = TEventList("eventList2", "eventList2")
                    nEvent = chain.Draw('>>eventList2', cut)
                    eventList1 = eventList2.Clone("eventList1")
                    chain.SetEventList(eventList1)
                    hCutFlowsRaw[ssName].AddBinContent(i, nEvent)

                    dout = fout.GetDirectory(cutName)
                    if dout == None: dout = fout.mkdir(cutName)

                    if self.prf != None: chain.SetProof(True)
                    for plotName in plots:
                        if plotName not in self.plots: continue
                        varExpr, axisTitles, hArgs = self.plots[plotName]
                        dout.cd()
                        h = TH1D("h%s_%s" % (plotName, ssName), "%s;%s" % (ssName, axisTitles), *hArgs)
                        print '@@@@ Projecting %s' % plotName
                        chain.Project(h.GetName(), varExpr, weight)
                        h.Write()

                eventList1, eventList2 = None, None

        for h in hCutFlows.values()+hCutFlowsRaw.values():
            fout.cd()
            h.Write()

        gROOT.SetBatch(isBatch)

class PlotMaker(HistInfo):
    def __init__(self, prefix, fName, hInfo, option):
        HistInfo.__init__(self, hInfo)
        doRatio = False
        if 'doRatio' in option: doRatio = option['doRatio']
        aspectRatio = 1
        if 'aspectRatio' in option: aspectRatio = option['aspectRatio']
        innerWidth = 400
        legEntryWidth = 0.25
        nCol = 1

        gROOT.ProcessLine(".L ../FlatTuple/tdrstyle.C")
        gROOT.ProcessLine("setTDRStyle();")
        gStyle.SetOptTitle(0)
        gStyle.SetOptStat(0)
        gStyle.SetPadLeftMargin(0.15)
        gStyle.SetPadRightMargin(0.05)
        gStyle.SetLabelSize(0.05, "X")
        gStyle.SetLabelSize(0.05, "Y")
        gStyle.SetTitleOffset(1, "X")
        gStyle.SetTitleOffset(1.2, "Y")

        self.canvases = {}
        self.objs = {}

        innerWidth = 400
        h2Ratio = 0.3

        self.leftMargin = 0.17
        self.pad1_topMargin = 0.2
        width = innerWidth*(1.+self.leftMargin+0.05)
        height, height2 = 0, 0
        if not doRatio:
            self.pad1_bottomMargin = 0.15
            height = innerWidth/aspectRatio*(1.+self.pad1_topMargin+self.pad1_bottomMargin)
        else:
            self.pad2_topMargin = 0.05
            self.pad1_bottomMargin, self.pad2_bottomMargin = 0.05, 0.42
            height = innerWidth/aspectRatio*(1.+self.pad1_topMargin+self.pad1_bottomMargin)
            height2 = innerWidth/aspectRatio*h2Ratio*(1.+self.pad1_topMargin+self.pad2_bottomMargin)
            height += height2

            self.pad2_xLabelSize = gStyle.GetLabelSize("X")/h2Ratio*0.8
            self.pad2_yLabelSize = gStyle.GetLabelSize("Y")/h2Ratio*0.8
            self.pad2_xTitleSize = gStyle.GetTitleSize("X")/h2Ratio*0.9
            self.pad2_yTitleSize = gStyle.GetTitleSize("Y")/h2Ratio*0.7
            self.pad2_xTitleOffset = gStyle.GetTitleOffset("X")
            self.pad2_yTitleOffset = gStyle.GetTitleOffset("Y")/h2Ratio*0.119

        if 'width' in option:
            scale = option['width']/width
            width *= scale
            height *= scale
            height2 *= scale

        height, width = int(height), int(width)
        height2 = int(height2)

        gROOT.cd()

        nLegItems = len(self.samples_bkg)+len(self.samples_sig)
        if len(self.samples_RD['fNames']) > 0: nLegItems += 1
        nRow = int(1.*nLegItems/nCol+0.5)
        self.leg = TLegend(1-gStyle.GetPadRightMargin()-nCol*legEntryWidth, 1-gStyle.GetPadTopMargin()-0.1-nRow*0.02,
                           1-gStyle.GetPadRightMargin()-0.05, 1-gStyle.GetPadTopMargin()-0.05)
        self.leg.SetFillStyle(0)
        self.leg.SetBorderSize(0)
        for title, info in self.samples_bkg.iteritems():
            entry = self.leg.AddEntry(0, title, "f")
            entry.SetFillColor(info['color'])
            entry.SetFillStyle(1111)
        for title, info in self.samples_sig.iteritems():
            entry = self.leg.AddEntry(0, title, "l")
            entry.SetLineColor(info['color'])
            entry.SetFillStyle(0)
        if len(self.samples_RD['fNames']) > 0:
            entry = self.leg.AddEntry(0, "Data", "lp")
            entry.SetLineColor(kBlack)
            entry.SetMarkerColor(kBlack)

        fin = TFile(fName)
        #if not os.path.isdir(prefix): os.makedirs(prefix)

        for stepName0, (cut, plots, weight) in self.cutsteps.iteritems():
            stepName = stepName0
            if prefix != "": stepName = "%s_%s" % (prefix, stepName0)
            self.canvases[stepName] = []
            self.objs[stepName] = []

            for plotName0 in plots:
                if plotName0 not in self.plots: continue
                hFrame = None
                h_dat, h_bkg, hs_sig, hs_bkg = None, None, None, None
                xTitle, yTitle = None, None
                maxY = 0

                plotName = plotName0
                if prefix != "": plotName = "%s_%s" % (prefix, plotName0)

                ## Get the data histogram
                h = fin.Get("%s/h%s_%s" % (stepName0, plotName, "RD"))
                if h != None:
                    gROOT.cd()
                    h_dat = h.Clone()
                    h_dat.SetName("h%s_%s_dat" % (stepName, plotName))
                    if xTitle == None: xTitle = h.GetXaxis().GetTitle()
                    if yTitle == None: yTitle = h.GetYaxis().GetTitle()
                    maxY = max(maxY, h_dat.GetMaximum())
                    if hFrame == None: hFrame = h.Clone()
                    h_dat.AddBinContent(h_dat.GetNbinsX(), h_dat.GetBinContent(h_dat.GetNbinsX()+1))

                for i, (sTitle, sInfo) in enumerate(self.samples_sig.iteritems()):
                    hout = None
                    for ssName, ssInfo in sInfo['subsamples'].iteritems():
                        h = fin.Get("%s/h%s_%s" % (stepName0, plotName0, ssName))
                        if xTitle == None: xTitle = h.GetXaxis().GetTitle()
                        if yTitle == None: yTitle = h.GetYaxis().GetTitle()
                        gROOT.cd()
                        if hout == None:
                            hout = h.Clone()
                            hout.SetName("h%s_%d" % (plotName, i))
                            hout.SetLineColor(sInfo['color'])
                            hout.SetFillStyle(0)
                            #hout.SetOption("hist")
                            hout.Reset()
                        if hs_sig == None:
                            hs_sig = THStack("hs%s_%s_sig" % (stepName, plotName), "%s sig" % plotName)
                        hout.Add(h, self.lumi*1000*ssInfo['xsec']/ssInfo['evt'])
                    if hout == None: continue
                    if hFrame == None: hFrame = h.Clone()
                    hout.AddBinContent(hout.GetNbinsX(), hout.GetBinContent(hout.GetNbinsX()+1))
                    hs_sig.Add(hout)
                    maxY = max(maxY, hout.GetMaximum())
                    self.objs[stepName].append(hout)

                for i, (sTitle, sInfo) in enumerate(self.samples_bkg.iteritems()):
                    hout = None
                    for ssName, ssInfo in sInfo['subsamples'].iteritems():
                        h = fin.Get("%s/h%s_%s" % (stepName0, plotName0, ssName))
                        if xTitle == None: xTitle = h.GetXaxis().GetTitle()
                        if yTitle == None: yTitle = h.GetYaxis().GetTitle()
                        gROOT.cd()
                        if hout == None:
                            hout = h.Clone()
                            hout.SetName("h%s_%d" % (plotName, i))
                            hout.SetLineColor(sInfo['color'])
                            hout.SetFillColor(sInfo['color'])
                            hout.SetFillStyle(1111)
                            #hout.SetOption("hist")
                            hout.Reset()
                        if h_bkg == None:
                            h_bkg = h.Clone()
                            h_bkg.SetName("h%s_%s_bkg" % (stepName, plotName))
                            h_bkg.Reset()
                        hout.Add(h, self.lumi*1000*ssInfo['xsec']/ssInfo['evt'])
                        h_bkg.Add(h, self.lumi*1000*ssInfo['xsec']/ssInfo['evt'])
                        hout.AddBinContent(hout.GetNbinsX(), hout.GetBinContent(hout.GetNbinsX()+1))
                        h_bkg.AddBinContent(h_bkg.GetNbinsX(), h_bkg.GetBinContent(h_bkg.GetNbinsX()+1))
                        if hs_bkg == None:
                            hs_bkg = THStack("hs%s_%s_bkg" % (stepName, plotName), "%s bkg" % plotName)
                    if hout == None: continue
                    if hFrame == None: hFrame = h.Clone()
                    if hs_bkg != None: hs_bkg.Add(hout)
                    self.objs[stepName].append(hout)

                if h_bkg != None: maxY = max(maxY, h_bkg.GetMaximum())

                c = TCanvas("c%s_%s" % (stepName, plotName), "%s %s" % (stepName, plotName), width, height)
                if not doRatio:
                    c.cd()
                    c.SetBottomMargin(self.pad1_bottomMargin)
                elif h_dat != None:
                    c.Divide(1,2)

                    pad2 = c.cd(2)
                    pad2.SetPad(0, 0, 1, 1.*height2/height)
                    pad2.SetTopMargin(0.05)
                    pad2.SetBottomMargin(self.pad2_bottomMargin)

                    h_ratio = h_dat.Clone("hr%s_%s" % (stepName, plotName))
                    h_ratio.SetTitle("ratio;%s;MC/Data" % (h_dat.GetXaxis().GetTitle()))
                    h_dat.SetTitle("Data;;%s" % (h_dat.GetYaxis().GetTitle()))
                    h_ratio.GetXaxis().SetLabelSize(self.pad2_xLabelSize)
                    h_ratio.GetXaxis().SetTitleSize(self.pad2_xTitleSize)
                    h_ratio.GetXaxis().SetTitleOffset(self.pad2_xTitleOffset)
                    h_ratio.GetYaxis().SetLabelSize(self.pad2_yLabelSize)
                    h_ratio.GetYaxis().SetTitleSize(self.pad2_yTitleSize)
                    h_ratio.GetYaxis().SetTitleOffset(self.pad2_yTitleOffset)
                    h_ratio.GetYaxis().SetNdivisions(505)
                    h_ratio.Divide(h_bkg, h_dat)
                    h_ratio.SetMinimum(0.0)
                    h_ratio.SetMaximum(3.0)
                    h_ratio.Draw()

                    self.objs[stepName].append(h_ratio)

                    pad1 = c.cd(1)
                    pad1.SetPad(0, 1.*height2/height, 1, 1)
                    pad1.SetBottomMargin(self.pad1_bottomMargin)

                ## Draw dummy histogram to with proper y range
                if xTitle == None: xTitle = ""
                if yTitle == None: yTitle = ""
                hFrame.Reset()
                hFrame.SetMaximum(1.2*maxY)
                hFrame.SetTitle("%s;%s;%s" % (h.GetTitle(), xTitle, yTitle))
                hFrame.Draw("")
                ### Draw all histograms
                if hs_bkg != None: hs_bkg.Draw("hist,same")
                if h_bkg  != None: h_bkg.Draw("hist,same")
                if hs_sig != None: hs_sig.Draw("hist,nostack,same")
                if h_dat  != None: h_dat.Draw("esame")
                self.leg.Draw()

                gPad.RedrawAxis()
                c.Modified(True)
                c.Update()

                self.objs[stepName].extend([hFrame, h_dat, h_bkg, hs_bkg, hs_sig])
                self.canvases[stepName].append(c)

