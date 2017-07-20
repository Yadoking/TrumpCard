#!/usr/bin/env python
import sys, os

if len(sys.argv) < 3:
    print sys.argv[0], "SAMPLETYPE SUFFIX"
    print '  sampleType = ["cmsTuple", "cmsTuple.withBtag", "delphes"]'
    print '  suffix     = [kin", "m3", "deltaR"]'
    sys.exit(1)

sampleType0, suffix = sys.argv[1], sys.argv[2]
if suffix not in ["kin", "m3", "deltaR"]:
    print "choose among kin, m3 or deltaR"
    sys.exit(1)
if sampleType0 not in ["cmsTuple", "cmsTuple.withBtag", "delphes"]:
    print "choose sampleType in cmsTuple, cmsTuple.withBtag, delphes"
    sys.exit(1)
if sampleType0 == "cmsTuple.withBtag":
    sampleType = "cmsTuple"
    doBtag = True
else:
    sampleType = sampleType0
    doBtag = False

from ROOT import *

TMVA.Tools.Instance()

fout = TFile("mva_%s_%s.root" % (sampleType0, suffix), "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification" )

loader = TMVA.DataLoader("mva_%s_%s" % (sampleType0, suffix))
loader.AddVariable("jets_n", "Jet multiplicity", "", "I")
loader.AddVariable("bjets_n", "B-jet multiplicity", "", "I")
loader.AddVariable("kin_nbjetInHad:=kin_bjetcode%10", "B-jet multiplicity in the hadronic top", "", "I")
loader.AddVariable("kin_theta1", "Theta1", "", "F")
loader.AddVariable("kin_theta2", "Theta2", "", "F")
#loader.AddVariable("kin_chi2", "Fit Chi2", "", "F")
loader.AddVariable("lepton_pt", "lepton pt", "GeV", "F")
loader.AddVariable("met_pt", "met pt", "GeV", "F")
loader.AddVariable("met_dphi", "met dphi", "GeV", "F")
for name in ["lepT", "lepB", "lepW", "hadT", "hadJ1", "hadJ2", "hadB"]:
    loader.AddVariable("kin_%s_m" % name, "%s mass" % name, "GeV", "F")
    loader.AddVariable("kin_%s_pt" % name, "%s pt" % name, "GeV", "F")
    loader.AddVariable("kin_%s_eta" % name, "%s eta" % name, "", "F")
    loader.AddVariable("kin_%s_dphi" % name, "%s dphi" % name, "", "F")
for name in ["hadW12", "hadW23", "hadW13"]:
    loader.AddVariable("kin_%s_m" % name, "%s mass" % name, "GeV", "F")
    loader.AddVariable("kin_%s_pt" % name, "%s pt" % name, "GeV", "F")
    loader.AddVariable("kin_%s_dR" % name, "%s deltaR" % name, "", "F")
    loader.AddVariable("kin_%s_dphi" % name, "%s dphi" % name, "", "F")
if doBtag and sampleType == "cmsTuple":
    for name in ["lepB", "hadJ1", "hadJ2", "hadB"]:
        loader.AddVariable("kin_%s_CSV:=max(0,kin_%s_CSV)" % (name, name), "%s CSV" % name, "", "F")
        loader.AddVariable("kin_%s_CvsB" % name, "%s CvsB" % name, "", "F")
        loader.AddVariable("kin_%s_CvsL" % name, "%s CvsL" % name, "", "F")

#loader.AddSpectator("kin_bjetcode", "event category by nbjet in the fit")
loader.AddSpectator("vertex_n", "nVertex")

## Load input files
fsig = TFile("/home/jhgoh/work/Delphes/TrumpCard/FlatTuple/%s/%s_FCNC.root" % (suffix, sampleType))
fbkg = TFile("/home/jhgoh/work/Delphes/TrumpCard/FlatTuple/%s/%s_ttbb.root" % (suffix, sampleType))
tsig = fsig.Get("tree")
tbkg = fbkg.Get("tree")

loader.AddSignalTree(tsig, 1.0)
loader.AddBackgroundTree(tbkg, 1.0)
#loader.SetBackgroundWeightExpression( "weight" );

sigCut = TCut("jets_n >= 4 && bjets_n >= 3")
bkgCut = TCut("jets_n >= 4 && bjets_n >= 3")
loader.PrepareTrainingAndTestTree(
    sigCut, bkgCut,
    #"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V"
    "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V"
)

methods = [
#    [TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart"],
#    [TMVA.Types.kCuts, "CutsPCA", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA"],
    [TMVA.Types.kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50"],
#    [TMVA.Types.kLikelihood, "LikelihoodPCA", "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"],

#    [TMVA.Types.kKNN, "KNN", "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"],
#    [TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm"],

#    [TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator"],
#    [TMVA.Types.kMLP, "MLP2N", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+N:TestRate=5:!UseRegulator"],

    [TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
#    [TMVA.Types.kBDT, "BDT5", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
#    [TMVA.Types.kBDT, "BDT7", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=7:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],

#    [TMVA.Types.kBDT, "BDTB", "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20"],
#    [TMVA.Types.kBDT, "BDTB850", "!H:!V:NTrees=850:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20"],


]
for m in methods: factory.BookMethod(loader, *m)

# For the DNN
dnnCommonOpt = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM"
#dnnCommonOpt += ":Architecture=CPU"
trainingCommonOpt = ["Repetitions=1", "ConvergenceSteps=20", "Multithreading=True", "Regularization=L2",
                     "WeightDecay=1e-4", "BatchSize=256", "TestRepetitions=10",]
dnnLayouts = [
    ["DNN", [
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
#    ["DNNReLU", [
#        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
#        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
#    ],
#    ["DNN_X50", [
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
#    ],
#    ["DNN_X50_Y5", [
#        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
#        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
#    ],
#    ["DNNReLU_X50_Y10", [
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
#    ],
#    ["DNNTANH_X50_Y10", [
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
#        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
#    ],
]

for name, dnnLayout in dnnLayouts:
    dnnOpts = [dnnCommonOpt,
        ("Layout="+("|".join([x[0] for x in dnnLayout]))+"|128,LINEAR"),
        ("TrainingStrategy="+("|".join([",".join(x[1]) for x in dnnLayout]))),
    ]
    factory.BookMethod(loader, TMVA.Types.kDNN, name, ":".join(dnnOpts))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
fout.Close()

TMVA.TMVAGui("mva_%s_%s.root" % (sampleType0, suffix))

