#!/usr/bin/env python
from ROOT import *

TMVA.Tools.Instance()

fout = TFile("mva.root", "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

loader = TMVA.DataLoader("dataset")
loader.AddVariable("jets_n", "Jet multiplicity", "", "I")
loader.AddVariable("bjets_n", "B-jet multiplicity", "", "I")
loader.AddVariable("kin_chi2", "Fit Chi2", "", "F")
loader.AddVariable("kin_lepW_m", "Leptonic W mass", "GeV", "F")
loader.AddVariable("kin_lepT_m", "Leptonic T mass", "GeV", "F")
loader.AddVariable("kin_hadT_m", "Hadronic T mass", "GeV", "F")
loader.AddVariable("kin_hadW12_m", "W12 mass", "GeV", "F")
loader.AddVariable("kin_hadW23_m", "W23 mass", "GeV", "F")
#loader.AddVariable("kin_lepB_CSV:=max(0,kin_lepB_CSV)", "lepB CSV", "", "F")
#loader.AddVariable("kin_hadB_CSV:=max(0,kin_hadB_CSV)", "hadB CSV", "", "F")
#loader.AddVariable("kin_hadJ1_CSV:=max(0,kin_hadJ1_CSV)", "hadJ1 CSV", "", "F") # F for the CMS Tree
#loader.AddVariable("kin_hadJ2_CSV:=max(0,kin_hadJ1_CSV)", "hadJ2 CSV", "", "F") # F for the CMS Tree
#loader.AddVariable("kin_lepB_CvsB", "Leptonic B CvsB", "", "F")
#loader.AddVariable("kin_hadB_CvsB", "Hadronic B CvsB", "", "F")
#loader.AddVariable("kin_hadJ1_CvsB", "Leptonic J1 CvsB", "", "F")
#loader.AddVariable("kin_hadJ2_CvsB", "Leptonic J2 CvsB", "", "F")
#loader.AddVariable("kin_lepB_CvsL", "Leptonic B CvsL", "", "F")
#loader.AddVariable("kin_hadB_CvsL", "Hadronic B CvsL", "", "F")
#loader.AddVariable("kin_hadJ1_CvsL", "Hadronic J1 CvsL", "", "F")
#loader.AddVariable("kin_hadJ2_CvsL", "Hadronic J2 CvsL", "", "F")
loader.AddVariable("kin_nbjetInFit:=kin_bjetcode/10+kin_bjetcode%10", "B-jet multiplicity in the fit", "", "I")

for name in ["lepB", "lepT", "hadB", "hadJ1", "hadJ2", "lep", "nu", "hadW12", "hadW23"]:
    for var in ["pt", "eta", "phi"]:
        loader.AddVariable("kin_%s_%s" % (name, var), "%s %s" % (name, var), "", "F")
for name in ["lepB", "hadB", "hadJ1", "hadJ2"]:
    loader.AddVariable("kin_%s_m" % name, "%s mass" % name, "GeV", "F")

loader.AddSpectator("kin_bjetcode", "event category by nbjet in the fit")
loader.AddSpectator("vertex_n", "nVertex")

## Load input files
fsig = TFile("../FlatTuple/ntuple_2Dfit_FCNC/delphes_FCNC.root")
fbkg = TFile("../FlatTuple/ntuple_2Dfit_FCNC/delphes_ttbb.root")
#fsig = TFile("../FlatTuple/ntuple_2Dfit_FCNC/cmsTuple_FCNC.root")
#fbkg = TFile("../FlatTuple/ntuple_2Dfit_FCNC/cmsTuple_ttbb.root")
tsig = fsig.Get("tree")
tbkg = fbkg.Get("tree")

loader.AddSignalTree(tsig, 1.0)
loader.AddBackgroundTree(tbkg, 1.0)
#loader.SetBackgroundWeightExpression( "weight" );

sigCut = TCut("kin_bjetcode%10 >= 1")
bkgCut = TCut("kin_bjetcode%10 >= 1")
loader.PrepareTrainingAndTestTree(
    sigCut, bkgCut,
    #"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V"
    "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V"
)

methods = [
    [TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart"],
    [TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator"],
    [TMVA.Types.kMLP, "MLP2N", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+N:TestRate=5:!UseRegulator"],
    [TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
    [TMVA.Types.kBDT, "BDT5", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
    [TMVA.Types.kBDT, "BDT7", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=7:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
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
    ["DNN5", [
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNNReLU", [
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNNReLU10", [
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],

        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ["DNNTANH10", [
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],

        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
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

TMVA.TMVAGui("mva.root")

