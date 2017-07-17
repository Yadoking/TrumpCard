#!/usr/bin/env python
from ROOT import *

TMVA.Tools.Instance()

fout = TFile("mva.root", "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

loader = TMVA.DataLoader("dataset")
loader.AddVariable("jets_n", "Jet multiplicity", "", "I")
loader.AddVariable("bjets_n", "B-jet multiplicity", "", "I")
#loader.AddVariable("m3_chi2", "Fit Chi2", "", "F")
loader.AddVariable("m3_lepW_m", "Leptonic W mass", "GeV", "F")
loader.AddVariable("m3_lepT_m", "Leptonic T mass", "GeV", "F")
loader.AddVariable("m3_hadT_m", "Hadronic T mass", "GeV", "F")
loader.AddVariable("m3_hadW12_m", "W12 mass", "GeV", "F")
loader.AddVariable("m3_hadW23_m", "W23 mass", "GeV", "F")
loader.AddVariable("m3_hadW13_m", "W13 mass", "GeV", "F")
loader.AddVariable("m3_hadW12_dR", "W12 deltaR", "", "F")
loader.AddVariable("m3_hadW23_dR", "W23 deltaR", "", "F")
loader.AddVariable("m3_hadW13_dR", "W13 deltaR", "", "F")
loader.AddVariable("m3_lepB_CSV:=max(0,m3_lepB_CSV)", "lepB CSV", "", "F")
loader.AddVariable("m3_hadB_CSV:=max(0,m3_hadB_CSV)", "hadB CSV", "", "F")
loader.AddVariable("m3_hadJ1_CSV:=max(0,m3_hadJ1_CSV)", "hadJ1 CSV", "", "F") # F for the CMS Tree
loader.AddVariable("m3_hadJ2_CSV:=max(0,m3_hadJ1_CSV)", "hadJ2 CSV", "", "F") # F for the CMS Tree
loader.AddVariable("m3_lepB_CvsB", "Leptonic B CvsB", "", "F")
loader.AddVariable("m3_hadB_CvsB", "Hadronic B CvsB", "", "F")
loader.AddVariable("m3_hadJ1_CvsB", "Leptonic J1 CvsB", "", "F")
loader.AddVariable("m3_hadJ2_CvsB", "Leptonic J2 CvsB", "", "F")
loader.AddVariable("m3_lepB_CvsL", "Leptonic B CvsL", "", "F")
loader.AddVariable("m3_hadB_CvsL", "Hadronic B CvsL", "", "F")
loader.AddVariable("m3_hadJ1_CvsL", "Hadronic J1 CvsL", "", "F")
loader.AddVariable("m3_hadJ2_CvsL", "Hadronic J2 CvsL", "", "F")
#loader.AddVariable("m3_nbjetInHad:=m3_bjetcode%10", "B-jet multiplicity in the hadronic top", "", "I")
loader.AddVariable("m3_theta1", "Theta1", "", "F")
loader.AddVariable("m3_theta2", "Theta2", "", "F")

for name in ["lepB", "lepT", "hadB", "hadJ1", "hadJ2", "hadW12", "hadW23", "hadW13"]:
    loader.AddVariable("m3_%s_pt" % name, "%s pt" % name, "GeV", "F")
for name in ["lepB", "hadB", "hadJ1", "hadJ2"]:
    loader.AddVariable("m3_%s_m" % name, "%s mass" % name, "GeV", "F")

#loader.AddSpectator("m3_bjetcode", "event category by nbjet in the fit")
loader.AddSpectator("vertex_n", "nVertex")

## Load input files
#fsig = TFile("../FlatTuple/delphes_FCNC.root")
#fbkg = TFile("../FlatTuple/delphes_ttbb.root")
fsig = TFile("../FlatTuple/cmsTuple_FCNC.root")
fbkg = TFile("../FlatTuple/cmsTuple_ttbb.root")
tsig = fsig.Get("tree")
tbkg = fbkg.Get("tree")

loader.AddSignalTree(tsig, 1.0)
loader.AddBackgroundTree(tbkg, 1.0)
#loader.SetBackgroundWeightExpression( "weight" );

#sigCut = TCut("")
#bkgCut = TCut("")
sigCut = TCut("m3_bjetcode%10 >= 2")
bkgCut = TCut("m3_bjetcode%10 >= 2")
loader.PrepareTrainingAndTestTree(
    sigCut, bkgCut,
    #"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V"
    "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V"
)

methods = [
    [TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart"],
    [TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator"],
    #[TMVA.Types.kMLP, "MLP2N", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+N:TestRate=5:!UseRegulator"],
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
    ["DNNReLU", [
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["ReLU|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNN_X50", [
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNN_X50_Y5", [
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNNReLU_X50_Y10", [
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],

        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
        ["ReLU|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
    ],
    ["DNNTANH_X50_Y10", [
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],

        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|50", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]],
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

