#!/usr/bin/env python
import sys, os
from collections import OrderedDict

if len(sys.argv) < 3:
    print sys.argv[0], "SAMPLETYPE SUFFIX"
    print '  sampleType = ["cmsTuple", "cmsTuple.withBtag", "delphes"]'
    print '  suffix     = ["kin", "m3", "deltaR"]'
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

if len(sys.argv) < 4: mvaType0 = "BDT"
else: mvaType0 = sys.argv[3]

from ROOT import *
import google.protobuf
import keras

TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

fout = TFile("mva_%s_%s_%s.root" % (sampleType0, suffix, mvaType0), "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification" )

loader = TMVA.DataLoader("mva_%s_%s" % (sampleType0, suffix))
loader.AddVariable("jets_n", "Jet multiplicity", "", "I")
loader.AddVariable("bjets_n", "B-jet multiplicity", "", "I")
loader.AddVariable("nbjetInHad:=bjetcode%10", "B-jet multiplicity in the hadronic top", "", "I")
loader.AddVariable("theta1", "Theta1", "", "F")
loader.AddVariable("theta2", "Theta2", "", "F")
#loader.AddVariable("chi2", "Fit Chi2", "", "F")
loader.AddVariable("lepton_pt", "lepton pt", "GeV", "F")
loader.AddVariable("met_pt", "met pt", "GeV", "F")
loader.AddVariable("met_dphi", "met dphi", "GeV", "F")
for name in ["lepT", "lepB", "lepW", "hadT", "hadJ1", "hadJ2", "hadB"]:
    loader.AddVariable("%s_m" % name, "%s mass" % name, "GeV", "F")
    loader.AddVariable("%s_pt" % name, "%s pt" % name, "GeV", "F")
    loader.AddVariable("%s_eta" % name, "%s eta" % name, "", "F")
    loader.AddVariable("%s_dphi" % name, "%s dphi" % name, "", "F")
for name in ["hadW12", "hadW23", "hadW13"]:
    loader.AddVariable("%s_m" % name, "%s mass" % name, "GeV", "F")
    loader.AddVariable("%s_pt" % name, "%s pt" % name, "GeV", "F")
    loader.AddVariable("%s_dR" % name, "%s deltaR" % name, "", "F")
    loader.AddVariable("%s_dphi" % name, "%s dphi" % name, "", "F")
if doBtag and sampleType == "cmsTuple":
    for name in ["lepB", "hadJ1", "hadJ2", "hadB"]:
        loader.AddVariable("%s_CSV:=max(0,%s_CSV)" % (name, name), "%s CSV" % name, "", "F")
        loader.AddVariable("%s_CvsB" % name, "%s CvsB" % name, "", "F")
        loader.AddVariable("%s_CvsL" % name, "%s CvsL" % name, "", "F")

#loader.AddSpectator("bjetcode", "event category by nbjet in the fit")
loader.AddSpectator("vertex_n", "nVertex")

## Load input files
fsig = TFile("/home/jhgoh/work/Delphes/TrumpCard/FlatTuple/%s/%s_FCNC.root" % (suffix, sampleType))
fbkg = TFile("/home/jhgoh/work/Delphes/TrumpCard/FlatTuple/%s/%s_tt.root" % (suffix, sampleType))
tsig = fsig.Get("tree")
tbkg = fbkg.Get("tree")

loader.AddSignalTree(tsig, 0.01)
loader.AddBackgroundTree(tbkg, 0.5)
#loader.SetBackgroundWeightExpression( "weight" );

sigCut = TCut("jets_n >= 4 && bjets_n >= 3")
bkgCut = TCut("jets_n >= 4 && bjets_n >= 3")
loader.PrepareTrainingAndTestTree(
    sigCut, bkgCut,
    #"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V"
    "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V"
)

if mvaType0 == "BDT":
    methods = [
    #    [TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart"],
    #    [TMVA.Types.kCuts, "CutsPCA", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA"],
    #    [TMVA.Types.kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50"],
    #    [TMVA.Types.kLikelihood, "LikelihoodPCA", "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"],

    #    [TMVA.Types.kKNN, "KNN", "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"],
    #    [TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm"],

    #    [TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator"],
    #    [TMVA.Types.kMLP, "MLP2N", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+N:TestRate=5:!UseRegulator"],

        [TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
        [TMVA.Types.kBDT, "BDT2k", "!H:!V:NTrees=2000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],
        [TMVA.Types.kBDT, "BDT300", "!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"],

    ]
    for m in methods: factory.BookMethod(loader, *m)
elif mvaType0.split('_', 1)[0] == "DNN":
    mvaType = mvaType0.split('_', 1)[-1]

    # For the DNN
    dnnCommonOpt = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM"
    dnnCommonOpt += ":Architecture=CPU"
    trainingCommonOpt = ["Repetitions=1", "ConvergenceSteps=20", "Multithreading=True", "Regularization=L2",
                         "WeightDecay=1e-4", "BatchSize=256", "TestRepetitions=10",]

    dnnLayouts = OrderedDict()
    dnnLayouts["DNN_Default"] = [
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-1","Momentum=0.9","DropConfig=0.0+0.5+0.5+0.5"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-2","Momentum=0.9","DropConfig=0.0+0.0+0.0+0.0"]],
        ["TANH|128", trainingCommonOpt+["LearningRate=1e-3","Momentum=0.0","DropConfig=0.0+0.0+0.0+0.0"]]
    ]
    for nX in [512, 256, 128, 64, 32, 16]:
        for nY in range(3,11):
            layers = []
            momConfig = "Momentum=0.9"
            rateConfig = "LearningRate=1e-1"
            dropConfig = "DropConfig=0.0+0.5+0.5+0.5"
            if mvaType != "Mix": ftnName = mvaType
            else: ftnName = "ReLU"

            for i in range(nY):
                if i == nY-1:
                    momConfig = "Momentum=0.9"
                    rateConfig = "LearningRate=1e-3"
                    dropConfig = "DropConfig=0.0+0.0+0.0+0.0"
                elif i == nY-2:
                    rateConfig = "LearningRate=1e-2"
                if mvaType == "Mix" and i > nY-3: ftnName = "TANH"

                layers.append(["%s|%d" % (ftnName, nX), trainingCommonOpt+[rateConfig, momConfig, dropConfig]])
            dnnLayouts["DNN_%s_X%d_Y%d" % (mvaType, nX, nY)] = layers

    for name, dnnLayout in dnnLayouts.iteritems():
        nNode = int(dnnLayout[-1][0].split('|')[-1])
        dnnOpts = [dnnCommonOpt,
            ("Layout="+("|".join([x[0] for x in dnnLayout]))+("|%d,LINEAR" % nNode)),
            ("TrainingStrategy="+("|".join([",".join(x[1]) for x in dnnLayout]))),
        ]
        factory.BookMethod(loader, TMVA.Types.kDNN, name, ":".join(dnnOpts))
elif mvaType0.split('_', 1)[0] == "Keras":
    #weightInit = keras.initializers.VarianceScaling(scale=0.05, mode='fan_in', distribution='normal', seed=None)
    weightInit = keras.initializers.glorot_uniform()

    model = keras.models.Sequential()
    model.add(keras.layers.core.Dense(300, kernel_initializer=weightInit, activation='relu', W_regularizer=keras.regularizers.l2(1e-5), input_dim=48))
    model.add(keras.layers.core.Dense(500, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(700, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(500, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(500, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(300, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(100, kernel_initializer=weightInit, activation='relu'))
    model.add(keras.layers.core.Dense(2, kernel_initializer=weightInit, activation='softmax'))

    model.compile(loss='categorical_crossentropy', optimizer=keras.optimizers.SGD(lr=0.01), metrics=['accuracy'])
    model.save('model.h5')
    model.summary()

    factory.BookMethod(loader, TMVA.Types.kPyKeras, 'Keras', "!H:V:VarTransform=N:FilenameModel=model.h5:NumEpochs=20:BatchSize=32")

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
fout.Close()

#TMVA.TMVAGui("mva_%s_%s_%s.root" % (sampleType0, suffix, mvaType0))

