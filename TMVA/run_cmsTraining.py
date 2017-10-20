#!/usr/bin/env python
import sys, os
from collections import OrderedDict

if len(sys.argv) < 2:
    print sys.argv[0], "SAMPLETYPE SUFFIX"
    print '  mvaType    = ["BDT", "DNN_TANH", "DNN_ReLU", "Keras_TANH", "Keras_ReLU"]'
    sys.exit(1)

mvaType0 = sys.argv[1]
mvaAlgo = mvaType0.split('_',1)[0]
if mvaAlgo == 'BDT':
    bdtParIndex = mvaType0.split('_')[1]
elif mvaAlgo in ('DNN', 'Keras'):
    ftnName, nX = mvaType0.split('_')[1:]
    nX = int(nX)

from ROOT import *
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

fout = TFile("mva_%s.root" % (mvaType0), "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification" )

loader = TMVA.DataLoader("mva")
intVars = [
    "njets", "nbjets_m", "ncjets_m",
]
floatVars = [
    #"lepDPhi", "transverseMass",
    "missingET",
    ##"bjetmDR","bjetmDEta","bjetmDPhi","dibjetsMass","bjetPt_dibjetsm","cjetPt",
    "jet1pt", "jet1eta",#"jet1phi",
    "jet1m","jet1csv","jet1cvsl","jet1cvsb",
    "jet2pt","jet2eta", #"jet2phi",
    "jet2m","jet2csv","jet2cvsl","jet2cvsb",
    "jet3pt","jet3eta",#"jet3phi",
    "jet3m","jet3csv","jet3cvsl","jet3cvsb",
    "jet4pt","jet4eta",#"jet4phi",
    "jet4m","jet4csv","jet4cvsl","jet4cvsb",
    "DRlepWpt","DRlepWeta","DRlepWdeta",#"DRlepWphi",
    "DRlepWdphi","DRlepWm",
    "DRjet0pt","DRjet0eta",#"DRjet0phi",
    "DRjet0m","DRjet0csv","DRjet0cvsl","DRjet0cvsb",
    "DRjet1pt","DRjet1eta",#"DRjet1phi",
    "DRjet1m","DRjet1csv","DRjet1cvsl","DRjet1cvsb",
    "DRjet2pt","DRjet2eta",#"DRjet2phi",
    "DRjet2m","DRjet2csv","DRjet2cvsl","DRjet2cvsb",
    "DRjet3pt","DRjet3eta",#"DRjet3phi",
    "DRjet3m","DRjet3csv","DRjet3cvsl","DRjet3cvsb",
    "DRjet12pt","DRjet12eta","DRjet12deta",#"DRjet12phi",
    "DRjet12dphi","DRjet12m","DRjet12DR",
    "DRjet23pt","DRjet23eta","DRjet23deta",#"DRjet23phi",
    "DRjet23dphi","DRjet23m",
    "DRjet31pt","DRjet31eta","DRjet31deta",#"DRjet31phi",
    "DRjet31dphi","DRjet31m",
    "DRlepTpt","DRlepTeta","DRlepTdeta",#"DRlepTphi",
    "DRlepTdphi","DRlepTm",
    "DRhadTpt","DRhadTeta",#"DRhadTphi",
    "DRhadTm",
    "DRhadTHbdeta","DRhadTWbdeta",
    "DRhadTHbdphi","DRhadTWbdphi",
]
for var in intVars:
    loader.AddVariable(var, "I")
for var in floatVars:
    loader.AddVariable(var, "F")
nVars = len(intVars)+len(floatVars)

## Load input files
sigs = [
    ("mkNtuple_v5/tmva_AntiTop_Hct.root", 0.1), 
    ("mkNtuple_v5/tmva_Top_Hct.root", 0.1),
]
bkgs = [
    ("mkNtuple_v5/tmva_tbarchannel.root", 0.024575262909),
    ("mkNtuple_v5/tmva_tbarWchannel.root", 0.193026936331),
    ("mkNtuple_v5/tmva_tchannel.root", 0.023844899675),
    ("mkNtuple_v5/tmva_tWchannel.root", 0.190335714074),
    ("mkNtuple_v5/tmva_ttbb.root", 0.0888153017294),
    ("mkNtuple_v5/tmva_ttbj.root", 0.0888153017294),
    ("mkNtuple_v5/tmva_ttcc.root", 0.0888153017294),
    ("mkNtuple_v5/tmva_ttLF.root", 0.0888153017294),
    ("mkNtuple_v5/tmva_tt.root", 0.0888153017294),
]
trees = []
for fName, scale in sigs:
    f = TFile(fName)
    t = f.Get("tmva_tree")
    loader.AddSignalTree(t, scale)
    trees.append([f, t])
for fName, scale in bkgs:
    f = TFile(fName)
    t = f.Get("tmva_tree")
    loader.AddBackgroundTree(t, scale)
    trees.append([f, t])
#loader.SetBackgroundWeightExpression( "weight" );

sigCut = TCut("missingET > 0 && cjetPt > 0 && jet1csv > 0 &&  jet2csv > 0 &&  jet3csv > 0 && jet4csv > 0 && DRlepWpt > 0 && DRjet0csv > 0 && DRjet1csv > 0 && DRjet2csv > 0 && DRjet3csv > 0")
bkgCut = TCut("missingET > 0 && cjetPt > 0 && jet1csv > 0 &&  jet2csv > 0 &&  jet3csv > 0 && jet4csv > 0 && DRlepWpt > 0 && DRjet0csv > 0 && DRjet1csv > 0 && DRjet2csv > 0 && DRjet3csv > 0")

loader.PrepareTrainingAndTestTree(
    sigCut, bkgCut,
    #"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V"
    "nTrain_Signal=30000:nTrain_Background=100000:SplitMode=Random:NormMode=NumEvents:!V"
)

if mvaAlgo == "BDT":
    methods = [
    #    [TMVA.Types.kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart"],
    #    [TMVA.Types.kCuts, "CutsPCA", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA"],
    #    [TMVA.Types.kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50"],
    #    [TMVA.Types.kLikelihood, "LikelihoodPCA", "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"],

    #    [TMVA.Types.kKNN, "KNN", "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"],
    #    [TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=D,G,N"],

    #    [TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=D,G,N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator"],
    #    [TMVA.Types.kMLP, "MLP2N", "H:!V:NeuronType=tanh:VarTransform=D,G,N:NCycles=600:HiddenLayers=N+N:TestRate=5:!UseRegulator"],
    ]

    bdtOption = "!H:!V:BoostType=AdaBoost:AdaBoostBeta=0.5"
    for sepType in ["GiniIndex"]:#, "CrossEntropy"]:
        for nCuts in [20, 5, 10, 30, 40]:#[50, 60, 80]:
            for maxDepth in [2,3,4,5,7,10]:
                for minNodeSize in [2.5, 1,2,3,5,10,20]:
                    for nTree in [850, 500, 1000, 10, 20, 30, 50, 100, 200]:
                        opt = [bdtOption, "NTrees=%d" % nTree, "SeparationType=%s" % sepType,
                               "nCuts=%d" % nCuts, "MaxDepth=%d" % maxDepth, "MinNodeSize=%g%%" % minNodeSize]
                        suffix = "%s_nCuts%d_maxDepth%d_minNode%g_nTree%d" % (sepType, nCuts, maxDepth, minNodeSize, nTree)
                        suffix = suffix.replace('.','p')
                        methods.append([TMVA.Types.kBDT, "BDT_%s" % suffix, ":".join(opt)])
                        methods.append([TMVA.Types.kBDT, "BaggedBDT_%s" % suffix, ":".join(opt)+":UseBaggedBoost:BaggedSampleFraction=0.5"])
    if bdtParIndex.isdigit():
        bdtIndex = int(bdtParIndex)
        factory.BookMethod(loader, *methods[bdtIndex])
    else:
        print "!!!Number of BDT methods to scan = ", len(methods)
        sys.exit()
    #except:
    #    for m in methods: factory.BookMethod(loader, *m)
    

elif mvaAlgo == "DNN":
    # For the DNN
    dnnCommonOpt = "!H:V:ErrorStrategy=CROSSENTROPY:WeightInitialization=XAVIERUNIFORM:VarTransform=D,G,N"
    hasCUDA = False#os.path.exists('/usr/bin/nvidia-smi') ## Can use GPU acceleration
    if hasCUDA:
        dnnCommonOpt += ":Architecture=GPU"
    else:
        dnnCommonOpt += ":Architecture=CPU"
    trainingCommonOpt = ["Repetitions=1", "ConvergenceSteps=20", "Multithreading=True", "Regularization=L2",
                         "WeightDecay=1e-4", "BatchSize=128", "TestRepetitions=10",]

    dnnLayouts = OrderedDict()
    for nY in range(20,0,-1):
        layers = []
        momConfig = "Momentum=0.9"
        rateConfig = "LearningRate=1e-1"
        dropConfig = "DropConfig=0.0+0.5+0.5+0.5"

        for i in range(nY):
            if i == nY-1:
                momConfig = "Momentum=0.9"
                rateConfig = "LearningRate=1e-3"
                dropConfig = "DropConfig=0.0+0.0+0.0+0.0"
            elif i == nY-2:
                rateConfig = "LearningRate=1e-2"

            layers.append(["%s|%d" % (ftnName, nX), trainingCommonOpt+[rateConfig, momConfig, dropConfig]])
        dnnLayouts["DNN_%s_X%d_Y%d" % (ftnName, nX, nY)] = layers

    for name, dnnLayout in dnnLayouts.iteritems():
        nNode = int(dnnLayout[-1][0].split('|')[-1])
        dnnOpts = [dnnCommonOpt,
            ("Layout="+("|".join([x[0] for x in dnnLayout]))+("|%d,LINEAR" % nNode)),
            ("TrainingStrategy="+("|".join([",".join(x[1]) for x in dnnLayout]))),
        ]
        factory.BookMethod(loader, TMVA.Types.kDNN, name, ":".join(dnnOpts))
elif mvaAlgo == "Keras":
    import google.protobuf
    import keras
    hasCUDA = os.path.exists('/usr/bin/nvidia-smi') ## Can use GPU acceleration
    if hasCUDA:
        import tensorflow as tf
        config = tf.ConfigProto()
        config.gpu_options.visible_device_list = "0"
        keras.backend.tensorflow_backend.set_session(tf.Session(config=config))

    if 'ReLU' == ftnName: activation = 'relu'
    else: activation = 'tanh'

    init='glorot_uniform'
    #init='glorot_normal'
    #init='random_normal'

    for nY in range(20,0,-1):
        model = keras.models.Sequential()
        model.add(keras.layers.core.Dense(nX, init=init, activation=activation, kernel_regularizer=keras.regularizers.l2(1e-6), input_dim=nVars))
        model.add(keras.layers.normalization.BatchNormalization())

        for i in range(nY):
            model.add(keras.layers.core.Dropout(0.5))
            model.add(keras.layers.core.Dense(nX, init=init, activation=activation)) 
            model.add(keras.layers.normalization.BatchNormalization())
        model.add(keras.layers.core.Dropout(0.5))
        model.add(keras.layers.core.Dense(nX, init=init, activation='linear'))
        model.add(keras.layers.normalization.BatchNormalization())
        model.add(keras.layers.core.Dense(2, activation='softmax'))

        #optimizer = keras.optimizers.SGD(lr=1e-3, decay=1e-9, momentum=0.5, nesterov=True)
        optimizer = keras.optimizers.Adam(lr=1e-3, beta_1=0.9, beta_2=0.999)
        model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['categorical_accuracy'])
        #model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
        modelFile = 'model_%s_X%d_Y%d.h5' % (ftnName, nX, nY)
        model.save(modelFile)
        #model.summary()

        factory.BookMethod(loader, TMVA.Types.kPyKeras, 'Keras_%s_X%d_Y%d' % (ftnName, nX, nY), "!H:V:FilenameModel=%s:NumEpochs=30:BatchSize=128:VarTransform=D,G,N" % modelFile)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
fout.Close()


