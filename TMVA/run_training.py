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
    ftnName, nX, nY = mvaType0.split('_')[1:]
    nX, nY = int(nX), int(nY)
hasCUDA = os.path.exists('/usr/bin/nvidia-smi') ## Can use GPU acceleration

from ROOT import *
TMVA.Tools.Instance()
if mvaAlgo in ('DNN', 'Keras'):
    TMVA.PyMethodBase.PyInitialize()

fout = TFile("mva_%s.root" % (mvaType0), "recreate")

factory = TMVA.Factory("TMVAClassification", fout,
                       "!V:Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification" )
                       #"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification" )

loader = TMVA.DataLoader("mva")
intVars = [
    "jets_n", "bjets_n",
]
floatVars = [
    "lepton_pt", "met_pt", "met_phi",
]
for name in ["lepT", "lepB", "lepW", "hadT", "hadJ1", "hadJ2", "hadJ3"]:
    floatVars.append("%s_m" % name)
    floatVars.append("%s_pt" % name)
    floatVars.append("%s_eta" % name)
    floatVars.append("%s_phi" % name)
#for name in ["hadJ1", "hadJ2", "hadJ3", "lepB"]:
#for name in ["hadJ3", "lepB"]:
#    floatVars.append("%s_bTag" % name)
for name in ["hadW12", "hadW23", "hadW13"]:
    floatVars.append("%s_m" % name)
    floatVars.append("%s_pt" % name)
    floatVars.append("%s_dR" % name)
for var in intVars:
    loader.AddVariable(var, "I")
for var in floatVars:
    loader.AddVariable(var, "F")
nVars = len(intVars)+len(floatVars)

loader.AddSpectator("theta1", "Theta1")
loader.AddSpectator("theta2", "Theta2")
#loader.AddSpectator("bjetcode", "event category by nbjet in the fit")
loader.AddSpectator("genMatch", "Gen match")
loader.AddSpectator("vertex_n", "nVertex")

## Load input files
fsig = TFile("../FlatTuple/delphes_FCNC.root")
fbkg = TFile("../FlatTuple/delphes_ttbb.root")
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
    "nTrain_Signal=30000:nTrain_Background=30000:SplitMode=Random:NormMode=NumEvents:!V"
)

if mvaType0.split('_')[0] == "BDT":
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

    bdtOption = "!H:!V:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex"
    bdtgOption = "!H:!V:BoostType=Grad:UseBaggedBoost:SeparationType=GiniIndex"
    for nCuts in reversed([20, 30, 40, 50]):#, 60, 80]):
        for maxDepth in [2,3,4]:
            for minNodeSize in [2,2.5,3,5,7]:
                for nTree in reversed([500, 850, 1000, 1500, 2000]):
                    opt = [bdtOption, "NTrees=%d" % nTree,
                           "nCuts=%d" % nCuts, "MaxDepth=%d" % maxDepth, "MinNodeSize=%g%%" % minNodeSize]
                    suffix = "nCuts%d_maxDepth%d_minNode%g_nTree%d" % (nCuts, maxDepth, minNodeSize, nTree)
                    suffix = suffix.replace('.','p')
                    methods.append([TMVA.Types.kBDT, "BDT_%s" % suffix, ":".join(opt)])
                    methods.append([TMVA.Types.kBDT, "BaggedBDT_%s" % suffix, ":".join(opt)+":UseBaggedBoost:BaggedSampleFraction=0.5"])

                    for shrink in [0.05, 0.1, 0.2]:
                        for bagFrac in [0.04, 0.1, 0.2]:
                            optg = [bdtgOption] + opt[1:] + ["Shrinkage=%g" % shrink, "BaggedSampleFraction=%g" %bagFrac]
                            suffix = "nCuts%d_maxDepth%d_minNode%g_nTree%d_shrink%g_bag%g" % (nCuts, maxDepth, minNodeSize, nTree, shrink, bagFrac)
                            suffix = suffix.replace('.','p')
                            methods.append([TMVA.Types.kBDT, "BDTG_%s" % suffix, ":".join(optg)])
    if bdtParIndex.isdigit():
        bdtIndex = int(bdtParIndex)
        factory.BookMethod(loader, *methods[bdtIndex])
    else:
        print "!!!Number of BDT methods to scan = ", len(methods)
        sys.exit()
    #except:
    #    for m in methods: factory.BookMethod(loader, *m)

elif mvaType0.split('_', 1)[0] == "DNN":
    # For the DNN
    name = "DNN_%s_X%d_Y%d" % (ftnName, nX, nY)

    dnnOpts = ["!H:V:ErrorStrategy=CROSSENTROPY:WeightInitialization=XAVIER:VarTransform=D,G,N"]
    if hasCUDA:
        dnnOpts.append("Architecture=GPU")
    else:
        dnnOpts.append("Architecture=CPU")

    ## Configure layout
    layers = "Layout="+(",".join(["%s|%d" % (ftnName, nX) for i in range(nY)]))+",LINEAR"
    dnnOpts.append(layers)

    ## Configure learning strategy in steps
    trainStrategy0 = ["Repetitions=3","ConvergenceSteps=20","Regularization=L2","BatchSize=100","TestRepetitions=5","Multithreading=True"]
    dropConfig = ["DropConfig=" + ("+".join(["0.0"]+["0.5" for i in range(nY-2)]+["0.0"]))]
    trainSteps = [
        ["LearningRate=1e-1","WeightDecay=1e-2","Momentum=0.9"],
        ["LearningRate=1e-2","WeightDecay=1e-3","Momentum=0.9"],
        ["LearningRate=1e-3","WeightDecay=1e-3","Momentum=0.8"],
        ["LearningRate=1e-3","WeightDecay=1e-3","Momentum=0.5"],
    ]
    dnnOpts.append("TrainingStrategy=%s" % ("|".join([",".join(trainStrategy0+strategy+dropConfig) for strategy in trainSteps])))
 
    factory.BookMethod(loader, TMVA.Types.kDNN, name, ":".join(dnnOpts))

elif mvaType0.split('_', 1)[0] == "Keras":
    import google.protobuf
    import keras
    import tensorflow as tf
    if hasCUDA:
        config = tf.ConfigProto()
        config.gpu_options.visible_device_list = "1"
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        config = tf.ConfigProto(
            device_count = {"GPU":0}
        )
    keras.backend.tensorflow_backend.set_session(tf.Session(config=config))

    if 'ReLU' == ftnName: activation = 'relu'
    else: activation = 'tanh'

    #init='glorot_uniform'
    init='glorot_normal'
    #init='random_normal'

    model = keras.models.Sequential()
    #model.add(keras.layers.core.Dense(nX, activation=activation, input_dim=nVars,
    #                                  kernel_regularizer=keras.regularizers.l2(1e-2),
    #                                  kernel_initializer=init, bias_initializer='zeros'))
    model.add(keras.layers.core.Dense(nX, input_dim=nVars,
                                      kernel_regularizer=keras.regularizers.l2(1e-2)))

    for i in range(nY):
        model.add(keras.layers.normalization.BatchNormalization())
        model.add(keras.layers.core.Dense(nX, activation=activation,
                                          kernel_initializer=init, bias_initializer='zeros'))
        if i != nY-1: model.add(keras.layers.core.Dropout(0.5))
    model.add(keras.layers.normalization.BatchNormalization())
    model.add(keras.layers.core.Dense(2, activation='softmax'))

    #optimizer = keras.optimizers.SGD(lr=1e-3, decay=1e-9, momentum=0.5, nesterov=True)
    optimizer = keras.optimizers.Adam(lr=1e-3, decay=1e-3, beta_1=0.9, beta_2=0.999, epsilon=1e-8)
    model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['binary_accuracy'])
    #model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    if not os.path.exists("model"): os.mkdir("model")
    modelFile = 'model/model_%s_X%d_Y%d.h5' % (ftnName, nX, nY)
    model.save(modelFile)
    #model.summary()

    factory.BookMethod(loader, TMVA.Types.kPyKeras, 'Keras_%s_X%d_Y%d' % (ftnName, nX, nY),
                       "!H:V:FilenameModel=%s:NumEpochs=50:BatchSize=100:VarTransform=D,G,N" % modelFile)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
fout.Close()

#TMVA.TMVAGui("mva_%s_%s_%s.root" % (sampleType0, suffix, mvaType0))

