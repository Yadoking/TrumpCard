#!/usr/bin/env python

import mxnet as mx
import RootIter

## Prepare to load ntuple
datasets = {
    'FCNC':["../FlatTuple/full_shift_rotation_flip/m3/delphes_FCNC.root",],
    "ttbb":["../FlatTuple/full_shift_rotation_flip/m3/delphes_ttbb.root",],
}

## Define variable groups
variables = {
    'ft_0':[ # Variables after full reconstruction
        "jets_n", "bjets_n",
        "kin_bjetcode", "kin_hadW12_m", "kin_hadW23_m", "kin_hadW13_m",
        "kin_hadT_m",
        "kin_lepW_m", "kin_lepT_m", 
        "kin_lepW_pt", "kin_lepW_eta", "kin_lepW_dphi",
        "kin_lepT_pt", "kin_lepT_eta", "kin_lepT_dphi",
        "kin_hadW12_pt", "kin_hadW12_eta", "kin_hadW12_dphi", "kin_hadW12_dR",
        "kin_hadW23_pt", "kin_hadW23_eta", "kin_hadW23_dphi", "kin_hadW23_dR",
        "kin_hadW13_pt", "kin_hadW13_eta", "kin_hadW13_dphi", "kin_hadW13_dR",
        "kin_hadT_pt", "kin_hadT_eta", "kin_hadT_dphi",
        "lepton_pt", "lepton_eta", "lepton_phi", "met_pt", "met_dphi",
        "kin_lepB_pt", "kin_lepB_eta", "kin_lepB_dphi", "kin_lepB_m",
        "kin_hadJ1_pt", "kin_hadJ1_eta", "kin_hadJ1_dphi", "kin_hadJ1_m",
        "kin_hadJ2_pt", "kin_hadJ2_eta", "kin_hadJ2_dphi", "kin_hadJ2_m",
        "kin_hadB_pt", "kin_hadB_eta", "kin_hadB_dphi", "kin_hadB_m",
        "kin_theta1", "kin_theta2",
    #    "kin_lepB_CSV", "kin_hadB_CSV", "kin_hadJ1_CSV", "kin_hadJ2_CSV",
    #    "kin_lepB_CvsB", "kin_hadB_CvsB", "kin_hadJ1_CvsB", "kin_hadJ2_CvsB",
    #    "kin_lepB_CvsL", "kin_hadB_CvsL", "kin_hadJ1_CvsL", "kin_hadJ2_CvsL",
    ],
    #'im_0':[
    #    "kin_hJetImage_ch_n", "kin_hJetImage_nh_n", "kin_hJetImage_ph_n",
    #    "kin_hJetImage_ch_pt", "kin_hJetImage_nh_pt", "kin_hJetImage_ph_pt",
    #],
}

## Start to setup the network
batch_size = 100
data_iter = RootIter.RootIter(datasets, variables, batch_size)
data_iter.next()

net = mx.sym.Variable('ft_0')
#net = mx.sym.Variable('im_0')
net = mx.sym.FullyConnected(data=net, name='fc1', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac1', act_type="relu")
net = mx.sym.FullyConnected(data=net, name='fc2', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac2', act_type="relu")
net = mx.sym.FullyConnected(data=net, name='fc3', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac3', act_type="relu")

net = mx.sym.SoftmaxOutput(data=net, name='proc')

print(net.list_arguments())
print(net.list_outputs())

import logging
logging.basicConfig(level=logging.INFO)

mod = mx.mod.Module(symbol=net, data_names=variables.keys(), label_names=['proc_label'])
mod.fit(data_iter,
        initializer=mx.init.Xavier(magnitude=2.),
        optimizer='sgd', optimizer_params=(('learning_rate', 0.1), ),
        eval_metric='acc',  # report accuracy during training
        batch_end_callback = mx.callback.Speedometer(batch_size, 10),
        num_epoch=5)

