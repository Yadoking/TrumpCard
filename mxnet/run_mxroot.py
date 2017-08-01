#!/usr/bin/env python

from ROOT import *
import mxnet as mx
import numpy as np

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

## Implementation of root iterator
class RootIter(mx.io.DataIter):
    def __init__(self, datasets, variables, batch_size):
        self.batch_size = batch_size
        self.cur_batch = -1

        self.im_size = {}
        self._provide_data = []
        self._provide_label = []

        self.datasets = {}
        for name in datasets:
            self.datasets[name] = {}
            self.datasets[name]['chain'] = TChain('tree')
            self.datasets[name]['entry'] = -1
            for fName in datasets[name]:
                self.datasets[name]['chain'].Add(fName)
            self.datasets[name]['nEntries'] = self.datasets[name]['chain'].GetEntries()
        maxEvent = min([self.datasets[name]['nEntries'] for name in datasets.iterkeys()])
        self.max_batch = int(maxEvent/self.batch_size)

        self.variables = variables
        for varClass in variables:
            nChannel = len(variables[varClass])
            if 'im_' not in varClass:
                self._provide_data.append(mx.io.DataDesc(varClass, (self.batch_size, nChannel), dtype='float'))
            else:
                chain = self.datasets.values()[0]['chain']
                chain.GetEntry(0)
                h = getattr(chain, variables[varClass][0])
                width, height = h.GetNbinsX()+2, h.GetNbinsY()+2
                self.im_size[varClass] = width, height
                self._provide_data.append(mx.io.DataDesc(varClass, (self.batch_size, nChannel, width, height), dtype='float'))

        self._provide_label = [mx.io.DataDesc("proc_label", (self.batch_size,), dtype='float')]

    def __iter__(self):
        return self

    def reset(self):
        if not hasattr(self, 'datasets'): return

        for val in self.datasets.itervalues():
            val['entry'] = 0
            val['chain'].GetEntry(0)

    def __next__(self):
        return self.next()

    @property
    def provide_data(self):
        return self._provide_data

    @property
    def provide_label(self):
        return self._provide_label

    def next(self):
        if self.cur_batch >= self.max_batch: raise StopIteration

        out_data = []
        out_label = mx.nd.zeros(self.batch_size*1, dtype='float').reshape((self.batch_size,))
        for varClass, varNames in self.variables.iteritems():
            if 'im_' not in varClass:
                l = len(varNames)
                arr = mx.nd.zeros(self.batch_size*l, dtype='float').reshape((self.batch_size, l))
            else:
                l, (w, h) = len(varNames), self.im_size[varClass]
                arr = mx.nd.zeros(self.batch_size*l*w*h, dtype='float').reshape((self.batch_size, l, w, h))
            out_data.append(arr)

        for i in range(self.batch_size):
            label = np.random.randint(0, len(self.datasets))
            labelName = self.datasets.keys()[label]

            ds = self.datasets[labelName]
            ds['entry'] += 1
            iEntry = ds['entry']
            if iEntry >= ds['nEntries']: raise StopIteration

            out_label[i] = label

            chain = ds['chain']
            chain.GetEntry(iEntry)
            for j, (varClass, varNames) in enumerate(self.variables.iteritems()):
                if 'im_' not in varClass:
                    for k, varName in enumerate(varNames):
                        out_data[j][i][k] = getattr(chain, varName)
                else:
                    for k, varName in enumerate(varNames):
                        h = getattr(chain, varName)
                        buf = h.GetArray()
                        buf.SetSize(h.GetSize())
                        buf = np.array(buf)
                        buf.resize(self.im_size[varClass])
                        out_data[j][i][k] = buf

        self.cur_batch += 1
        return mx.io.DataBatch(out_data, [out_label])

## Start to setup the network
batch_size = 100
data_iter = RootIter(datasets, variables, batch_size)
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

