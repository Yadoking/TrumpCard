#!/usr/bin/env python

## Example: two helix with shifted by pi/2
## Step1: generate xor-like distributions using ROOT
## Step2: save signal TTree and background TTree
## Step3: Load TTrees using RootIter, custom mxnet-iterator
## Step4: Run the mxnet

import mxnet as mx
import RootIter
from math import sin, cos, pi

###########################################################
## Generate ntuples
import ctypes
from ROOT import *

zmax = 2*(2*pi)
smear = 0.01

hsig = TH3F("hsig", "hsig", 100, -1.5, 1.5, 100, -1.5, 1.5, 100, 0, zmax)
hbkg = TH3F("hbkg", "hbkg", 100, -1.5, 1.5, 100, -1.5, 1.5, 100, 0, zmax)

x, y, z = ctypes.c_float(), ctypes.c_float(), ctypes.c_float()
fsig = TFile("sig.root", "recreate")
tsig = TTree("tree", "tree")
tsig.SetDirectory(fsig)
tsig.Branch("x", x, "x/F");
tsig.Branch("y", y, "y/F");
tsig.Branch("z", z, "z/F");
fbkg = TFile("bkg.root", "recreate")
tbkg = TTree("tree", "tree")
tbkg.SetDirectory(fbkg)
tbkg.Branch("x", x, "x/F");
tbkg.Branch("y", y, "y/F");
tbkg.Branch("z", z, "z/F");

for i in range(100000):
    z.value = gRandom.Uniform(0, zmax)
    x.value =  cos(z.value)+sin(z.value)+gRandom.Gaus(0, smear)
    y.value = -sin(z.value)+cos(z.value)+gRandom.Gaus(0, smear)
    tsig.Fill()
    hsig.Fill(x.value,y.value,z.value)

    phase = pi
    z.value = gRandom.Uniform(0, zmax)
    x.value =  cos(z.value+phase)+sin(z.value+phase)+gRandom.Gaus(0, smear)
    y.value = -sin(z.value+phase)+cos(z.value+phase)+gRandom.Gaus(0, smear)
    tbkg.Fill()
    hbkg.Fill(x.value,y.value,z.value)

hsig.SetMarkerColor(kBlue)
hbkg.SetMarkerColor(kRed)
c = TCanvas("c", "c", 500, 500)
hsig.Draw()
hbkg.Draw("same")

dummy = raw_input('Press enter to continue')

fsig.Write()
fbkg.Write()
fsig.Close()
fbkg.Close()
###################################################

## Prepare to load ntuple
datasets = {'sig':["sig.root"], "bkg":["bkg.root"]}

## Define variable groups
variables = {
    'ft_0':['x', 'y', 'z'],
}

## Start to setup the network
batch_size = 100
data_iter = RootIter.RootIter('tree', datasets, variables, batch_size)

net = mx.sym.Variable('ft_0')
#net = mx.sym.Variable('im_0')
net = mx.sym.FullyConnected(data=net, name='fc1', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac1', act_type="relu")
net = mx.sym.FullyConnected(data=net, name='fc2', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac2', act_type="relu")
net = mx.sym.FullyConnected(data=net, name='fc3', num_hidden=128)
net = mx.sym.Activation(data=net, name='ac3', act_type="relu")

net = mx.sym.SoftmaxOutput(data=net, name='out')

print(net.list_arguments())
print(net.list_outputs())

import logging
logging.basicConfig(level=logging.INFO)

mod = mx.mod.Module(symbol=net, data_names=variables.keys(), label_names=['out_label'])
mod.fit(data_iter,
        initializer=mx.init.Xavier(magnitude=2.),
        optimizer='sgd', optimizer_params=(('learning_rate', 0.1), ),
        eval_metric='acc',  # report accuracy during training
        batch_end_callback = mx.callback.Speedometer(batch_size, 10),
        num_epoch=5)

