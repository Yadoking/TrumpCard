#!/usr/bin/env python

from ROOT import *
import mxnet as mx
import numpy as np

## Implementation of root iterator
class RootIter(mx.io.DataIter):
    def __init__(self, treeName, datasets, variables, batch_size):
        self.batch_size = batch_size
        self.cur_batch = -1

        self.im_size = {}
        self._provide_data = []
        self._provide_label = []

        self.datasets = {}
        for name in datasets:
            self.datasets[name] = {}
            self.datasets[name]['chain'] = TChain(treeName)
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

