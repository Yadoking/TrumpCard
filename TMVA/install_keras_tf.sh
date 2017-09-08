#!/bin/bash

if [ !-d python27 ];
  virtualenv --system-site-packages ./python27
fi
source python27/bin/activate

pip install --update pip
pip install --ugrade pip
pip install --upgrade pip
pip install --upgrade numpy protobuf wheel 
pip install /share/apps/tensorflow-1.3.0-cp27-cp27mu-linux_x86_64.whl 

pip install h5py
pip install graphviz
pip install pydot keras 

pip install -U scikit-learn

