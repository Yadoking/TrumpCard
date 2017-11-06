#!/bin/bash

if [ -d sw ]; then
    echo "Basic software setup is already done. Skip"
    exit 1
fi

mkdir -p sw/install
virtualenv --system-site-packages ./sw/python27

source sw/python27/bin/activate

pip install --upgrade pip
pip install --upgrade numpy protobuf wheel 
if [ -f /share/apps/tensorflow-1.3.0-cp27-cp27mu-linux_x86_64.whl ]; then
    pip install /share/apps/tensorflow-1.3.0-cp27-cp27mu-linux_x86_64.whl 
fi

pip install h5py
pip install graphviz
pip install pydot keras 

pip install -U scikit-learn

