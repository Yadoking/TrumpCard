#!/bin/bash

if [ -d sw ]; then
    echo "Basic software setup is already done. Skip"
    source sw/python27/bin/activate
else
    mkdir -p sw/install
    virtualenv --system-site-packages ./sw/python27

    source sw/python27/bin/activate

    pip install --upgrade pip

    pip install --upgrade numpy protobuf wheel 
    pip install h5py
    pip install graphviz

    ## Install TensorFlow
    if [ -f /share/apps/tensorflow-1.3.0-cp27-cp27mu-linux_x86_64.whl ]; then
        pip install /share/apps/tensorflow-1.3.0-cp27-cp27mu-linux_x86_64.whl 
    fi

    ## Install mxnet
    cd sw/install
    git clone https://github.com/dmlc/mxnet --branch 0.12.0 --recursive
    cd mxnet
    cp make/config.mk ./
    if [ -d /usr/local/cuda ]; then
        echo 'USE_CUDA = 1' >> config.mk
        echo 'USE_CUDA_PATH = /usr/local/cuda' >> config.mk
        echo 'USE_CUDNN = 1' >> config.mk
    fi
    make -j $(nproc)
    cd python
    python setup.py install

    cd ../../../..

    pip install pydot keras 
    pip install -U scikit-learn
fi
