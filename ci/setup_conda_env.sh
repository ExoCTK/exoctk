#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n exoctk python=$PYTHON_VERSION || exit 1
source activate exoctk
conda env list

echo "Installing packages..."
pip install numpy astropy
conda install numpy astropy scipy cython matplotlib numba mock bokeh h5py sphinx pandas flask
pip install bibtexparser astroquery svo_filters==0.2.16 batman-package lmfit
echo $PATH