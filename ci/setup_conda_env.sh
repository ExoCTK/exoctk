#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n exoctk python=$PYTHON_VERSION || exit 1
source activate exoctk

echo "Installing packages..."
conda install numpy astropy scipy cython matplotlib numba mock bokeh h5py sphinx pandas lmfit
pip install numpy astropy bibtexparser astroquery svo_filters
pip install batman-package