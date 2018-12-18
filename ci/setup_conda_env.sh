#!/bin/bash

echo "Creating conda environment for $PYTHON_VERSION"
conda env create -f "environment-${PYTHON_VERSION}.yml" || exit 1
export CONDA_ENV=exoctk-$PYTHON_VERSION
source activate $CONDA_ENV

git clone https://github.com/lkreidberg/batman.git
python batman/setup.py install