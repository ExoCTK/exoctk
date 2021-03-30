#!/bin/bash
echo "Creating base conda environment for Python $PYTHON_VERSION"
conda create --yes --prefix /home/travis/envs python=$PYTHON_VERSION
conda activate /home/travis/envs
conda install sphinx numpy flask

echo "Creating ExoCTK conda environment for Python $PYTHON_VERSION"
conda env update -f "env/environment-${PYTHON_VERSION}.yml" || exit 1
export CONDA_ENV=exoctk-$PYTHON_VERSION
source activate $CONDA_ENV

echo "The installed environment:"
conda env export