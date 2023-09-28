#!/bin/bash
echo "Creating ExoCTK conda environment for Python $PYTHON_VERSION"
conda env update -f "env/environment-${PYTHON_VERSION}.yml" || exit 1
export CONDA_ENV=exoctk-$PYTHON_VERSION
source activate $CONDA_ENV

echo "The installed environment:"
conda env export