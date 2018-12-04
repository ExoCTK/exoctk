#!/bin/bash

echo "Creating conda environment for $PYTHON_VERSION"
conda env create -f "environment-${PYTHON_VERSION}.yml" || exit 1
source activate "exoctk-${PYTHON_VERSION}"