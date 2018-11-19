#!/bin/bash

echo "Creating conda environment"
conda env create -f environment.yml || exit 1
source activate exoctk