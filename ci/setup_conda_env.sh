#!/bin/bash
echo "Creating ExoCTK conda environment"
conda env update -f "exoctk-env.yml" || exit 1
export CONDA_ENV=exoctk-env
source activate $CONDA_ENV

echo "The installed environment:"
conda env export
