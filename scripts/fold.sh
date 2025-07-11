#!/usr/bin/env bash
set -e
conda activate folding_env
export FASTFOLD_SKIP_TEMPLATES=1
export FASTFOLD_NO_MSA=1
export FASTFOLD_PARAMS_PATH=$HOME/.fastfold/fastfold_params
python ~/fastfold/inference.py "$@" 