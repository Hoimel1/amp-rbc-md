#!/usr/bin/env bash
set -e
conda activate md_env
amp-rbc-md run_sim --pdb "$1" --n-replica 1 -j 1 