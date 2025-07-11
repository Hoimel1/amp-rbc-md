# AMP-RBC-MD (Clean 2025)

Eine zweistufige Pipeline zur Simulation antimikrobieller Peptide in RBC-Membranen.

## Architektur
```
folding_env (Py3.9) ──► FastFold  →  *.pdb
                                   ▼
md_env (Py3.10)     ──► amp-rbc-md  →  CG-MD-Simulation
```

## Installation (Kurzfassung)
```bash
# Folding-Env
conda env create -f envs/folding.yml
conda activate folding_env
mkdir -p ~/.fastfold && cd ~/.fastfold
wget -q https://github.com/hpcaitech/FastFold/releases/download/v0.2.0/fastfold_params.tar.gz
 tar -xf fastfold_params.tar.gz ; cd -

# MD-Env
conda env create -f envs/md.yml
```

## Verwendung
```bash
# 1. Falten
conda activate folding_env
scripts/fold.sh seq.fasta ./fold_out --gpus 1

# 2. Simulation
conda activate md_env
scripts/simulate.sh ./fold_out/model_1.pdb
```

## Abhängigkeiten
* Folding-Env: FastFold, PyTorch, CUDA 11.8 
* MD-Env: GROMACS 2024.5, martinize2, numpy-Stack

## Lizenz
MIT
