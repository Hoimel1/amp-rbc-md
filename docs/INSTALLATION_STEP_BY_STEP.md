# Schritt-für-Schritt Installation für amp-rbc-md

Diese Anleitung führt Sie durch die manuelle Installation von amp-rbc-md auf Ubuntu 22.04 mit NVIDIA L4 GPU.

## Voraussetzungen

- Ubuntu 22.04
- NVIDIA L4 GPU mit CUDA 12.1
- 12 CPU-Kerne, 48 GB RAM
- 400 GB Speicherplatz

## Schritt 1: System vorbereiten

```bash
# System aktualisieren
sudo apt update && sudo apt upgrade -y

# NVIDIA-Treiber prüfen
nvidia-smi

# CUDA-Version prüfen
nvcc --version
```

## Schritt 2: Miniconda installieren

```bash
# Miniconda herunterladen und installieren
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# Conda initialisieren
$HOME/miniconda3/bin/conda init zsh
source ~/.zshrc

# Conda-Version prüfen
conda --version
```

## Schritt 3: Repository klonen

```bash
# Repository mit Submodulen klonen
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# Submodule-Status prüfen
git submodule status
```

## Schritt 4: Conda-Umgebung erstellen

```bash
# Umgebung aus environment.yml erstellen
conda env create -f environment.yml

# Umgebung aktivieren
conda activate amp-rbc-md

# Installation prüfen
python --version
conda list | grep -E "(python|numpy|pandas)"
```

## Schritt 5: GROMACS installieren

```bash
# GROMACS aus conda-forge installieren (GPU-fähig)
conda install -c conda-forge gromacs=2024.5

# GROMACS-Version prüfen
gmx --version
```

## Schritt 6: PyTorch installieren

```bash
# PyTorch mit CUDA 12.1 aus conda-forge installieren
conda install -c conda-forge pytorch=2.5.1 pytorch-cuda=12.1

# PyTorch-Test
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"
```

## Schritt 7: ColabFold installieren

```bash
# ColabFold mit AlphaFold-Unterstützung installieren
pip install "colabfold[alphafold]"

# ColabFold-Test
python -c "import colabfold; print('ColabFold installiert')"
```

## Schritt 8: FastFold installieren

```bash
# FastFold aus dem Submodul installieren
cd external/fastfold
pip install -e . --no-build-isolation
cd ../..

# FastFold-Test
python -c "import fastfold; print('FastFold installiert')"
```

## Schritt 9: Datenbanken herunterladen (Lightweight)

```bash
# Lightweight-Datenbanken herunterladen (~300 GB)
mkdir -p $HOME/alphafold_data
export OPENFOLD_DATA=$HOME/alphafold_data

# Dummy-Verzeichnisse für Template-Skip erstellen
mkdir -p $OPENFOLD_DATA/pdb_mmcif/mmcif_files

# ColabFold-Parameter herunterladen (~3.6 GB)
mkdir -p $HOME/.colabfold
cd $HOME/.colabfold
wget https://github.com/sokrypton/ColabFold/releases/download/v1.5.5/params.tar.gz
tar -xzf params.tar.gz
cd -

# FastFold-Parameter herunterladen (~3.6 GB)
mkdir -p $HOME/.fastfold
cd $HOME/.fastfold
wget https://github.com/hpcaitech/FastFold/releases/download/v0.2.0/fastfold_params.tar.gz
tar -xzf fastfold_params.tar.gz
cd -
```

## Schritt 10: Pipeline installieren

```bash
# Pipeline im Entwicklungsmodus installieren
pip install -e .

# Installation testen
amp-rbc-md --help
```

## Schritt 11: Testen

```bash
# Dry-Run mit ColabFold-Batch
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run

# Dry-Run mit FastFold (ohne MSA/Templates)
export FASTFOLD_SKIP_TEMPLATES=1
export FASTFOLD_NO_MSA=1
amp-rbc-md --seq GLSILGKLL --dry-run
```

## Schritt 12: Erste Simulation

```bash
# Einzel-Simulation mit ColabFold
amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1 --profile default -j 1

# Batch-Simulation
amp-rbc-md -f examples/batch.fasta --backend colabfold --n-replica 1 --profile default -j 3
```

## Troubleshooting

### FastFold-Installation schlägt fehl
```bash
# Alternative: FastFold mit anderen Flags
cd external/fastfold
pip install -e . --no-build-isolation --config-settings editable_mode=compat
cd ../..
```

### GROMACS-Probleme
```bash
# GROMACS neu installieren
conda remove gromacs
conda install -c conda-forge gromacs=2024.5
```

### PyTorch-Probleme
```bash
# PyTorch neu installieren
conda remove pytorch pytorch-cuda
conda install -c conda-forge pytorch=2.5.1 pytorch-cuda=12.1
```

### ColabFold-Probleme
```bash
# ColabFold neu installieren
pip uninstall colabfold
pip install "colabfold[alphafold]"
```

### Speicherplatz-Probleme
```bash
# Speicherplatz prüfen
df -h

# Temporäre Dateien löschen
conda clean -a
pip cache purge
```

## Nächste Schritte

1. **Erste Simulation**: Führen Sie eine einfache Simulation aus
2. **Batch-Verarbeitung**: Testen Sie Batch-Simulationen
3. **Analyse**: Verwenden Sie die Analyse-Tools
4. **HPC**: Konfigurieren Sie für Slurm-Cluster

## Support

Bei Problemen:
1. Prüfen Sie die Logs in `results/`
2. Verwenden Sie `--dry-run` für Debugging
3. Kontaktieren Sie das Team über GitHub Issues 