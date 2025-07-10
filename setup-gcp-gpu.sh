#!/bin/bash
# Setup-Skript für Google Cloud GPU-VM
# Ausführung: bash setup-gcp-gpu.sh

set -e  # Bei Fehler stoppen

echo "=== AMP-RBC-MD GPU SETUP FÜR GOOGLE CLOUD ==="

# System-Pakete aktualisieren
echo "Aktualisiere System-Pakete..."
sudo apt-get update
sudo apt-get install -y build-essential cmake git wget curl

# CUDA-Treiber prüfen
echo "Prüfe CUDA-Treiber..."
nvidia-smi || {
    echo "CUDA-Treiber nicht gefunden! Installiere CUDA 11.8..."
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
    sudo dpkg -i cuda-keyring_1.0-1_all.deb
    sudo apt-get update
    sudo apt-get install -y cuda-toolkit-11-8
    echo 'export PATH=/usr/local/cuda-11.8/bin:$PATH' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
    source ~/.bashrc
}

# Miniconda installieren
echo "Installiere Miniconda..."
if [ ! -d "$HOME/miniconda3" ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
fi

# Conda-Umgebung erstellen
echo "Erstelle amp-rbc-md Umgebung..."
conda env create -f environment-linux.yml || {
    echo "Conda-Umgebung existiert bereits, aktualisiere..."
    conda env update -f environment-linux.yml
}

# Umgebung aktivieren
echo "Aktiviere Umgebung..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Projekt installieren
echo "Installiere amp-rbc-md..."
pip install -e .

# CUDA-kompatible JAX installieren
echo "Installiere CUDA-kompatible JAX..."
bash install-jax-cuda.sh

# GPU-Test
echo "Teste GPU-Zugriff..."
python -c "
import torch
print(f'PyTorch Version: {torch.__version__}')
print(f'CUDA verfügbar: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA Version: {torch.version.cuda}')
    print(f'GPU Gerät: {torch.cuda.get_device_name(0)}')
    print(f'GPU Speicher: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')
"

# ColabFold-Test
echo "Teste ColabFold..."
python -c "
import colabfold
print('ColabFold verfügbar:', hasattr(colabfold, 'batch'))
if hasattr(colabfold, 'batch'):
    print('ColabFold Batch-Modul verfügbar!')
    print('ColabFold Version:', colabfold.__version__)
"

echo "=== SETUP ABGESCHLOSSEN ==="
echo "Starten Sie eine Simulation mit:"
echo "amp-rbc-md --seq YYLKRLIKRYKTLKNR --n-replica 1 --profile default -j 1" 