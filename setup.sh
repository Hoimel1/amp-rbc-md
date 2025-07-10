#!/bin/bash
# amp-rbc-md Setup Script für Linux
# PyTorch 2.3.0+cu121 + nvidia-cudnn-cu12==8.9.2.26 + JAX 0.4.25
# Ausführung: bash setup.sh

set -e

echo "=== AMP-RBC-MD LINUX SETUP ==="

# Prüfe CUDA-Installation
if ! command -v nvidia-smi &> /dev/null; then
    echo "❌ NVIDIA GPU nicht gefunden. Bitte installieren Sie NVIDIA-Treiber."
    exit 1
fi

echo "✅ NVIDIA GPU gefunden:"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader

# Aktiviere conda
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
elif [ -f ~/anaconda3/etc/profile.d/conda.sh ]; then
    source ~/anaconda3/etc/profile.d/conda.sh
else
    echo "❌ Conda nicht gefunden. Bitte installieren Sie Miniconda oder Anaconda."
    exit 1
fi

# Entferne alte Umgebung falls vorhanden
echo "Entferne alte Umgebung..."
conda deactivate
conda env remove -n amp-rbc-md -y 2>/dev/null || true

# Erstelle neue Umgebung
echo "Erstelle neue Umgebung..."
conda env create -f environment.yml

# Aktiviere neue Umgebung
conda activate amp-rbc-md

# Installiere JAX 0.4.25 mit CUDA-Unterstützung
echo "Installiere JAX 0.4.25 mit CUDA..."
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Installiere PyTorch 2.3.0+cu121
echo "Installiere PyTorch 2.3.0+cu121..."
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Installiere spezifische cudnn-Version
echo "Installiere nvidia-cudnn-cu12==8.9.2.26..."
pip install nvidia-cudnn-cu12==8.9.2.26

# Installiere ColabFold OHNE Abhängigkeiten zu überschreiben
echo "Installiere ColabFold ohne Abhängigkeiten zu überschreiben..."
pip install colabfold==1.5.5 --no-deps

# Installiere fehlende ColabFold-Abhängigkeiten
echo "Installiere ColabFold-Abhängigkeiten..."
pip install absl-py appdirs biopython matplotlib numpy pandas py3Dmol requests tqdm dm-haiku

# Installiere Projekt
echo "Installiere amp-rbc-md..."
pip install -e .

# Teste Installation
echo "Teste Installation..."
python verify-installation.py

echo "=== SETUP ERFOLGREICH ABGESCHLOSSEN ==="
echo ""
echo "Nächste Schritte:"
echo "1. Testen Sie eine Simulation:"
echo "   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run"
echo "2. Führen Sie eine echte Simulation aus:"
echo "   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1" 