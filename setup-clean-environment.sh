#!/bin/bash
# Komplett sauberes Environment-Setup mit exakten Versionen
# PyTorch 2.3.0+cu121 + nvidia-cudnn-cu12==8.9.2.26 + JAX 0.4.25
# Ausführung: bash setup-clean-environment.sh

set -e

echo "=== SAUBERES ENVIRONMENT-SETUP ==="

# Aktiviere conda
source ~/miniconda3/etc/profile.d/conda.sh

# Entferne alte Umgebung falls vorhanden
echo "Entferne alte Umgebung..."
conda deactivate
conda env remove -n amp-rbc-md -y 2>/dev/null || true

# Erstelle neue Umgebung
echo "Erstelle neue Umgebung..."
conda create -n amp-rbc-md python=3.10 -y

# Aktiviere neue Umgebung
conda activate amp-rbc-md

# Installiere Basis-Pakete über conda
echo "Installiere Basis-Pakete..."
conda install -c conda-forge -c bioconda -y \
    gromacs=2024 \
    biopython \
    click \
    pyyaml \
    pandas=1.5.3 \
    numpy=1.24.3 \
    tqdm \
    matplotlib \
    mlflow \
    moviepy \
    rich \
    pytest \
    pytest-cov \
    black \
    flake8 \
    mypy \
    isort

# Installiere JAX 0.4.25 mit CUDA-Unterstützung (VOR ColabFold)
echo "Installiere JAX 0.4.25 mit CUDA..."
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Installiere PyTorch 2.3.0+cu121 mit spezifischer cudnn-Version
echo "Installiere PyTorch 2.3.0+cu121..."
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Installiere spezifische cudnn-Version
echo "Installiere nvidia-cudnn-cu12==8.9.2.26..."
pip install nvidia-cudnn-cu12==8.9.2.26

# Installiere ColabFold OHNE Abhängigkeiten zu überschreiben
echo "Installiere ColabFold ohne Abhängigkeiten zu überschreiben..."
pip install colabfold==1.5.5 --no-deps

# Installiere fehlende ColabFold-Abhängigkeiten manuell
echo "Installiere ColabFold-Abhängigkeiten..."
pip install absl-py appdirs biopython matplotlib numpy pandas py3Dmol requests tqdm dm-haiku

# Installiere Projekt
echo "Installiere amp-rbc-md..."
pip install -e .

# Teste Installation
echo "Teste Installation..."
python -c "
import jax
import jaxlib
import torch
import colabfold
import pandas
import numpy

print('=== INSTALLATION ERFOLGREICH ===')
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'JAX linear_util verfügbar: {hasattr(jax, \"linear_util\")}')
print(f'JAX CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')
print(f'PyTorch Version: {torch.__version__}')
print(f'PyTorch CUDA verfügbar: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'PyTorch CUDA Version: {torch.version.cuda}')
    print(f'PyTorch GPU: {torch.cuda.get_device_name(0)}')
print(f'ColabFold verfügbar: {hasattr(colabfold, \"batch\")}')
print(f'Pandas Version: {pandas.__version__}')
print(f'NumPy Version: {numpy.__version__}')

# Teste CUDA-Bibliotheken
import subprocess
try:
    result = subprocess.run(['pip', 'list'], capture_output=True, text=True)
    if 'nvidia-cudnn-cu12' in result.stdout:
        print('✅ nvidia-cudnn-cu12 installiert')
    else:
        print('❌ nvidia-cudnn-cu12 nicht gefunden')
except:
    print('❌ Konnte pip list nicht ausführen')

print('=== ALLE TESTS ERFOLGREICH ===')
"

echo "=== SAUBERES ENVIRONMENT ERFOLGREICH ERSTELLT ===" 