#!/bin/bash
# Komplette Neuinstallation der amp-rbc-md Umgebung ohne Konflikte
# Ausführung: bash reset-environment.sh

set -e

echo "=== KOMPLETTE UMWELT-NEUINSTALLATION ==="

# Aktiviere conda
source ~/miniconda3/etc/profile.d/conda.sh

# Entferne alte Umgebung
echo "Entferne alte Umgebung..."
conda deactivate
conda env remove -n amp-rbc-md -y

# Erstelle neue Umgebung mit minimalen Abhängigkeiten
echo "Erstelle neue Umgebung..."
conda create -n amp-rbc-md python=3.10 -y

# Aktiviere neue Umgebung
conda activate amp-rbc-md

# Installiere Basis-Pakete
echo "Installiere Basis-Pakete..."
conda install -c conda-forge -c bioconda -y gromacs=2024 biopython click pyyaml tqdm matplotlib mlflow moviepy rich pytest pytest-cov black flake8 mypy isort

# Installiere kompatible NumPy/Pandas
echo "Installiere kompatible NumPy/Pandas..."
conda install -y numpy=1.24.3 pandas=1.5.3

# Installiere JAX-CUDA (vor ColabFold)
echo "Installiere JAX-CUDA..."
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Prüfe JAX-CUDA Installation
echo "Prüfe JAX-CUDA Installation..."
python -c "
import jax
import jaxlib
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')
if len([d for d in jax.devices() if d.platform == \"gpu\"]) == 0:
    print('WARNUNG: CUDA nicht erkannt! Installiere CUDA-kompatible Version...')
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--force-reinstall', 'jax==0.4.25', 'jaxlib==0.4.25+cuda12.cudnn89', '-f', 'https://storage.googleapis.com/jax-releases/jax_cuda_releases.html'])
    print('JAX-CUDA neu installiert!')
"

# Installiere ColabFold (ohne TensorFlow-Konflikte)
echo "Installiere ColabFold..."
pip install colabfold==1.5.5

# Installiere PyTorch separat
echo "Installiere PyTorch..."
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

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
print(f'JAX CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')
print(f'PyTorch Version: {torch.__version__}')
print(f'PyTorch CUDA verfügbar: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'PyTorch CUDA Version: {torch.version.cuda}')
    print(f'PyTorch GPU: {torch.cuda.get_device_name(0)}')
print(f'ColabFold verfügbar: {hasattr(colabfold, \"batch\")}')
print(f'Pandas Version: {pandas.__version__}')
print(f'NumPy Version: {numpy.__version__}')
print('=== ALLE TESTS ERFOLGREICH ===')
"

echo "=== UMWELT ERFOLGREICH NEU INSTALLIERT ===" 