#!/bin/bash
# Behebt Abhängigkeitskonflikte für Google Cloud VM
# Ausführung: bash fix-dependencies.sh

set -e

echo "=== BEHEBE ABHÄNGIGKEITSKONFLIKTE ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Entferne problematische Pakete
echo "Entferne problematische Pakete..."
pip uninstall -y numpy pandas tensorflow-macos torch

# Installiere kompatible Versionen
echo "Installiere kompatible Versionen..."
conda install -y numpy=1.24.3 pandas=1.5.3

# Installiere kompatible PyTorch-Version
echo "Installiere kompatible PyTorch-Version..."
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# JAX-CUDA installieren
echo "Installiere JAX-CUDA..."
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# ColabFold neu installieren
echo "Installiere ColabFold neu..."
pip install --force-reinstall colabfold

# Teste Installation
echo "Teste Installation..."
python -c "
import jax
import jaxlib
import colabfold
import pandas
import numpy
import torch

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
"

echo "=== ABHÄNGIGKEITEN BEHOBEN ===" 