#!/bin/bash
# Behebt PyTorch-CUDA-Abhängigkeitskonflikte
# Ausführung: bash fix-pytorch-cuda.sh

set -e

echo "=== BEHEBE PYTORCH-CUDA KONFLIKTE ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Entferne alte PyTorch-Installation
echo "Entferne alte PyTorch-Installation..."
pip uninstall -y torch torchvision torchaudio

# Installiere kompatible PyTorch-Version mit CUDA 12.1
echo "Installiere PyTorch 2.3.0 mit CUDA 12.1..."
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Teste PyTorch-CUDA
echo "Teste PyTorch-CUDA..."
python -c "
import torch
print('=== PYTORCH-CUDA TEST ===')
print(f'PyTorch Version: {torch.__version__}')
print(f'CUDA verfügbar: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA Version: {torch.version.cuda}')
    print(f'GPU Gerät: {torch.cuda.get_device_name(0)}')
    print(f'GPU Speicher: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')
    
    # Teste GPU-Berechnung
    x = torch.randn(1000, 1000).cuda()
    y = torch.mm(x, x.t())
    print('GPU-Berechnung erfolgreich!')
else:
    print('CUDA nicht verfügbar!')
"

echo "=== PYTORCH-CUDA KONFLIKTE BEHOBEN ===" 