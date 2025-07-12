#!/usr/bin/env bash
set -e

echo "🔧 Installing FastFold without CUDA extensions..."

# Aktiviere die Folding-Umgebung
conda activate folding-env

# Installiere FastFold ohne CUDA-Kernel
pip install git+https://github.com/hpcaitech/FastFold@0.1.0#egg=fastfold \
    --no-build-isolation \
    --config-settings build_ext.define="CUDA_EXT=0"

echo "✅ FastFold installation completed!"
echo "📝 Note: FastFold is installed without CUDA extensions but with PyTorch CUDA support" 