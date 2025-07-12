#!/usr/bin/env bash
set -e

echo "🔧 Installing FastFold without CUDA extensions..."

# Aktiviere die Folding-Umgebung
conda activate folding-env

# Setze Umgebungsvariablen für FastFold ohne CUDA-Kernel
export CUDA_EXT=0
export FASTFOLD_SKIP_CUDA=1
export TORCH_CUDA_ARCH_LIST=""

echo "📝 Installing FastFold with CUDA extensions disabled..."

# Installiere FastFold ohne CUDA-Kernel
pip install git+https://github.com/hpcaitech/FastFold@0.1.0#egg=fastfold \
    --no-build-isolation \
    --config-settings build_ext.define="CUDA_EXT=0" \
    --config-settings build_ext.define="FASTFOLD_SKIP_CUDA=1"

echo "✅ FastFold installation completed!"
echo "📝 Note: FastFold is installed without CUDA extensions but with PyTorch CUDA support"
echo "🔍 Testing FastFold installation..."
python -c "import fastfold; print('✅ FastFold imported successfully!')" 