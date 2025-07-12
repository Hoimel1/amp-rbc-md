#!/usr/bin/env bash
set -e

echo "🔧 Installing FastFold..."

# Aktiviere die Folding-Umgebung
conda activate folding-env

# Einfache FastFold-Installation
pip install git+https://github.com/hpcaitech/FastFold@0.1.0

echo "✅ FastFold installation completed!"
echo "🔍 Testing FastFold installation..."
python -c "import fastfold; print('✅ FastFold imported successfully!')" 