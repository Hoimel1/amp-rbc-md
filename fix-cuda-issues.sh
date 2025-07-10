#!/bin/bash
# Fix-Skript für CUDA- und GROMACS-Probleme
# Ausführung: bash fix-cuda-issues.sh

set -e

echo "=== FIX CUDA UND GROMACS PROBLEME ==="

# Aktiviere conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

echo "1. Installiere gmxapi für GROMACS Python-Interface..."
conda install -c conda-forge gmxapi -y

echo "2. Prüfe CUDA-Version..."
CUDA_VERSION=$(nvcc --version | grep "release" | sed 's/.*release \([0-9]\+\.[0-9]\+\).*/\1/')
echo "System CUDA Version: $CUDA_VERSION"

echo "3. Deinstalliere alte JAX-Version..."
pip uninstall jax jaxlib -y

echo "4. Installiere JAX mit korrekter CUDA-Version..."
if [[ "$CUDA_VERSION" == "12.1" ]]; then
    echo "Installiere JAX für CUDA 12.1..."
    pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
elif [[ "$CUDA_VERSION" == "12.3" ]]; then
    echo "Installiere JAX für CUDA 12.3..."
    pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
else
    echo "Installiere JAX für CUDA 12.x..."
    pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
fi

echo "5. Setze CUDA-Umgebungsvariablen..."
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

echo "6. Teste Installation..."
python verify-installation.py

echo "=== FIX ABGESCHLOSSEN ===" 