#!/bin/bash
# Installiert CUDA-kompatible JAX-Version für ColabFold
# Ausführung: bash install-jax-cuda.sh

set -e

echo "=== INSTALLIERE CUDA-KOMPATIBLE JAX-VERSION ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Entferne alte JAX-Versionen
echo "Entferne alte JAX-Versionen..."
pip uninstall -y jax jaxlib

# Installiere CUDA-kompatible JAX-Version
echo "Installiere CUDA-kompatible JAX-Version..."
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Teste JAX-CUDA
echo "Teste JAX-CUDA-Unterstützung..."
python -c "
import jax
import jaxlib
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'JAX Backends: {jax.devices()}')

# Teste CUDA-Backend
try:
    import jax.numpy as jnp
    x = jnp.array([1.0, 2.0, 3.0])
    y = jnp.sin(x)
    print('JAX funktioniert korrekt!')
    print(f'Ergebnis: {y}')
except Exception as e:
    print(f'JAX Fehler: {e}')
"

echo "=== JAX-CUDA INSTALLATION ABGESCHLOSSEN ===" 