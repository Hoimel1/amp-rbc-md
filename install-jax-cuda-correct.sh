#!/bin/bash
# Installiert die korrekte CUDA-kompatible JAX-Version
# Ausführung: bash install-jax-cuda-correct.sh

set -e

echo "=== INSTALLIERE KORREKTE JAX-CUDA VERSION ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Entferne alte JAX-Versionen
echo "Entferne alte JAX-Versionen..."
pip uninstall -y jax jaxlib

# Installiere CUDA-kompatible JAX-Version (0.4.25)
echo "Installiere JAX 0.4.25 mit CUDA-Unterstützung..."
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Teste JAX-CUDA
echo "Teste JAX-CUDA..."
python -c "
import jax
import jaxlib
print('=== JAX-CUDA TEST ===')
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')

if len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0:
    print('✅ JAX-CUDA funktioniert korrekt!')
    
    # Teste GPU-Berechnung
    import jax.numpy as jnp
    x = jnp.array([1.0, 2.0, 3.0])
    y = jnp.sin(x)
    print('✅ GPU-Berechnung erfolgreich!')
    print(f'Ergebnis: {y}')
else:
    print('❌ CUDA nicht verfügbar!')
    print('Mögliche Ursachen:')
    print('- NVIDIA-Treiber nicht installiert')
    print('- CUDA-Toolkit nicht installiert')
    print('- Falsche JAX-Version installiert')
"

echo "=== JAX-CUDA INSTALLATION ABGESCHLOSSEN ===" 