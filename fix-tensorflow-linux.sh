#!/bin/bash
# Installiert TensorFlow für Linux mit CUDA-Unterstützung
# Ausführung: bash fix-tensorflow-linux.sh

set -e

echo "=== INSTALLIERE TENSORFLOW FÜR LINUX ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Entferne alte TensorFlow-Installation
echo "Entferne alte TensorFlow-Installation..."
pip uninstall -y tensorflow tensorflow-macos

# Installiere TensorFlow für Linux
echo "Installiere TensorFlow 2.13.1 für Linux..."
pip install tensorflow==2.13.1

# Teste TensorFlow-CUDA
echo "Teste TensorFlow-CUDA..."
python -c "
import tensorflow as tf
print('=== TENSORFLOW-CUDA TEST ===')
print(f'TensorFlow Version: {tf.__version__}')
print(f'CUDA verfügbar: {len(tf.config.list_physical_devices(\"GPU\")) > 0}')
if len(tf.config.list_physical_devices(\"GPU\")) > 0:
    print(f'GPU Geräte: {tf.config.list_physical_devices(\"GPU\")}')
    
    # Teste GPU-Berechnung
    with tf.device('/GPU:0'):
        a = tf.constant([[1.0, 2.0], [3.0, 4.0]])
        b = tf.constant([[1.0, 1.0], [0.0, 1.0]])
        c = tf.matmul(a, b)
        print('GPU-Berechnung erfolgreich!')
        print(f'Ergebnis: {c.numpy()}')
else:
    print('CUDA nicht verfügbar!')
"

echo "=== TENSORFLOW FÜR LINUX INSTALLIERT ===" 