#!/bin/bash
# Behebt importlib-metadata Abhängigkeitskonflikt
# Ausführung: bash fix-importlib-metadata.sh

set -e

echo "=== BEHEBE IMPORTLIB-METADATA KONFLIKT ==="

# Aktiviere conda-Umgebung
source ~/miniconda3/etc/profile.d/conda.sh
conda activate amp-rbc-md

# Prüfe aktuelle Version
echo "Aktuelle importlib-metadata Version:"
pip show importlib-metadata

# Behebe Konflikt
echo "Behebe importlib-metadata Konflikt..."
pip install --upgrade "importlib-metadata>=6.0,<8.8.0"

# Teste Installation
echo "Teste Installation..."
python -c "
import importlib_metadata
print('=== IMPORTLIB-METADATA TEST ===')
print(f'Version: {importlib_metadata.__version__}')
print('Import erfolgreich!')

# Teste opentelemetry-api
try:
    import opentelemetry.api
    print('opentelemetry-api funktioniert!')
except ImportError as e:
    print(f'opentelemetry-api Fehler: {e}')
"

echo "=== IMPORTLIB-METADATA KONFLIKT BEHOBEN ===" 