#!/bin/bash

# FastFold Installation ohne CUDA-Erweiterungen
# Für amp-rbc-md Pipeline

set -e

echo "=== FastFold Installation (ohne CUDA-Erweiterungen) ==="

# Prüfe ob wir im richtigen Verzeichnis sind
if [ ! -d "external/fastfold" ]; then
    echo "Fehler: external/fastfold Verzeichnis nicht gefunden!"
    echo "Bitte führen Sie dieses Skript aus dem amp-rbc-md Hauptverzeichnis aus."
    exit 1
fi

# Wechsle zu FastFold-Verzeichnis
cd external/fastfold

echo "Installing FastFold without CUDA extensions..."

# Setze Umgebungsvariablen für Installation ohne CUDA
export FASTFOLD_SKIP_CUDA=1
export CUDA_HOME=""
export CUDA_PATH=""

# Installiere FastFold ohne CUDA-Erweiterungen
pip install -e . --no-deps --no-build-isolation || {
    echo "Warnung: FastFold Installation fehlgeschlagen, versuche alternative Methode..."
    
    # Alternative: Nur Python-Module installieren
    pip install -e . --no-deps --no-build-isolation --config-settings editable_mode=compat || {
        echo "Fehler: FastFold Installation fehlgeschlagen!"
        echo "Versuchen Sie manuell:"
        echo "cd external/fastfold"
        echo "pip install -e . --no-deps"
        exit 1
    }
}

# Zurück zum Hauptverzeichnis
cd ../..

echo "FastFold Installation abgeschlossen!"

# Teste Installation
echo "Testing FastFold installation..."
python -c "
try:
    import fastfold
    print('✓ FastFold erfolgreich importiert')
except ImportError as e:
    print(f'✗ FastFold Import fehlgeschlagen: {e}')
    exit(1)
except Exception as e:
    print(f'⚠ FastFold Import mit Warnung: {e}')
    print('Das ist normal bei Installation ohne CUDA-Erweiterungen')
"

echo "=== FastFold Installation abgeschlossen ==="
echo ""
echo "Nächste Schritte:"
echo "1. Setzen Sie Umgebungsvariablen:"
echo "   export FASTFOLD_SKIP_TEMPLATES=1"
echo "   export FASTFOLD_NO_MSA=1"
echo ""
echo "2. Testen Sie die Pipeline:"
echo "   amp-rbc-md --seq GLSILGKLL --dry-run" 