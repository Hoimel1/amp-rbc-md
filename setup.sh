#!/bin/bash
# AMP-RBC-MD Setup Script
# Vollständige Installation für Ubuntu 22.04 mit GPU-Support

set -e  # Exit on error

echo "=== AMP-RBC-MD Setup ==="
echo "Dieses Skript installiert alle Abhängigkeiten für AMP-RBC-MD"

# Prüfe ob wir in einem Git-Repository sind
if [ ! -d ".git" ]; then
    echo "❌ Fehler: Bitte führe dieses Skript im AMP-RBC-MD Repository aus"
    exit 1
fi

# Prüfe ob conda verfügbar ist
if ! command -v conda &> /dev/null; then
    echo "❌ Fehler: Conda ist nicht installiert. Bitte installiere Miniconda/Anaconda zuerst."
    echo "Download: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Prüfe GPU
if command -v nvidia-smi &> /dev/null; then
    echo "✅ NVIDIA GPU gefunden"
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader,nounits
else
    echo "⚠️  Keine NVIDIA GPU gefunden - CPU-only Modus"
fi

# Git-Submodule initialisieren
echo "📥 Initialisiere Git-Submodule..."
git submodule update --init --recursive

# Conda-Environment erstellen
echo "🐍 Erstelle Conda-Environment..."
if conda env list | grep -q "amp-rbc-md"; then
    echo "⚠️  Environment 'amp-rbc-md' existiert bereits. Überspringe..."
else
    conda env create -f environment.yml
fi

# Environment aktivieren
echo "🔧 Aktiviere Environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate amp-rbc-md

# FastFold installieren (editable)
echo "⚡ Installiere FastFold..."
cd external/fastfold
pip install -e . --no-build-isolation
cd ../..

# OpenFold installieren (editable)
echo "🧬 Installiere OpenFold..."
cd external/openfold
pip install -e . --no-build-isolation
cd ../..

# Verzeichnisse erstellen
echo "📁 Erstelle Verzeichnisse..."
mkdir -p $HOME/alphafold_dbs
mkdir -p $HOME/openfold_data

# Umgebungsvariablen setzen
echo "🔧 Setze Umgebungsvariablen..."
echo 'export OPENFOLD_DATA=$HOME/openfold_data' >> ~/.bashrc
echo 'export ALPHAFOLD_DATA=$HOME/alphafold_dbs' >> ~/.bashrc

# Installation testen
echo "🧪 Teste Installation..."
python verify-installation.py --engine fastfold

echo ""
echo "✅ Setup abgeschlossen!"
echo ""
echo "=== NÄCHSTE SCHRITTE ==="
echo "1. Starte eine neue Shell oder führe aus: source ~/.bashrc"
echo "2. Lade die Datenbanken herunter:"
echo "   cd external/fastfold"
echo "   ./scripts/download_all_data.sh \$HOME/alphafold_dbs/"
echo "3. Teste eine Simulation:"
echo "   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run"
echo ""
echo "💡 Tipp: Verwende 'tmux' für lange Downloads im Hintergrund" 