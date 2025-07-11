#!/bin/bash

# amp-rbc-md Lightweight Setup für 400GB VM
# Verwendet ColabFold-Batch ohne lokale MSA-Datenbanken

set -e

echo "🚀 amp-rbc-md Lightweight Setup (400GB VM)"
echo "=========================================="

# Farben für Output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Funktionen
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Prüfe verfügbaren Speicherplatz
check_disk_space() {
    log_info "Prüfe verfügbaren Speicherplatz..."
    available_space=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if [ "$available_space" -lt 50 ]; then
        log_error "Nicht genug Speicherplatz! Benötigt: 50GB, Verfügbar: ${available_space}GB"
        exit 1
    fi
    
    log_success "Speicherplatz OK: ${available_space}GB verfügbar"
}

# Conda Environment erstellen
setup_conda() {
    log_info "Erstelle Conda Environment..."
    
    if command -v conda &> /dev/null; then
        conda env create -f environment.yml
        log_success "Conda Environment erstellt"
    else
        log_error "Conda nicht gefunden. Bitte installiere Miniconda/Anaconda."
        exit 1
    fi
}

# ColabFold installieren
setup_colabfold() {
    log_info "Installiere ColabFold-Batch (ohne lokale DBs)..."
    
    # Aktiviere Environment
    source ~/.bashrc
    conda activate amp-rbc-md
    
    # Installiere ColabFold
    pip install colabfold batchfold
    
    log_success "ColabFold installiert"
}

# Nur AlphaFold-Parameter herunterladen
download_alphafold_params() {
    log_info "Lade AlphaFold-Parameter herunter (~3.6 GB)..."
    
    # Erstelle Download-Verzeichnis
    mkdir -p ~/alphafold_dbs
    
    # Wechsle zu FastFold-Verzeichnis
    cd external/fastfold
    
    # Lade nur Parameter (keine MSA-DBs)
    if [ -f "scripts/download_alphafold_params.sh" ]; then
        ./scripts/download_alphafold_params.sh ~/alphafold_dbs/
        log_success "AlphaFold-Parameter heruntergeladen"
    else
        log_error "download_alphafold_params.sh nicht gefunden"
        exit 1
    fi
    
    cd ../..
}

# Umgebungsvariablen setzen
setup_environment() {
    log_info "Setze Umgebungsvariablen..."
    
    # Erstelle .env Datei
    cat > .env << EOF
# AlphaFold-Parameter (ohne MSA-DBs)
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs

# ColabFold Remote-MMSeqs2
export COLABFOLD_REMOTE=1

# GPU-Einstellungen
export CUDA_VISIBLE_DEVICES=0

# GROMACS GPU
export GMX_GPU=0
EOF
    
    # Füge zu .bashrc hinzu
    if ! grep -q "ALPHAFOLD_DATA_DIR" ~/.bashrc; then
        echo "" >> ~/.bashrc
        echo "# amp-rbc-md Lightweight Setup" >> ~/.bashrc
        echo "export ALPHAFOLD_DATA_DIR=~/alphafold_dbs" >> ~/.bashrc
        echo "export COLABFOLD_REMOTE=1" >> ~/.bashrc
    fi
    
    log_success "Umgebungsvariablen gesetzt"
}

# Test-Installation
test_installation() {
    log_info "Teste Installation..."
    
    # Aktiviere Environment
    source ~/.bashrc
    conda activate amp-rbc-md
    
    # Teste ColabFold
    if python -c "import colabfold; print('ColabFold OK')" 2>/dev/null; then
        log_success "ColabFold funktioniert"
    else
        log_error "ColabFold Import fehlgeschlagen"
        exit 1
    fi
    
    # Teste GROMACS
    if command -v gmx &> /dev/null; then
        log_success "GROMACS verfügbar"
    else
        log_error "GROMACS nicht gefunden"
        exit 1
    fi
    
    # Teste amp-rbc-md
    if command -v amp-rbc-md &> /dev/null; then
        log_success "amp-rbc-md CLI verfügbar"
    else
        log_warning "amp-rbc-md CLI nicht gefunden - installiere mit 'pip install -e .'"
    fi
}

# Hauptfunktion
main() {
    log_info "Starte Lightweight Setup für 400GB VM..."
    
    check_disk_space
    setup_conda
    setup_colabfold
    download_alphafold_params
    setup_environment
    
    echo ""
    log_success "Lightweight Setup abgeschlossen!"
    echo ""
    log_info "Speicherplatz-Verbrauch:"
    echo "  - AlphaFold-Parameter: ~3.6 GB"
    echo "  - Repository + Code: ~1 GB"
    echo "  - Gesamt: < 5 GB"
    echo ""
    log_info "Nächste Schritte:"
    echo "1. Aktiviere Environment: conda activate amp-rbc-md"
    echo "2. Installiere Package: pip install -e ."
    echo "3. Teste: amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run"
    echo "4. Echte Simulation: amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1"
    echo ""
    log_warning "Hinweis: Diese Installation verwendet ColabFold-Batch mit Remote-MMSeqs2"
    log_warning "Keine lokalen MSA-Datenbanken nötig - MSA kommt vom Server"
    
    test_installation
}

# Führe Setup aus
main "$@" 