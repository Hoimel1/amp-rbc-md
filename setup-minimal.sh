#!/bin/bash

# amp-rbc-md Minimal Setup für 400GB VM
# Verwendet reduzierte Datenbanken für AlphaFold/FastFold

set -e

echo "🚀 amp-rbc-md Minimal Setup (400GB VM)"
echo "======================================"

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
    
    if [ "$available_space" -lt 350 ]; then
        log_error "Nicht genug Speicherplatz! Benötigt: 350GB, Verfügbar: ${available_space}GB"
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

# FastFold Submodule initialisieren
setup_fastfold() {
    log_info "Initialisiere FastFold Submodule..."
    
    if [ ! -d "external/fastfold" ]; then
        log_error "FastFold Submodule nicht gefunden. Führe 'git clone --recursive' aus."
        exit 1
    fi
    
    log_success "FastFold Submodule bereit"
}

# Minimale Datenbanken herunterladen
download_minimal_databases() {
    log_info "Lade minimale Datenbanken herunter (~300GB)..."
    
    # Erstelle Download-Verzeichnis
    mkdir -p ~/alphafold_dbs
    
    # Wechsle zu FastFold-Verzeichnis
    cd external/fastfold
    
    # Verwende reduced_dbs Modus
    log_info "Starte Download mit reduced_dbs Modus..."
    ./scripts/download_all_data.sh ~/alphafold_dbs reduced_dbs
    
    cd ../..
    
    log_success "Minimale Datenbanken heruntergeladen"
}

# Umgebungsvariablen setzen
setup_environment() {
    log_info "Setze Umgebungsvariablen..."
    
    # Erstelle .env Datei
    cat > .env << EOF
# AlphaFold/FastFold Datenbanken
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs

# GPU-Einstellungen
export CUDA_VISIBLE_DEVICES=0

# GROMACS GPU
export GMX_GPU=0
EOF
    
    log_success "Umgebungsvariablen gesetzt"
}

# Test-Installation
test_installation() {
    log_info "Teste Installation..."
    
    # Aktiviere Environment
    source ~/.bashrc
    conda activate amp-rbc-md
    
    # Teste FastFold
    if python -c "import fastfold; print('FastFold OK')" 2>/dev/null; then
        log_success "FastFold funktioniert"
    else
        log_warning "FastFold Import fehlgeschlagen - wird bei Bedarf installiert"
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
    log_info "Starte Minimal Setup für 400GB VM..."
    
    check_disk_space
    setup_conda
    setup_fastfold
    setup_environment
    
    echo ""
    log_info "Nächste Schritte:"
    echo "1. Aktiviere Environment: conda activate amp-rbc-md"
    echo "2. Installiere Package: pip install -e ."
    echo "3. Lade Datenbanken: ./download-minimal-dbs.sh"
    echo "4. Teste: amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run"
    echo ""
    log_warning "Hinweis: Download der Datenbanken dauert mehrere Stunden!"
    log_warning "Verwende tmux für Hintergrund-Download: tmux new-session -d './download-minimal-dbs.sh'"
    
    test_installation
}

# Führe Setup aus
main "$@" 