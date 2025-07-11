#!/bin/bash

# amp-rbc-md Ultra-Fast Setup für 400GB VM
# Minimale Installation mit nur essentiellen Paketen

set -e

echo "⚡ amp-rbc-md Ultra-Fast Setup (400GB VM)"
echo "========================================"

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
    
    if [ "$available_space" -lt 20 ]; then
        log_error "Nicht genug Speicherplatz! Benötigt: 20GB, Verfügbar: ${available_space}GB"
        exit 1
    fi
    
    log_success "Speicherplatz OK: ${available_space}GB verfügbar"
}

# Schnelles Conda Environment erstellen
setup_conda() {
    log_info "Erstelle minimales Conda Environment..."
    
    if command -v conda &> /dev/null; then
        # Erstelle Environment mit minimalen Paketen
        conda create -n amp-rbc-md python=3.10 -y
        log_success "Conda Environment erstellt"
        
        # Aktiviere und installiere Pakete einzeln
        source ~/.bashrc
        conda activate amp-rbc-md
        
        log_info "Installiere essentielle Pakete..."
        conda install -c conda-forge numpy pandas scipy click pyyaml -y
        conda install -c conda-forge gromacs=2024.5 -y
        conda install -c conda-forge pytorch=2.5.1 -y
        
        log_success "Essentielle Pakete installiert"
    else
        log_error "Conda nicht gefunden. Bitte installiere Miniconda/Anaconda."
        exit 1
    fi
}

# ColabFold installieren
setup_colabfold() {
    log_info "Installiere ColabFold..."
    
    # Aktiviere Environment
    source ~/.bashrc
    conda activate amp-rbc-md
    
    # Installiere ColabFold mit AlphaFold-Unterstützung
pip install "colabfold[alphafold]"
    
    log_success "ColabFold installiert"
}

# Projekt installieren
install_project() {
    log_info "Installiere amp-rbc-md Projekt..."
    
    # Aktiviere Environment
    source ~/.bashrc
    conda activate amp-rbc-md
    
    # Installiere Projekt
    pip install -e .
    
    log_success "Projekt installiert"
}

# Umgebungsvariablen setzen
setup_environment() {
    log_info "Setze Umgebungsvariablen..."
    
    # Erstelle .env Datei
    cat > .env << EOF
# ColabFold Remote-MMSeqs2
export COLABFOLD_REMOTE=1

# GPU-Einstellungen
export CUDA_VISIBLE_DEVICES=0

# GROMACS GPU
export GMX_GPU=0
EOF
    
    # Füge zu .bashrc hinzu
    if ! grep -q "COLABFOLD_REMOTE" ~/.bashrc; then
        echo "" >> ~/.bashrc
        echo "# amp-rbc-md Ultra-Fast Setup" >> ~/.bashrc
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
        log_warning "amp-rbc-md CLI nicht gefunden"
    fi
}

# Hauptfunktion
main() {
    log_info "Starte Ultra-Fast Setup für 400GB VM..."
    
    check_disk_space
    setup_conda
    setup_colabfold
    install_project
    setup_environment
    
    echo ""
    log_success "Ultra-Fast Setup abgeschlossen!"
    echo ""
    log_info "Speicherplatz-Verbrauch:"
    echo "  - Repository + Code: ~1 GB"
    echo "  - Conda Environment: ~2 GB"
    echo "  - Gesamt: < 5 GB"
    echo ""
    log_info "Nächste Schritte:"
    echo "1. Aktiviere Environment: conda activate amp-rbc-md"
    echo "2. Teste: amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run"
    echo "3. Echte Simulation: amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1"
    echo ""
    log_warning "Hinweis: Diese Installation verwendet ColabFold-Batch mit Remote-MMSeqs2"
    log_warning "Keine lokalen MSA-Datenbanken nötig - MSA kommt vom Server"
    
    test_installation
}

# Führe Setup aus
main "$@" 