#!/bin/bash

# Download minimale AlphaFold/FastFold Datenbanken für 400GB VM
# Verwendet reduced_dbs Modus (~300GB statt ~2TB)

set -e

echo "📥 Download minimale AlphaFold Datenbanken"
echo "=========================================="

# Farben für Output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

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

# Prüfe Voraussetzungen
check_prerequisites() {
    log_info "Prüfe Voraussetzungen..."
    
    # Prüfe aria2c
    if ! command -v aria2c &> /dev/null; then
        log_error "aria2c nicht gefunden. Installiere mit: sudo apt install aria2"
        exit 1
    fi
    
    # Prüfe Speicherplatz
    available_space=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "$available_space" -lt 350 ]; then
        log_error "Nicht genug Speicherplatz! Benötigt: 350GB, Verfügbar: ${available_space}GB"
        exit 1
    fi
    
    log_success "Voraussetzungen erfüllt"
}

# Erstelle Download-Verzeichnis
setup_directories() {
    log_info "Erstelle Download-Verzeichnisse..."
    
    mkdir -p ~/alphafold_dbs
    export ALPHAFOLD_DATA_DIR=~/alphafold_dbs
    
    log_success "Download-Verzeichnis: ~/alphafold_dbs"
}

# Lade minimale Datenbanken
download_databases() {
    log_info "Starte Download mit reduced_dbs Modus..."
    log_warning "Dies dauert mehrere Stunden (~300GB Download)"
    
    # Wechsle zu FastFold-Verzeichnis
    cd external/fastfold
    
    # Verwende reduced_dbs Modus für kleinere Datenbanken
    ./scripts/download_all_data.sh ~/alphafold_dbs reduced_dbs
    
    cd ../..
    
    log_success "Download abgeschlossen!"
}

# Setze Umgebungsvariablen
setup_environment() {
    log_info "Setze Umgebungsvariablen..."
    
    # Füge zu .bashrc hinzu
    if ! grep -q "ALPHAFOLD_DATA_DIR" ~/.bashrc; then
        echo "" >> ~/.bashrc
        echo "# AlphaFold/FastFold Datenbanken" >> ~/.bashrc
        echo "export ALPHAFOLD_DATA_DIR=~/alphafold_dbs" >> ~/.bashrc
    fi
    
    log_success "Umgebungsvariablen gesetzt"
}

# Teste Installation
test_databases() {
    log_info "Teste Datenbanken..."
    
    if [ -d "~/alphafold_dbs" ]; then
        db_size=$(du -sh ~/alphafold_dbs | cut -f1)
        log_success "Datenbanken installiert: ${db_size}"
        
        # Liste verfügbare Datenbanken
        echo ""
        log_info "Verfügbare Datenbanken:"
        ls -la ~/alphafold_dbs/
    else
        log_error "Datenbanken nicht gefunden"
        exit 1
    fi
}

# Hauptfunktion
main() {
    log_info "Starte Download minimaler Datenbanken..."
    
    check_prerequisites
    setup_directories
    download_databases
    setup_environment
    test_databases
    
    echo ""
    log_success "Minimale Datenbanken erfolgreich installiert!"
    echo ""
    log_info "Nächste Schritte:"
    echo "1. Aktiviere Environment: conda activate amp-rbc-md"
    echo "2. Teste Pipeline: amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run"
    echo "3. Starte echte Simulation: amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1"
    echo ""
    log_warning "Hinweis: Mit reduzierten Datenbanken ist die Qualität der Strukturvorhersage etwas geringer,"
    log_warning "aber für die meisten Peptide ausreichend."
}

# Führe Download aus
main "$@" 