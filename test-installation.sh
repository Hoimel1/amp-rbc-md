#!/bin/bash

# amp-rbc-md Installation Test Script
# Führt nach jedem Schritt aus, um Fortschritt zu prüfen

set -e

echo "🧪 amp-rbc-md Installation Test"
echo "================================"

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

# Test 1: System
test_system() {
    log_info "Teste System..."
    
    # Ubuntu Version
    if lsb_release -d | grep -q "Ubuntu 22.04"; then
        log_success "Ubuntu 22.04 erkannt"
    else
        log_warning "Ubuntu Version: $(lsb_release -d)"
    fi
    
    # NVIDIA GPU
    if command -v nvidia-smi &> /dev/null; then
        if nvidia-smi &> /dev/null; then
            log_success "NVIDIA GPU verfügbar"
            nvidia-smi --query-gpu=name,memory.total --format=csv,noheader,nounits
        else
            log_error "NVIDIA GPU nicht erkannt"
        fi
    else
        log_error "nvidia-smi nicht gefunden"
    fi
    
    # Speicherplatz
    available_space=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
    log_info "Verfügbarer Speicherplatz: ${available_space}GB"
    
    if [ "$available_space" -lt 50 ]; then
        log_warning "Wenig Speicherplatz verfügbar"
    else
        log_success "Ausreichend Speicherplatz"
    fi
}

# Test 2: Conda
test_conda() {
    log_info "Teste Conda..."
    
    if command -v conda &> /dev/null; then
        log_success "Conda verfügbar: $(conda --version)"
        
        # Prüfe Environment
        if conda env list | grep -q "amp-rbc-md"; then
            log_success "amp-rbc-md Environment existiert"
        else
            log_warning "amp-rbc-md Environment nicht gefunden"
        fi
    else
        log_error "Conda nicht gefunden"
    fi
}

# Test 3: Python Packages
test_python_packages() {
    log_info "Teste Python Packages..."
    
    # Teste essentielle Packages (ohne conda activate im Skript)
    packages=("numpy" "pandas" "scipy" "matplotlib" "click" "pyyaml" "rich" "tqdm")
    
    for package in "${packages[@]}"; do
        if python -c "import $package; print('$package OK')" 2>/dev/null; then
            log_success "$package verfügbar"
        else
            log_warning "$package nicht verfügbar"
        fi
    done
}

# Test 4: GROMACS
test_gromacs() {
    log_info "Teste GROMACS..."
    
    if command -v gmx &> /dev/null; then
        log_success "GROMACS verfügbar"
        gmx --version | head -1
    else
        log_warning "GROMACS nicht gefunden"
    fi
}

# Test 5: PyTorch + CUDA
test_pytorch() {
    log_info "Teste PyTorch + CUDA..."
    
    if python -c "import torch; print(f'PyTorch: {torch.__version__}')" 2>/dev/null; then
        log_success "PyTorch verfügbar"
        
        if python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')" 2>/dev/null; then
            if python -c "import torch; print(torch.cuda.is_available())" | grep -q "True"; then
                log_success "CUDA verfügbar"
                python -c "import torch; print(f'CUDA devices: {torch.cuda.device_count()}')"
            else
                log_warning "CUDA nicht verfügbar"
            fi
        else
            log_warning "CUDA-Test fehlgeschlagen"
        fi
    else
        log_warning "PyTorch nicht verfügbar"
    fi
}

# Test 6: Submodules
test_submodules() {
    log_info "Teste Git Submodules..."
    
    if [ -d "external/fastfold" ]; then
        log_success "FastFold Submodule vorhanden"
    else
        log_warning "FastFold Submodule nicht gefunden"
    fi
    
    if [ -d "external/openfold" ]; then
        log_success "OpenFold Submodule vorhanden"
    else
        log_warning "OpenFold Submodule nicht gefunden"
    fi
}

# Test 7: ColabFold
test_colabfold() {
    log_info "Teste ColabFold..."
    
    if python -c "import colabfold; print('ColabFold OK')" 2>/dev/null; then
        log_success "ColabFold verfügbar"
    else
        log_warning "ColabFold nicht verfügbar"
    fi
}

# Test 8: amp-rbc-md
test_amp_rbc_md() {
    log_info "Teste amp-rbc-md..."
    
    if command -v amp-rbc-md &> /dev/null; then
        log_success "amp-rbc-md CLI verfügbar"
        amp-rbc-md --help | head -5
    else
        log_warning "amp-rbc-md CLI nicht gefunden"
    fi
}

# Test 9: Environment Variables
test_env_vars() {
    log_info "Teste Umgebungsvariablen..."
    
    if [ "$COLABFOLD_REMOTE" = "1" ]; then
        log_success "COLABFOLD_REMOTE gesetzt"
    else
        log_warning "COLABFOLD_REMOTE nicht gesetzt"
    fi
    
    if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
        log_success "CUDA_VISIBLE_DEVICES gesetzt: $CUDA_VISIBLE_DEVICES"
    else
        log_warning "CUDA_VISIBLE_DEVICES nicht gesetzt"
    fi
}

# Hauptfunktion
main() {
    echo ""
    test_system
    echo ""
    test_conda
    echo ""
    test_python_packages
    echo ""
    test_gromacs
    echo ""
    test_pytorch
    echo ""
    test_submodules
    echo ""
    test_colabfold
    echo ""
    test_amp_rbc_md
    echo ""
    test_env_vars
    echo ""
    
    log_info "Installation Test abgeschlossen!"
    echo ""
    log_info "Nächste Schritte:"
    echo "1. Führe fehlende Schritte aus der Anleitung aus"
    echo "2. Führe diesen Test erneut aus"
    echo "3. Starte erste Simulation: amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run"
}

# Führe Tests aus
main "$@" 