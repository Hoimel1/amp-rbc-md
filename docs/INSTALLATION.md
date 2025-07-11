# Installation Guide

## Übersicht

AMP-RBC-MD ist eine vollständige Pipeline für die Simulation antimikrobieller Peptide. Diese Anleitung führt Sie durch die komplette Installation.

## Systemanforderungen

### Hardware
- **CPU**: 8+ Cores (empfohlen)
- **RAM**: 32+ GB (empfohlen)
- **GPU**: NVIDIA GPU mit CUDA 12.1 Support
- **Speicher**: 2+ TB freier Speicherplatz

### Software
- **OS**: Ubuntu 22.04 (empfohlen)
- **Conda**: Miniconda oder Anaconda
- **Git**: Für Repository-Klon

## Schritt-für-Schritt Installation

### 1. System vorbereiten

```bash
# System aktualisieren
sudo apt update && sudo apt upgrade -y

# NVIDIA-Treiber installieren (falls nicht vorhanden)
sudo apt install nvidia-driver-535

# CUDA Toolkit installieren
wget https://developer.download.nvidia.com/compute/cuda/12.1.0/local_installers/cuda_12.1.0_530.30.02_linux.run
sudo sh cuda_12.1.0_530.30.02_linux.run
```

### 2. Miniconda installieren

```bash
# Miniconda herunterladen
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Installieren
bash Miniconda3-latest-Linux-x86_64.sh

# Shell neu laden
source ~/.bashrc
```

### 3. Repository klonen

```bash
# Repository mit Submodulen klonen
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md
```

### 4. Automatisches Setup

```bash
# Setup-Skript ausführen
./setup.sh
```

Das Setup-Skript:
- ✅ Erstellt das Conda-Environment
- ✅ Installiert alle Python-Abhängigkeiten
- ✅ Installiert FastFold & OpenFold
- ✅ Konfiguriert Umgebungsvariablen
- ✅ Testet die Installation

### 5. Datenbanken herunterladen

```bash
# Nach der Installation
cd external/fastfold
./scripts/download_all_data.sh $HOME/alphafold_dbs/
```

**Wichtig**: Der Download dauert mehrere Stunden und benötigt ~2 TB Speicherplatz.

### 6. Installation testen

```bash
# Zurück ins Projekt-Verzeichnis
cd ~/amp-rbc-md

# Dry-Run testen
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run
```

## Manuelle Installation (Alternative)

Falls das automatische Setup nicht funktioniert:

```bash
# 1. Environment erstellen
conda env create -f environment.yml
conda activate amp-rbc-md

# 2. Git-Submodule initialisieren
git submodule update --init --recursive

# 3. FastFold installieren
cd external/fastfold
pip install -e . --no-build-isolation
cd ../..

# 4. OpenFold installieren
cd external/openfold
pip install -e . --no-build-isolation
cd ../..

# 5. Projekt installieren
pip install -e .
```

## Verifikation

### GPU-Support prüfen

```bash
# CUDA verfügbar?
nvidia-smi

# PyTorch CUDA-Support
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

### Installation testen

```bash
# Verifikationsskript
python verify-installation.py --engine fastfold
```

### Erste Simulation

```bash
# Dry-Run
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run

# Echte Simulation
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1
```

## Troubleshooting

### Häufige Probleme

**GPU-Probleme:**
```bash
# CUDA-Version prüfen
nvcc --version

# PyTorch neu installieren
pip uninstall torch torchvision torchaudio
pip install torch==2.3.0+cu121 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

**Speicherplatz-Probleme:**
```bash
# Speicherplatz prüfen
df -h

# Nur essentielle Datenbanken
cd external/fastfold
./scripts/download_pdb_mmcif.sh $HOME/alphafold_dbs/
```

**Conda-Probleme:**
```bash
# Environment neu erstellen
conda env remove -n amp-rbc-md
conda env create -f environment.yml
```

### Support

Bei Problemen:
1. Prüfe die [Troubleshooting-Seite](../TROUBLESHOOTING.md)
2. Öffne ein [GitHub Issue](https://github.com/Hoimel1/amp-rbc-md/issues)
3. Prüfe die Logs in `results/`

## Nächste Schritte

Nach erfolgreicher Installation:

1. **Tutorial**: Siehe [Tutorial](TUTORIAL.md)
2. **Beispiele**: Teste die Beispiel-Sequenzen in `examples/`
3. **Konfiguration**: Passe die Profile in `config/` an
4. **Dokumentation**: Lese die [API-Dokumentation](API.md) 