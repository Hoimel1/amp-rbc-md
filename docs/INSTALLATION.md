# Installation Guide

## Übersicht

Diese Anleitung führt Sie durch die Installation von amp-rbc-md in verschiedenen Umgebungen:

1. **Vollständige Installation** (2TB+ Speicherplatz) - Beste Qualität
2. **Minimale Installation** (400GB Speicherplatz) - Reduzierte Qualität, aber funktionsfähig
3. **Docker Installation** - Isolierte Umgebung
4. **HPC Installation** - Für Cluster/Server

## Option 1: Vollständige Installation (Empfohlen)

### Voraussetzungen

- **Ubuntu 22.04** (empfohlen)
- **NVIDIA GPU** mit CUDA 12.1 Support
- **Miniconda/Anaconda**
- **Mindestens 2 TB freier Speicherplatz**

### Installation

```bash
# 1. Repository klonen (mit Submodulen)
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Setup-Skript ausführen
./setup.sh
```

Das Setup-Skript installiert automatisch:
- ✅ Conda-Environment mit allen Abhängigkeiten
- ✅ PyTorch mit CUDA 12.1 Support
- ✅ FastFold & OpenFold als Git-Submodule
- ✅ GROMACS 2024 mit GPU-Support
- ✅ Alle Python-Abhängigkeiten

### Datenbanken herunterladen

```bash
# Nach der Installation
cd external/fastfold
./scripts/download_all_data.sh $HOME/alphafold_dbs/
```

**Hinweis:** Der Download dauert mehrere Stunden und benötigt ~2 TB Speicherplatz.

## Option 2: Minimale Installation (400GB VM)

### Voraussetzungen

- **Ubuntu 22.04**
- **NVIDIA GPU** mit CUDA 12.1 Support
- **Miniconda/Anaconda**
- **Mindestens 400GB freier Speicherplatz**

### Installation

```bash
# 1. Repository klonen (mit Submodulen)
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Minimales Setup ausführen
./setup-minimal.sh
```

### Minimale Datenbanken herunterladen

```bash
# Nach der Installation
./download-minimal-dbs.sh
```

**Hinweis:** 
- Download dauert mehrere Stunden (~300GB)
- Verwende tmux für Hintergrund-Download: `tmux new-session -d './download-minimal-dbs.sh'`
- Qualität der Strukturvorhersage ist etwas geringer, aber für die meisten Peptide ausreichend

### Unterschiede zur vollständigen Installation

| Komponente | Vollständig | Minimal |
|------------|-------------|---------|
| BFD Datenbank | ~2TB | ~50GB (Small BFD) |
| Uniref90 | Vollständig | Reduziert |
| Uniref30 | Vollständig | Reduziert |
| MGnify | Vollständig | Reduziert |
| PDB70 | Vollständig | Reduziert |
| Gesamtspeicher | ~2TB | ~300GB |
| Vorhersagequalität | Beste | Gut |

## Option 3: Docker Installation

### Voraussetzungen

- **Docker** mit NVIDIA Container Runtime
- **NVIDIA GPU** mit CUDA 12.1 Support

### Installation

```bash
# 1. Repository klonen
git clone https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Docker Image bauen
docker build -t amp-rbc-md .

# 3. Container starten
docker run --gpus all -it amp-rbc-md
```

## Option 4: HPC Installation (Slurm)

### Voraussetzungen

- **HPC Cluster** mit Slurm
- **NVIDIA GPUs** verfügbar
- **Conda/Anaconda** verfügbar

### Installation

```bash
# 1. Repository klonen
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Setup ausführen
./setup.sh

# 3. Snakemake-Profil konfigurieren
# Bearbeite workflow/profile/slurm/config.yaml
```

### Ausführung auf HPC

```bash
# Batch-Job starten
snakemake --profile workflow/profile/slurm -j 64
```

## Nach der Installation

### Environment aktivieren

```bash
conda activate amp-rbc-md
```

### Package installieren

```bash
pip install -e .
```

### Erste Simulation

```bash
# Dry-Run testen
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run

# Echte Simulation
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1
```

## Troubleshooting

### GPU-Probleme

```bash
# Prüfe CUDA-Installation
nvidia-smi

# Prüfe PyTorch CUDA-Support
python -c "import torch; print(torch.cuda.is_available())"

# Setze GPU-ID
export CUDA_VISIBLE_DEVICES=0
```

### Speicherplatz-Probleme

```bash
# Prüfe verfügbaren Speicherplatz
df -h

# Lösche temporäre Dateien
rm -rf /tmp/*
```

### Conda-Probleme

```bash
# Environment neu erstellen
conda env remove -n amp-rbc-md
conda env create -f environment.yml

# Abhängigkeiten aktualisieren
conda update --all
```

### Datenbank-Probleme

```bash
# Prüfe Datenbank-Pfad
echo $ALPHAFOLD_DATA_DIR

# Setze Datenbank-Pfad
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs
```

## Support

Bei Problemen:

1. Prüfe die [Troubleshooting-Sektion](#troubleshooting)
2. Schaue in die [Logs](docs/TROUBLESHOOTING.md)
3. Erstelle ein [Issue](https://github.com/Hoimel1/amp-rbc-md/issues)

## Nächste Schritte

Nach erfolgreicher Installation:

1. Lese das [Tutorial](docs/TUTORIAL.md)
2. Teste mit [Beispiel-Sequenzen](examples/)
3. Konfiguriere [Batch-Verarbeitung](docs/BATCH.md)
4. Optimiere für [HPC](docs/HPC.md) 