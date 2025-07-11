# AMP-RBC-MD

**Antimicrobial Peptide RBC Membrane Disruption - Molecular Dynamics Pipeline**

Eine vollständige Pipeline zur Simulation von antimikrobiellen Peptiden und deren Interaktion mit roten Blutkörperchen (RBC) unter Verwendung von Martini-3 Coarse-Grained Molecular Dynamics.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/)

## 🚀 Schnellstart

### Voraussetzungen

- **Ubuntu 22.04** (empfohlen)
- **NVIDIA GPU** mit CUDA 12.1 Support
- **Miniconda/Anaconda**
- **Mindestens 2 TB freier Speicherplatz** (für alle Datenbanken)

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

### Erste Simulation

```bash
# Dry-Run testen
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run

# Echte Simulation
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1
```

## 📋 Features

- **Strukturvorhersage**: FastFold/OpenFold für AlphaFold2-basierte Strukturvorhersage
- **Coarse-Grained MD**: Martini-3 Kraftfeld für effiziente Simulationen
- **GPU-Beschleunigung**: Vollständige CUDA-Unterstützung
- **Batch-Verarbeitung**: Parallele Verarbeitung mehrerer Sequenzen
- **Automatisierte Pipeline**: Von FASTA zu MD-Trajektorien in einem Schritt

## 🏗️ Architektur

```
amp-rbc-md/
├── src/amp_rbc_md/          # Hauptcode
│   ├── fastfold_wrap.py     # FastFold-Integration
│   ├── fasta_to_pdb.py      # Strukturvorhersage
│   └── ...
├── external/                # Git-Submodule
│   ├── fastfold/           # FastFold (optimiertes AlphaFold2)
│   ├── openfold/           # OpenFold (Referenz-Implementation)
│   └── martinize2/         # Martini-3 Topologie-Generator
├── workflow/               # Snakemake-Pipeline
├── mdp_templates/          # GROMACS-Parameter
└── environment.yml         # Conda-Environment
```

## 🔧 Konfiguration

### Profile

Verschiedene Simulationsprofile sind verfügbar:

```bash
# Standard-Profil
amp-rbc-md --seq <sequence> --profile default

# GPU-optimiertes Profil
amp-rbc-md --seq <sequence> --profile gpu

# Lange Trajektorien
amp-rbc-md --seq <sequence> --profile long
```

### Batch-Verarbeitung

```bash
# Mehrere Sequenzen aus FASTA-Datei
amp-rbc-md -f sequences.fasta --n-replica 3 -j 4

# Einzelne Sequenz mit mehreren Replikaten
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 5 -j 2
```

## 📊 Ausgabe

Die Pipeline erzeugt:

- **Strukturen**: PDB-Dateien (AlphaFold2-Vorhersagen)
- **Topologien**: Martini-3 CG-Topologien
- **Trajektorien**: GROMACS-Trajektorien (.xtc, .trr)
- **Analysen**: Automatische Analyse der Membraninteraktionen
- **Reports**: Zusammenfassende Berichte und Visualisierungen

## 🐛 Troubleshooting

### Häufige Probleme

**GPU-Probleme:**
```bash
# CUDA-Version prüfen
nvidia-smi
nvcc --version

# PyTorch CUDA-Support testen
python -c "import torch; print(torch.cuda.is_available())"
```

**Speicherplatz-Probleme:**
```bash
# Speicherplatz prüfen
df -h

# Nur essentielle Datenbanken herunterladen
cd external/fastfold
./scripts/download_pdb_mmcif.sh $HOME/alphafold_dbs/
```

**Installations-Probleme:**
```bash
# Environment neu erstellen
conda env remove -n amp-rbc-md
./setup.sh
```

## 📚 Dokumentation

- [Tutorial](docs/TUTORIAL.md) - Schritt-für-Schritt Anleitung
- [API-Dokumentation](docs/API.md) - Code-Referenz
- [Troubleshooting](TROUBLESHOOTING.md) - Häufige Probleme und Lösungen

## 🤝 Beitragen

1. Fork das Repository
2. Erstelle einen Feature-Branch (`git checkout -b feature/amazing-feature`)
3. Committe deine Änderungen (`git commit -m 'Add amazing feature'`)
4. Push zum Branch (`git push origin feature/amazing-feature`)
5. Öffne einen Pull Request

## 📄 Lizenz

Dieses Projekt ist unter der MIT-Lizenz lizenziert - siehe [LICENSE](LICENSE) Datei für Details.

## 🙏 Danksagungen

- **FastFold**: Optimierte AlphaFold2-Implementation
- **OpenFold**: Referenz-Implementation von AlphaFold2
- **Martini**: Coarse-Grained Kraftfeld
- **GROMACS**: Molecular Dynamics Engine

## 📖 Zitierung

Wenn du AMP-RBC-MD in deiner Forschung verwendest, zitiere bitte:

```bibtex
@software{amp_rbc_md,
  title={AMP-RBC-MD: Antimicrobial Peptide RBC Membrane Disruption Pipeline},
  author={Hoimel, Michel},
  year={2024},
  url={https://github.com/Hoimel1/amp-rbc-md}
}
```
