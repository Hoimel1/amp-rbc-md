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
│   ├── fasta_to_pdb.py      # Strukturvorhersage (FastFold/OpenFold)
│   ├── martinize_wrap.py    # CG-Konvertierung (Martini-3)
│   ├── build_membrane.py    # Membran-Aufbau
│   ├── gen_mdp.py          # GROMACS-Parameter
│   ├── gmx_runner.py       # GROMACS-Ausführung
│   └── analyse.py          # Trajektorien-Analyse
├── external/                # Git-Submodule
│   ├── fastfold/           # FastFold (AlphaFold2-Implementation)
│   └── openfold/           # OpenFold (Alternative)
├── config/                 # Konfigurationsdateien
├── mdp_templates/          # GROMACS-Templates
└── workflow/               # Snakemake-Workflow (HPC)
```

## 🔧 Installation-Optionen

### Vollständige Installation (Empfohlen)
- **Speicherplatz**: 2TB+
- **Qualität**: Beste Strukturvorhersage
- **Verwendung**: Produktion, Forschung
- **Datenbanken**: Alle AlphaFold-Datenbanken

### Docker Installation
```bash
docker run --gpus all mydockerhub/amp-rbc-md:latest
```

### HPC Installation (Slurm)
```bash
snakemake --profile workflow/profile/slurm -j 64
```

## 📖 Dokumentation

- [Installation Guide](docs/INSTALLATION.md) - Detaillierte Installationsanleitung
- [Tutorial](docs/TUTORIAL.md) - Schritt-für-Schritt Tutorial
- [API Documentation](docs/API.md) - Code-Dokumentation
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Häufige Probleme

## 🎯 Beispiele

### Einzel-Sequenz
```bash
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 3 --profile default -j 3
```

### Batch-Verarbeitung
```bash
amp-rbc-md -f examples/batch.fasta --n-replica 1 --profile default -j 4
```

### Dry-Run (Test)
```bash
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run
```

## 🤝 Beitragen

1. Fork das Repository
2. Erstelle einen Feature-Branch (`git checkout -b feature/AmazingFeature`)
3. Committe deine Änderungen (`git commit -m 'Add some AmazingFeature'`)
4. Push zum Branch (`git push origin feature/AmazingFeature`)
5. Öffne einen Pull Request

## 📄 Lizenz

Dieses Projekt ist unter der MIT-Lizenz lizenziert - siehe [LICENSE](LICENSE) Datei für Details.

## 🙏 Danksagungen

- **AlphaFold2**: DeepMind für die revolutionäre Strukturvorhersage
- **FastFold**: NVIDIA für die GPU-optimierte Implementation
- **Martini**: Marrink Lab für das Coarse-Grained Kraftfeld
- **GROMACS**: GROMACS Development Team für die MD-Software

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/Hoimel1/amp-rbc-md/issues)
- **Documentation**: [docs/](docs/)
- **Tutorial**: [docs/TUTORIAL.md](docs/TUTORIAL.md)

---

**Happy Simulating! 🧬⚡**
