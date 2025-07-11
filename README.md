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

### Installation

#### Option 1: Schritt-für-Schritt Installation (empfohlen) ⭐

```bash
# 1. Repository klonen (mit Submodulen)
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Folge der detaillierten Anleitung
docs/INSTALLATION_STEP_BY_STEP.md

# 3. Teste nach jedem Schritt
./test-installation.sh
```

**Vorteile:**
- ✅ **Kontrolle über jeden Schritt**
- ✅ **Einfaches Troubleshooting**
- ✅ **Verständnis der Abhängigkeiten**
- ✅ **Flexibel für verschiedene Systeme**

#### Option 2: Automatische Setup-Skripte

**Lightweight Installation (400GB VM):**
```bash
./setup-lightweight.sh
```

**Ultra-Fast Setup (falls Setup hängt):**
```bash
./setup-ultra-fast.sh
```

**Vollständige Installation (2TB+ VM):**
```bash
./setup.sh
```

### Erste Simulation

```bash
# Lightweight (ColabFold-Batch)
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run

# Vollständig (AlphaFold)
amp-rbc-md --seq GLSILGKLL --dry-run

# Echte Simulation
amp-rbc-md --seq GLSILGKLL --n-replica 1 --profile default -j 1
```

## 📋 Features

- **Strukturvorhersage**: Multiple Backends (ColabFold, FastFold, ESMFold, AlphaFold)
- **Coarse-Grained MD**: Martini-3 Kraftfeld für effiziente Simulationen
- **GPU-Beschleunigung**: Vollständige CUDA-Unterstützung
- **Batch-Verarbeitung**: Parallele Verarbeitung mehrerer Sequenzen
- **Automatisierte Pipeline**: Von FASTA zu MD-Trajektorien in einem Schritt
- **Flexible Installation**: Unterstützt 400GB (Lightweight) und 2TB+ (Vollständig)

## 🏗️ Architektur

```
amp-rbc-md/
├── src/amp_rbc_md/          # Hauptcode
│   ├── fasta_to_pdb.py      # Strukturvorhersage (Multi-Backend)
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

## 📊 Speicherplatz-Anforderungen

| Installation | Speicherplatz | Backend | Qualität | Empfehlung |
|--------------|---------------|---------|----------|------------|
| **Lightweight** | < 5 GB | ColabFold-Batch | Sehr gut | ✅ **400GB VM** |
| **FastFold No-MSA** | < 5 GB | FastFold | Sehr gut | ✅ **400GB VM** |
| **ESMFold** | < 3 GB | ESMFold | Gut | ✅ **Kurze Peptide** |
| **Vollständig** | 2TB+ | AlphaFold | Beste | ✅ **2TB+ VM** |

## 🔧 Installation-Optionen

### Lightweight Installation (400GB VM) - Empfohlen
- **Speicherplatz**: < 5 GB
- **Backend**: ColabFold-Batch mit Remote-MMSeqs2
- **Qualität**: Sehr gut für Peptide
- **Vorteil**: Sofort einsatzbereit, keine lokalen DBs

### Vollständige Installation (2TB+ VM)
- **Speicherplatz**: 2TB+
- **Backend**: AlphaFold mit allen MSA-Datenbanken
- **Qualität**: Beste Vorhersagequalität
- **Vorteil**: Offline-fähig, maximale Qualität

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
- [Optimization Guide](docs/OPTIMIZATION.md) - Speicherplatz-Optimierung & Backends
- [Tutorial](docs/TUTORIAL.md) - Schritt-für-Schritt Tutorial
- [API Documentation](docs/API.md) - Code-Dokumentation
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Häufige Probleme

## 🎯 Beispiele

### Lightweight (ColabFold-Batch)
```bash
# Einzel-Sequenz
amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 3 --profile default -j 3

# Batch-Verarbeitung
amp-rbc-md -f examples/batch.fasta --backend colabfold --n-replica 1 --profile default -j 4

# Dry-Run (Test)
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run
```

### Vollständig (AlphaFold)
```bash
# Einzel-Sequenz
amp-rbc-md --seq GLSILGKLL --n-replica 3 --profile default -j 3

# Batch-Verarbeitung
amp-rbc-md -f examples/batch.fasta --n-replica 1 --profile default -j 4
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
- **ColabFold**: Sergey Ovchinnikov für die optimierte Implementation
- **FastFold**: NVIDIA für die GPU-optimierte Implementation
- **Martini**: Marrink Lab für das Coarse-Grained Kraftfeld
- **GROMACS**: GROMACS Development Team für die MD-Software

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/Hoimel1/amp-rbc-md/issues)
- **Documentation**: [docs/](docs/)
- **Tutorial**: [docs/TUTORIAL.md](docs/TUTORIAL.md)

---

**Happy Simulating! 🧬⚡**
