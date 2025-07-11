# AMP-RBC-MD

**Antimicrobial Peptide RBC Membrane Disruption - Molecular Dynamics Pipeline**

Eine vollständige Pipeline zur Simulation von antimikrobiellen Peptiden und deren Interaktion mit roten Blutkörperchen (RBC) unter Verwendung von Martini-3 Coarse-Grained Molecular Dynamics.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/)

## 🚀 Schnellstart

### Voraussetzungen

- **Ubuntu 22.04** (empfohlen)
- **NVIDIA GPU** mit CUDA 11.5+ Support
- **Miniconda/Anaconda**

### Installation

```bash
# 1. Repository klonen (mit Submodulen)
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Conda-Umgebung erstellen (schlank, nur FastFold)
conda env create -f environment.yml

# 3. Umgebung aktivieren
conda activate amp-rbc-md

# 4. FastFold-Parameter herunterladen
mkdir -p $HOME/.fastfold
cd $HOME/.fastfold
wget https://github.com/hpcaitech/FastFold/releases/download/v0.2.0/fastfold_params.tar.gz
tar -xzf fastfold_params.tar.gz
cd -

# 5. Umgebungsvariablen setzen
echo 'export FASTFOLD_SKIP_TEMPLATES=1' >> ~/.zshrc
echo 'export FASTFOLD_NO_MSA=1' >> ~/.zshrc
echo 'export FASTFOLD_PARAMS_PATH=$HOME/.fastfold/fastfold_params' >> ~/.zshrc
source ~/.zshrc

# 6. Testen
amp-rbc-md --seq GLSILGKLL --dry-run
```

### Erste Simulation

```bash
# Einzel-Simulation
amp-rbc-md --seq GLSILGKLL --n-replica 1 --profile default -j 1

# Batch-Simulation
amp-rbc-md -f examples/batch.fasta --n-replica 1 --profile default -j 3
```

## 📋 Features

- **Strukturvorhersage**: FastFold ohne MSA (schnell und effizient)
- **Coarse-Grained MD**: Martini-3 Kraftfeld für effiziente Simulationen
- **GPU-Beschleunigung**: Vollständige CUDA-Unterstützung
- **Batch-Verarbeitung**: Parallele Verarbeitung mehrerer Sequenzen
- **Automatisierte Pipeline**: Von FASTA zu MD-Trajektorien in einem Schritt
- **Schlanke Installation**: Nur ~5GB Speicherplatz benötigt

## 🏗️ Architektur

```
amp-rbc-md/
├── src/amp_rbc_md/          # Hauptcode
│   ├── fasta_to_pdb.py      # Strukturvorhersage (FastFold)
│   ├── martinize_wrap.py    # CG-Konvertierung (Martini-3)
│   ├── build_membrane.py    # Membran-Aufbau
│   ├── gen_mdp.py          # GROMACS-Parameter
│   ├── gmx_runner.py       # GROMACS-Ausführung
│   └── analyse.py          # Trajektorien-Analyse
├── external/                # Git-Submodule
│   └── fastfold/           # FastFold (AlphaFold2-Implementation)
├── config/                 # Konfigurationsdateien
├── mdp_templates/          # GROMACS-Templates
└── workflow/               # Snakemake-Workflow (HPC)
```

## 📊 Speicherplatz-Anforderungen

| Komponente | Speicherplatz | Beschreibung |
|------------|---------------|--------------|
| **FastFold-Parameter** | ~3.6 GB | AlphaFold2-Parameter |
| **Conda-Umgebung** | ~2 GB | Python-Pakete |
| **Repository + Code** | ~1 GB | Quellcode |
| **Gesamt** | **~6.6 GB** | Vollständige Installation |

## 🔧 Installation-Optionen

### Standard Installation (empfohlen)
- **Speicherplatz**: ~6.6 GB
- **Backend**: FastFold ohne MSA
- **Qualität**: Sehr gut für Peptide bis ~50 Aminosäuren
- **Vorteil**: Schnell, schlank, offline-fähig

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
# Dry-Run (Test)
amp-rbc-md --seq GLSILGKLL --dry-run

# Echte Simulation
amp-rbc-md --seq GLSILGKLL --n-replica 3 --profile default -j 3
```

### Batch-Verarbeitung
```bash
# Batch-Simulation
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

- **FastFold**: NVIDIA für die GPU-optimierte AlphaFold2-Implementation
- **Martini**: Marrink Lab für das Coarse-Grained Kraftfeld
- **GROMACS**: GROMACS Development Team für die MD-Software

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/Hoimel1/amp-rbc-md/issues)
- **Documentation**: [docs/](docs/)
- **Tutorial**: [docs/TUTORIAL.md](docs/TUTORIAL.md)

---

**Happy Simulating! 🧬⚡**
