# amp-rbc-md  
[![CI](https://github.com/Hoimel1/amp-rbc-md/actions/workflows/ci.yml/badge.svg)](https://github.com/Hoimel1/amp-rbc-md/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)

Batch-fähige Martini-3-Simulations­pipeline für antimikrobielle Peptide an roten Blutkörperchen – inkl. Replika-Support, WHAM-Bootstrap-Analyse & automatischer Tox-Bewertung.

---

## 🚀 Quick Start

### Voraussetzungen
- **Linux** (Ubuntu 20.04+ empfohlen)
- **NVIDIA GPU** mit CUDA 12.1+
- **Conda** oder Miniconda
- **Python 3.10**

### Installation

```bash
# Repository klonen
git clone https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# Einfaches Setup
conda env create -f environment.yml
conda activate amp-rbc-md

# Editable install (entwickeln & testen)
pip install -e .

# Datenbanken für FastFold/OpenFold (einmalig)
# export OPENFOLD_DATA=/path/to/alphafold_dbs
```

### Verifikation

```bash
python verify-installation.py
```

### Erste Simulation

```bash
# Dry-run
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run

# Echte Simulation
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1
```

### PDB-mmCIF-Daten für FastFold herunterladen

FastFold / OpenFold benötigt die **experimentellen PDB-Strukturdaten** (mmCIF-Format). Einmalig herunterladen ⇒ ~95 GB :

```bash
# Zielverzeichnis festlegen
export OPENFOLD_DATA=$HOME/alphafold_dbs
mkdir -p "$OPENFOLD_DATA/pdb_mmcif/mmcif_files"

# Kompletten mmCIF-Mirror vom PDBe-rsync holen
rsync -rlpt -v -z --delete \
    rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ \
    "$OPENFOLD_DATA/pdb_mmcif/mmcif_files"

# (Empfohlen) SEQRES-Datei für Template-Suche
wget -P "$OPENFOLD_DATA/pdb_mmcif" \
    https://ftp.ebi.ac.uk/pub/databases/pdb/derived_data/pdb_seqres.txt
```

Danach die Variable dauerhaft setzen (z. B. in `~/.bashrc`):

```bash
echo 'export OPENFOLD_DATA="$HOME/alphafold_dbs"' >> ~/.bashrc
```

## 🧬 Theorie-Hintergrund

Die Insertions-Freie-Energie (ΔG_insert) wird aus Umbrella-Sampling-Fenstern mittels WHAM bestimmt. Ein Bootstrap (N=200) liefert Konfidenzintervalle, die mit heuristischen Schwellen aus `config/judge.yaml` verglichen werden, um toxische Kandidaten zu markieren.

**Strukturvorhersage**: Nutzt FastFold/OpenFold (PyTorch-basierte AlphaFold2-Reimplementation) – GPU-beschleunigt ohne JAX-Abhängigkeit.

## 🐳 Docker

```bash
# Build Image
docker build -t amp-rbc-md .

# Run mit GPU
docker run --gpus all -it amp-rbc-md
```

## ⚡ HPC/Slurm

```bash
#!/bin/bash
#SBATCH --gpus=1 --time=24:00:00
module load cuda/12.3 gromacs/2024-gpu

conda activate amp-rbc-md
export GMX_GPU=0
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 3 --profile chol_high --gpu 0
```

## 📊 MLflow-Tracking

```bash
export MLFLOW_TRACKING_URI=file:$(pwd)/mlruns
amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default
```

## 🔧 CLI-Features

| Flag | Beschreibung |
|------|--------------|
| `--seq` | Peptid-Sequenz |
| `-f / --fasta` | Batch-FASTA-Datei |
| `--n-replica` | Anzahl Replika |
| `--profile` | Lipid-Profil (default, chol_high, chol_low) |
| `--gpu` | GPU-ID |
| `-j / --workers` | Parallele Ausführung |
| `--resume` | Überspringe berechnete Replika |
| `--dry-run` | Nur Planung, keine Ausführung |

## 📁 Projekt-Struktur

```
amp-rbc-md/
├── src/                    # Python-Module
│   ├── fasta_to_pdb.py    # Strukturvorhersage
│   ├── martinize_wrap.py  # Martini-3-Wrapper
│   ├── build_membrane.py  # Membran-Aufbau
│   ├── gmx_runner.py      # GROMACS-Integration
│   ├── analyse.py         # WHAM-Bootstrap
│   └── judge.py           # Tox-Logik
├── workflow/              # Snakemake-Workflow
├── config/               # Konfigurationsdateien
├── examples/             # Beispiel-FASTA
├── tests/                # Unit-Tests
└── docs/                 # Dokumentation
```

## 🧪 Tests

```bash
pytest tests/
pytest --cov=src tests/
```

## 📚 Dokumentation

- `docs/TUTORIAL.md` - Schritt-für-Schritt-Anleitung
- `notebooks/quickstart.ipynb` - Interaktiver Workflow

## 🤝 Beitragen

1. Fork das Repository
2. Erstelle einen Feature-Branch
3. Committe deine Änderungen
4. Push zum Branch
5. Erstelle einen Pull Request

## 📄 Lizenz

Dieses Projekt ist unter der MIT-Lizenz lizenziert - siehe [LICENSE](LICENSE) für Details.

## 🙏 Danksagungen

- **ColabFold** für AlphaFold2-Implementierung
- **GROMACS** für MD-Simulationen
- **Martini** für Coarse-Grained-Kraftfeld

## 📜 Zitieren

Bitte zitieren Sie diese Software wie folgt (Platzhalter-DOI wird bei Release ersetzt):

```
Michel Hüller et al. amp-rbc-md: Martini-3 Workflow for Red Blood Cell Peptide Insertion. Version 1.0.0, 2025. DOI:10.5281/zenodo.1234567
```

---
