# amp-rbc-md

Batch-fähige Martini-3-Simulation roter Blutkörperchen-Peptide (≤10) mit Replika-Support, WHAM-Bootstrap-Analyse und konfigurierbarer Tox-Logik.

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

# Automatisches Setup (empfohlen)
bash setup.sh

# Oder manuell:
conda env create -f environment.yml
conda activate amp-rbc-md
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install nvidia-cudnn-cu12==8.9.2.26
pip install colabfold==1.5.5 --no-deps
pip install absl-py appdirs biopython matplotlib numpy pandas py3Dmol requests tqdm dm-haiku
pip install -e .
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

## 🧬 Theorie-Hintergrund

Die Insertions-Freie-Energie (ΔG_insert) wird aus Umbrella-Sampling-Fenstern mittels WHAM bestimmt. Ein Bootstrap (N=200) liefert Konfidenzintervalle, die mit heuristischen Schwellen aus `config/judge.yaml` verglichen werden, um toxische Kandidaten zu markieren.

**Strukturvorhersage**: Verwendet ColabFold (AlphaFold2-Implementierung) für State-of-the-Art Protein-Strukturvorhersage.

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
