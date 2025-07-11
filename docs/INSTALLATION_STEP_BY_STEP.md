# amp-rbc-md: Schritt-für-Schritt Installation

Diese Anleitung führt Sie durch die manuelle Installation von amp-rbc-md auf Ubuntu 22.04 mit NVIDIA GPU.

## 🎯 Übersicht

Wir installieren:
1. **Basis-Environment** (Conda + Python)
2. **Core-Dependencies** (GROMACS, PyTorch, etc.)
3. **Submodule** (FastFold, OpenFold)
4. **ColabFold** (für Strukturvorhersage)
5. **amp-rbc-md** (unser Projekt)

## 📋 Voraussetzungen

- Ubuntu 22.04 LTS
- NVIDIA GPU (L4 oder besser)
- Mindestens 400GB Speicherplatz
- Internetverbindung

## 🚀 Schritt 1: System vorbereiten

```bash
# System aktualisieren
sudo apt update && sudo apt upgrade -y

# NVIDIA Treiber prüfen
nvidia-smi

# Miniconda installieren (falls nicht vorhanden)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## 🐍 Schritt 2: Conda Environment erstellen

```bash
# Repository klonen
git clone https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# Basis-Environment erstellen
conda create -n amp-rbc-md python=3.10 -y
conda activate amp-rbc-md
```

## 🔧 Schritt 3: Core-Dependencies installieren

```bash
# Wissenschaftliche Basis-Pakete
conda install -c conda-forge numpy pandas scipy matplotlib -y

# CLI-Tools
conda install -c conda-forge click pyyaml rich tqdm -y

# GROMACS (MD-Simulation) - aus conda-forge, nicht bioconda
conda install -c conda-forge gromacs=2024.5 -y

# PyTorch (aus conda-forge, nicht nvidia)
conda install -c conda-forge pytorch=2.5.1 -y

# System-Tools
conda install -c conda-forge wget curl git -y
```

## 📦 Schritt 4: Submodule initialisieren

```bash
# Git Submodule initialisieren
git submodule update --init --recursive

# Prüfe Submodule-Status
git submodule status
```

## 🧬 Schritt 5: ColabFold installieren

```bash
# ColabFold für Remote-Strukturvorhersage
pip install colabfold

# Teste Installation
python -c "import colabfold; print('ColabFold OK')"
```

## 🎯 Schritt 6: amp-rbc-md installieren

```bash
# Projekt installieren
pip install -e .

# Teste CLI
amp-rbc-md --help
```

## ⚙️ Schritt 7: Umgebungsvariablen setzen

```bash
# Erstelle .env Datei
cat > .env << EOF
# ColabFold Remote-MMSeqs2
export COLABFOLD_REMOTE=1

# GPU-Einstellungen
export CUDA_VISIBLE_DEVICES=0

# GROMACS GPU
export GMX_GPU=0
EOF

# Füge zu .bashrc hinzu
echo "" >> ~/.bashrc
echo "# amp-rbc-md Environment" >> ~/.bashrc
echo "export COLABFOLD_REMOTE=1" >> ~/.bashrc
echo "export CUDA_VISIBLE_DEVICES=0" >> ~/.bashrc
```

## 🧪 Schritt 8: Installation testen

```bash
# Aktiviere Environment
conda activate amp-rbc-md

# Teste alle Komponenten
echo "=== Testing Components ==="

# GROMACS
gmx --version

# PyTorch + CUDA
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}')"

# ColabFold
python -c "import colabfold; print('ColabFold: OK')"

# amp-rbc-md
amp-rbc-md --help

echo "=== All tests completed ==="
```

## 🎉 Schritt 9: Erste Simulation

```bash
# Dry-Run Test
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run

# Echte Simulation (klein)
amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1
```

## 🔍 Troubleshooting

### Problem: Conda hängt beim "Solving environment"
**Lösung:** Verwende `setup-ultra-fast.sh` oder installiere Pakete einzeln:
```bash
conda install -c conda-forge numpy -y
conda install -c conda-forge pandas -y
# ... weitere Pakete einzeln
```

### Problem: GROMACS nicht in bioconda verfügbar
**Lösung:** GROMACS ist in conda-forge verfügbar:
```bash
conda install -c conda-forge gromacs=2024.5 -y
```

### Problem: CUDA nicht erkannt
**Lösung:** Prüfe NVIDIA-Treiber und PyTorch-Installation:
```bash
nvidia-smi
python -c "import torch; print(torch.cuda.is_available())"
```

### Problem: PyTorch nicht in nvidia Channel verfügbar
**Lösung:** PyTorch ist in conda-forge verfügbar:
```bash
conda install -c conda-forge pytorch=2.5.1 -y
```

### Problem: GROMACS nicht gefunden
**Lösung:** Reinstalliere GROMACS:
```bash
conda install -c conda-forge gromacs=2024.5 -y
```

### Problem: ColabFold Import-Fehler
**Lösung:** Reinstalliere ColabFold:
```bash
pip uninstall colabfold -y
pip install colabfold
```

### Problem: batchfold nicht verfügbar
**Lösung:** batchfold ist nicht nötig, nur colabfold installieren:
```bash
pip install colabfold
```

## 📊 Speicherplatz-Übersicht

- **Repository + Code:** ~1 GB
- **Conda Environment:** ~3 GB
- **ColabFold Parameter:** ~4 GB (werden automatisch heruntergeladen)
- **Gesamt:** ~8 GB

## 🎯 Nächste Schritte

Nach erfolgreicher Installation:

1. **Erste Simulation:** `amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1`
2. **Batch-Simulation:** `amp-rbc-md -f examples/batch.fasta --backend colabfold --n-replica 1`
3. **Analyse:** Ergebnisse in `output/` Verzeichnis
4. **Dokumentation:** Siehe `docs/` für weitere Details

## 💡 Tipps

- **Environment aktivieren:** Immer `conda activate amp-rbc-md` vor der Arbeit
- **GPU-Monitoring:** `watch -n 1 nvidia-smi` für GPU-Auslastung
- **Logs:** Alle Logs in `logs/` Verzeichnis
- **Backup:** Regelmäßige Backups der `output/` Verzeichnisse

---

**Fertig!** 🎉 Ihre amp-rbc-md Installation ist bereit für RBC MD Simulationen. 