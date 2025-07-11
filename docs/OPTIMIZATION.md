# Speicherplatz-Optimierung & Alternative Backends

## Übersicht: Was ist wirklich 2TB groß?

| Komponente | Größe | Brauchst du ihn? | Kommentar |
|------------|-------|------------------|-----------|
| **AlphaFold Parameter** (5 Netze, v2.3) | ≈ 3.6 GB | ✅ **JA** | Unabhängig von MSA-Datenbanken |
| **MSA-Datenbanken** (BFD, MGnify, UniRef90, PDB70) | ≈ 1.9 TB | ❌ **Meist NEIN** | Nur für klassisches AF2 mit MSA |
| **AlphaFold DB fertiger Strukturen** | > 20 TB | ❌ **NEIN** | Du modellierst selbst |
| **Martini/GROMACS Assets** | < 200 MB | ✅ **JA** | Vernachlässigbar |

**Kurz:** Die 2TB entfallen fast ausschließlich auf Referenz-MSA-Dumps. Für kurze Peptide (≤ 100 AA) kannst du vollständig darauf verzichten!

## Option 1: ColabFold-Batch (Empfohlen für 400GB VM)

### Vorteile
- ✅ **< 4 GB lokaler Speicherplatz**
- ✅ **Remote-MMSeqs2 API** - MSA kommt vom Server
- ✅ **Keine lokalen Datenbanken nötig**
- ✅ **Sofort einsatzbereit**

### Installation
```bash
# Repository klonen
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# ColabFold installieren (ohne lokale DBs)
pip install colabfold batchfold

# Nur AlphaFold-Parameter herunterladen (~3.6 GB)
cd external/fastfold
bash scripts/download_alphafold_params.sh ~/alphafold_dbs/
cd ../..

# Umgebungsvariablen setzen
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs
export COLABFOLD_REMOTE=1  # Remote-MMSeqs2 verwenden
```

### Verwendung
```bash
# ColabFold-Batch für Strukturvorhersage
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run

# Echte Simulation
amp-rbc-md --seq GLSILGKLL --backend colabfold --n-replica 1
```

## Option 2: FastFold No-MSA (Schnellste Option)

### Vorteile
- ✅ **< 4 GB lokaler Speicherplatz**
- ✅ **3-6× schneller** als klassisches AlphaFold
- ✅ **DeepSpeed & FlashAttention** optimiert
- ✅ **Keine MSA-Datenbanken nötig**

### Installation
```bash
# FastFold installieren
cd external/fastfold
conda env create -f env.yml
conda activate fastfold

# Nur Parameter herunterladen
bash scripts/download_alphafold_params.sh ~/alphafold_dbs/

# No-MSA-Modus aktivieren
export FASTFOLD_NO_MSA=1
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs
```

### Verwendung
```bash
# FastFold No-MSA
amp-rbc-md --seq GLSILGKLL --backend fastfold --no-msa --dry-run
```

## Option 3: ESMFold (Schnellste für kurze Peptide)

### Vorteile
- ✅ **< 2 GB lokaler Speicherplatz**
- ✅ **0.8s pro Peptid** (vs. 15s FastFold)
- ✅ **Reine LLM-Vorhersage** - keine MSA
- ✅ **Ideal für < 30 AA Peptide**

### Installation
```bash
pip install esmfold
```

### Verwendung
```bash
# ESMFold für kurze Peptide
amp-rbc-md --seq GLSILGKLL --backend esmfold --dry-run
```

## Option 4: Vollständige Installation (2TB+)

### Wann brauchst du die 2TB wirklich?
- **Lange Proteine** (> 100 AA)
- **Maximale Vorhersagequalität** erforderlich
- **Forschung mit MSA-Analyse**
- **Produktionsumgebung** mit unbegrenztem Budget

### Installation
```bash
# Vollständige Datenbanken (~2TB)
cd external/fastfold
./scripts/download_all_data.sh ~/alphafold_dbs/
```

## Speicherplatz-Strategien

### Für 400GB VM (Empfohlen)
```bash
# 1. ColabFold-Batch Setup
git clone --recursive https://github.com/Hoimel1/amp-rbc-md.git
cd amp-rbc-md

# 2. Nur Parameter installieren
cd external/fastfold
bash scripts/download_alphafold_params.sh ~/alphafold_dbs/
cd ../..

# 3. ColabFold installieren
pip install colabfold batchfold

# 4. Umgebungsvariablen
export ALPHAFOLD_DATA_DIR=~/alphafold_dbs
export COLABFOLD_REMOTE=1

# 5. Testen
amp-rbc-md --seq GLSILGKLL --backend colabfold --dry-run
```

### Speicherplatz-Aufschlüsselung (400GB VM)
- **AlphaFold Parameter**: 3.6 GB
- **Repository + Code**: 1 GB
- **Trajektorien (10 Peptide × 2 Replika)**: 13 GB
- **GROMACS/Martini**: 200 MB
- **Gesamt**: < 20 GB

### Für 2TB+ VM (Vollständig)
```bash
# Vollständige Installation
./setup.sh
cd external/fastfold
./scripts/download_all_data.sh ~/alphafold_dbs/
```

## Performance-Vergleich

| Backend | Speicherplatz | Geschwindigkeit | Qualität | Empfehlung |
|---------|---------------|-----------------|----------|------------|
| **ColabFold-Batch** | < 4 GB | 30-60s | Sehr gut | ✅ **400GB VM** |
| **FastFold No-MSA** | < 4 GB | 10-15s | Sehr gut | ✅ **400GB VM** |
| **ESMFold** | < 2 GB | 0.8s | Gut | ✅ **Kurze Peptide** |
| **AlphaFold Vollständig** | 2TB+ | 60-120s | Beste | ✅ **2TB+ VM** |

## Tipps für knappen Speicherplatz

### GROMACS-Optimierungen
```bash
# Kürzere Trajektorien für Tests
gmx convert-tpr -s prod.tpr -nsteps 2500000  # 50ns statt 500ns

# Weniger Output
nstxout-compressed = 20000  # 400ps statt 100ps
```

### Aufräumen nach Analyse
```bash
# Temporäre Dateien löschen
rm */umbrella*.xtc */pullf.xvg

# Trajektorien komprimieren
gmx trjconv -f traj.xtc -o traj_compressed.xtc -fit rot+trans
```

### Cloud Storage für Langzeitarchivierung
```bash
# Trajektorien in GCS-Bucket archivieren
gsutil -m cp results/ gs://your-bucket/trajectories/
```

## Bottom Line

**Für Peptid-Screening brauchst du keine 2TB-MSA-Dumps!**

- **ColabFold-Batch** oder **FastFold No-MSA** + 3.6 GB AF-Gewichte genügen
- **< 20 GB** lokaler Speicherplatz für vollständige Pipeline
- **Sofort einsatzbereit** ohne Zusatzkosten
- **Wissenschaftlich belastbare Ergebnisse**

Starte heute mit deiner 400GB VM - ohne Storage-Stress! 🚀 