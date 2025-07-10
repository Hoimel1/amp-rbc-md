# AMP-RBC-MD Troubleshooting

## Häufige Probleme und Lösungen

### 1. Abhängigkeitskonflikte mit ColabFold

**Problem**: ColabFold installiert inkompatible Versionen von JAX, biopython, dm-haiku, etc.

**Lösung**: Verwenden Sie `setup-fixed.sh`:
```bash
bash setup-fixed.sh
```

**Details**: Das Skript installiert:
- JAX 0.4.25 VOR ColabFold
- biopython 1.82 (kompatibel mit ColabFold 1.5.5)
- dm-haiku 0.0.10 (exakte Version für ColabFold)
- absl-py 1.4.0 (kompatibel)
- importlib-metadata 4.8.2 (kompatibel)

### 2. CUDA-Versionskonflikte

**Problem**: JAX wurde gegen CUDA 12.3 kompiliert, aber CUDA 12.1 ist installiert

**Lösung**: 
```bash
# JAX mit korrekter CUDA-Version installieren
pip uninstall jax jaxlib -y
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### 3. gmxapi Kompilierungsfehler

**Problem**: gmxapi kann nicht kompiliert werden wegen fehlender Compiler

**Lösung 1**: Über conda (empfohlen)
```bash
conda install -c conda-forge gmxapi -y
```

**Lösung 2**: Build-Tools installieren
```bash
sudo apt-get update
sudo apt-get install build-essential
pip install gmxapi
```

**Lösung 3**: Ohne gmxapi (Fallback)
```bash
bash setup-no-gmxapi.sh
```

### 4. JAX linear_util nicht verfügbar

**Problem**: JAX 0.4.38+ hat linear_util entfernt

**Lösung**: JAX 0.4.25 verwenden
```bash
pip uninstall jax jaxlib -y
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### 5. ColabFold-Abhängigkeiten inkorrekt

**Problem**: ColabFold installiert automatisch inkompatible Versionen

**Lösung**: Manuelle Installation der Abhängigkeiten
```bash
pip install absl-py==1.4.0 dm-haiku==0.0.10 biopython==1.82
pip install colabfold==1.5.5 --no-deps
```

### 6. GROMACS nicht gefunden

**Problem**: gromacs Python-Modul nicht installiert

**Lösung**: Über conda installieren
```bash
conda install -c conda-forge gromacs=2024 -y
```

### 7. NVIDIA-Treiber Probleme

**Problem**: CUDA funktioniert nicht trotz installierter Treiber

**Lösung**: Treiber neu installieren
```bash
sudo apt-get update
sudo apt-get install nvidia-driver-535
sudo reboot
```

### 8. Speicher-Probleme

**Problem**: Nicht genug RAM für ColabFold/AlphaFold

**Lösung**: 
- Mindestens 16GB RAM verwenden
- Swap-Space erweitern
- Kleinere Proteine simulieren

### 9. Disk-Space Probleme

**Problem**: Nicht genug Speicherplatz für AlphaFold-Modelle

**Lösung**:
```bash
# AlphaFold-Modelle in /tmp speichern
export COLABFOLD_DATA_DIR=/tmp/colabfold
```

### 10. Netzwerk-Probleme

**Problem**: AlphaFold-Modelle können nicht heruntergeladen werden

**Lösung**:
```bash
# Proxy-Einstellungen (falls nötig)
export HTTP_PROXY=http://proxy.company.com:8080
export HTTPS_PROXY=http://proxy.company.com:8080
```

## Setup-Skripte Übersicht

| Skript | Beschreibung | Verwendung |
|--------|-------------|------------|
| `setup-fixed.sh` | **Empfohlen** - Behebt alle bekannten Konflikte | Neue Installation |
| `setup-with-gmxapi.sh` | Mit gmxapi über conda | Wenn gmxapi benötigt |
| `setup-no-gmxapi.sh` | Ohne gmxapi | Fallback bei Kompilierungsproblemen |

## Verifikation

Nach der Installation:
```bash
python verify-installation.py
```

## Support

Bei weiteren Problemen:
1. Prüfen Sie die Logs: `tail -f ~/.conda/envs/amp-rbc-md/conda-meta/history`
2. Führen Sie `python verify-installation.py` aus
3. Erstellen Sie ein Issue mit den Ausgaben 