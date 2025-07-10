# Saubere Neuinstallation ohne Abhängigkeitskonflikte

Das Problem: ColabFold und TensorFlow haben unvereinbare Abhängigkeitsanforderungen, die zu komplexen Konflikten führen.

## Lösung: Komplette Neuinstallation

### Option 1: Automatische Neuinstallation (Empfohlen)

```bash
# Repository aktualisieren
cd amp-rbc-md
git pull

# Komplette Neuinstallation
bash reset-environment.sh
```

### Option 2: Manuelle Neuinstallation

```bash
# 1. Alte Umgebung entfernen
conda deactivate
conda env remove -n amp-rbc-md -y

# 2. Neue Umgebung erstellen
conda create -n amp-rbc-md python=3.10 -y
conda activate amp-rbc-md

# 3. Basis-Pakete installieren
conda install -y gromacs=2024 biopython click pyyaml tqdm matplotlib mlflow moviepy rich pytest pytest-cov black flake8 mypy isort

# 4. Kompatible NumPy/Pandas
conda install -y numpy=1.24.3 pandas=1.5.3

# 5. JAX-CUDA installieren (korrekte Version)
pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# 6. ColabFold installieren
pip install colabfold==1.5.5

# 7. PyTorch separat installieren
pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 8. Projekt installieren
pip install -e .
```

## Warum diese Lösung funktioniert

### Das Problem
- **ColabFold 1.5.5** benötigt: `numpy<2.0.0`, `pandas<2.0.0`, `importlib-metadata<5.0.0`
- **TensorFlow 2.13.1** benötigt: `numpy<=1.24.3`, `typing-extensions<4.6.0`
- **JAX 0.6.2** benötigt: `numpy>=1.26`

Diese Anforderungen schließen sich gegenseitig aus!

### Die Lösung
1. **Kein TensorFlow**: Wir verwenden nur ColabFold für AlphaFold2
2. **Kompatible Versionen**: NumPy 1.24.3, Pandas 1.5.3
3. **Separate Installation**: JAX und PyTorch werden separat installiert
4. **Saubere Reihenfolge**: JAX vor ColabFold, um Konflikte zu vermeiden

## Erwartete Ausgabe

```
=== INSTALLATION ERFOLGREICH ===
JAX Version: 0.4.25
JAXlib Version: 0.4.25+cuda12.cudnn89
JAX Devices: [cuda(id=0)]
JAX CUDA verfügbar: True
PyTorch Version: 2.3.0+cu121
PyTorch CUDA verfügbar: True
PyTorch CUDA Version: 12.1
PyTorch GPU: NVIDIA L4
ColabFold verfügbar: True
Pandas Version: 1.5.3
NumPy Version: 1.24.3
=== ALLE TESTS ERFOLGREICH ===
```

## Vorteile

✅ **Keine Konflikte**: Alle Pakete sind kompatibel  
✅ **GPU-Unterstützung**: JAX und PyTorch nutzen CUDA  
✅ **AlphaFold2**: ColabFold funktioniert korrekt  
✅ **Stabilität**: Keine unerwarteten Fehler  
✅ **Performance**: Optimale GPU-Nutzung  

## Nächste Schritte

Nach erfolgreicher Installation:

```bash
# Teste die Pipeline
amp-rbc-md --seq "AAHHIIGGLFSAGKAIHRLIRRRRR" --n-replica 1 --profile default -j 1
```

Die Pipeline sollte jetzt stabil und schnell laufen! 🚀 