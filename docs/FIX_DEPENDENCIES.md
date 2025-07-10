# Vollständige Behebung von Abhängigkeitskonflikten

Diese Anleitung zeigt, wie du alle Abhängigkeitskonflikte in der amp-rbc-md Pipeline behebst.

## Problem

Die Pipeline hat mehrere Abhängigkeitskonflikte:
- NumPy/Pandas-Kompatibilitätsprobleme
- PyTorch-CUDA-Versionskonflikte
- JAX-CUDA-Installationsprobleme
- ColabFold-Abhängigkeitskonflikte

## Lösung

### Option 1: Automatische Behebung (Empfohlen)

```bash
# Repository aktualisieren
cd amp-rbc-md
git pull

# Alle Abhängigkeitskonflikte automatisch beheben
conda activate amp-rbc-md
bash fix-dependencies.sh
```

### Option 2: Schrittweise Behebung

#### Schritt 1: NumPy/Pandas-Konflikte beheben

```bash
conda activate amp-rbc-md
pip uninstall -y numpy pandas
conda install -y numpy=1.24.3 pandas=1.5.3
```

#### Schritt 2: PyTorch-CUDA-Konflikte beheben

```bash
conda activate amp-rbc-md
bash fix-pytorch-cuda.sh
```

#### Schritt 3: JAX-CUDA installieren

```bash
conda activate amp-rbc-md
bash install-jax-cuda.sh
```

#### Schritt 4: ColabFold neu installieren

```bash
conda activate amp-rbc-md
pip install --force-reinstall colabfold
```

## Verifikation

Nach der Behebung solltest du folgende Ausgabe sehen:

```bash
python -c "
import jax, jaxlib, torch, colabfold, pandas, numpy
print('=== VERIFIKATION ===')
print(f'JAX: {jax.__version__} (CUDA: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0})')
print(f'PyTorch: {torch.__version__} (CUDA: {torch.cuda.is_available()})')
print(f'ColabFold: verfügbar ({hasattr(colabfold, \"batch\")})')
print(f'Pandas: {pandas.__version__}')
print(f'NumPy: {numpy.__version__}')
"
```

Erwartete Ausgabe:
```
=== VERIFIKATION ===
JAX: 0.4.25 (CUDA: True)
PyTorch: 2.3.0+cu121 (CUDA: True)
ColabFold: verfügbar (True)
Pandas: 1.5.3
NumPy: 1.24.3
```

## Troubleshooting

### Problem: "importlib-metadata" Konflikt

```bash
pip install --upgrade importlib-metadata
```

### Problem: CUDA-Bibliotheken Konflikt

```bash
# Entferne alle CUDA-Bibliotheken und installiere neu
pip uninstall -y nvidia-cublas-cu12 nvidia-cuda-cupti-cu12 nvidia-cuda-runtime-cu12 nvidia-cudnn-cu12 nvidia-cufft-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12
bash fix-pytorch-cuda.sh
```

### Problem: ColabFold funktioniert nicht

```bash
# ColabFold komplett neu installieren
pip uninstall -y colabfold
pip install colabfold>=1.5.5
```

## Nächste Schritte

Nach erfolgreicher Behebung:

```bash
# Teste die Pipeline
amp-rbc-md --seq "AAHHIIGGLFSAGKAIHRLIRRRRR" --n-replica 1 --profile default -j 1
```

Die Pipeline sollte jetzt:
- ✅ Alle Abhängigkeitskonflikte behoben haben
- ✅ CUDA für JAX und PyTorch nutzen
- ✅ ColabFold für echte AlphaFold2-Vorhersage verwenden
- ✅ Keine Warnungen mehr anzeigen
- ✅ Optimal performant laufen

## Prävention

Um zukünftige Konflikte zu vermeiden:

1. Verwende immer die `fix-dependencies.sh` nach Updates
2. Installiere Pakete in der richtigen Reihenfolge
3. Teste die Installation nach jedem Update
4. Dokumentiere Änderungen in der environment-linux.yml 