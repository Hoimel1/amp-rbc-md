# Manuelle JAX-CUDA Installation auf Google Cloud VM

Falls das automatische Setup nicht funktioniert, hier die manuelle Installation der CUDA-kompatiblen JAX-Version.

## Voraussetzungen

- Ubuntu 22.04 LTS VM mit NVIDIA GPU
- Conda-Umgebung `amp-rbc-md` aktiviert
- CUDA-Treiber installiert (`nvidia-smi` funktioniert)

## Schritt 1: Alte JAX-Versionen entfernen

```bash
conda activate amp-rbc-md
pip uninstall -y jax jaxlib
```

## Schritt 2: CUDA-kompatible JAX installieren

```bash
# Installiere JAX mit CUDA 12 Support
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

## Schritt 3: Installation verifizieren

```bash
# Prüfe installierte Versionen
pip list | grep -E "(jax|jaxlib)"

# Teste CUDA-Erkennung
python -c "
import jax
import jaxlib
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')
"
```

## Erwartete Ausgabe

```
jax                                0.6.2
jaxlib                             0.6.2

JAX Version: 0.6.2
JAXlib Version: 0.6.2
JAX Devices: [gpu:0, cpu:0]
CUDA verfügbar: True
```

## Troubleshooting

### Problem: "An NVIDIA GPU may be present, but a CUDA-enabled jaxlib is not installed"

**Lösung**: Stelle sicher, dass du auf der Linux-VM bist und nicht auf einem Mac:

```bash
# Prüfe Plattform
uname -a

# Sollte Linux zeigen, nicht Darwin
```

### Problem: Falsche jaxlib-Version installiert

**Lösung**: Manuell die korrekte Version herunterladen:

```bash
# Entferne alle JAX-Versionen
pip uninstall -y jax jaxlib

# Installiere explizit die Linux-CUDA-Version
pip install jax==0.6.2 jaxlib==0.6.2+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### Problem: CUDA-Treiber nicht erkannt

**Lösung**: NVIDIA-Treiber neu installieren:

```bash
# Prüfe CUDA-Treiber
nvidia-smi

# Falls nicht verfügbar:
sudo apt-get update
sudo apt-get install -y nvidia-driver-525
sudo reboot
```

## Alternative: Conda-basierte Installation

Falls pip-Probleme auftreten:

```bash
# Conda-Installation versuchen
conda install -c conda-forge jax jaxlib cudatoolkit=12.0
```

## Teste ColabFold mit CUDA

Nach erfolgreicher JAX-CUDA-Installation:

```bash
# Teste ColabFold
python -c "
import colabfold
print('ColabFold Version:', colabfold.__version__)
print('ColabFold verfügbar:', hasattr(colabfold, 'batch'))

# Teste AlphaFold2-Vorhersage
import colabfold.batch
print('AlphaFold2-Module verfügbar:', hasattr(colabfold.batch, 'run'))
"
```

## Nächste Schritte

Nach erfolgreicher Installation:

```bash
# Starte Simulation mit GPU-Unterstützung
amp-rbc-md --seq "AAHHIIGGLFSAGKAIHRLIRRRRR" --n-replica 1 --profile default -j 1
```

Die Simulation sollte jetzt deutlich schneller laufen, da ColabFold die GPU für die AlphaFold2-Strukturvorhersage nutzen kann. 