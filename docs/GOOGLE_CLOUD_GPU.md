# Google Cloud GPU Setup mit CUDA-Unterstützung

Diese Anleitung zeigt, wie du die amp-rbc-md Pipeline auf Google Cloud mit NVIDIA L4 GPU und CUDA-Unterstützung einrichtest.

## Voraussetzungen

- Google Cloud Account mit GPU-Quota
- Ubuntu 22.04 LTS VM mit NVIDIA L4 GPU
- Mindestens 16 GB RAM
- 100 GB Festplattenspeicher

## 1. VM erstellen

```bash
# Erstelle VM mit NVIDIA L4 GPU
gcloud compute instances create amp-rbc-gpu \
  --zone=us-central1-a \
  --machine-type=n1-standard-4 \
  --maintenance-policy=TERMINATE \
  --accelerator="type=nvidia-l4,count=1" \
  --image-family=ubuntu-2204-lts \
  --image-project=ubuntu-os-cloud \
  --boot-disk-size=100GB \
  --metadata="install-nvidia-driver=true"
```

## 2. VM verbinden und Setup starten

```bash
# SSH zur VM
gcloud compute ssh amp-rbc-gpu --zone=us-central1-a

# Repository klonen
git clone https://github.com/your-username/amp-rbc-md.git
cd amp-rbc-md

# Setup ausführen
bash setup-gcp-gpu.sh
```

## 3. CUDA-Treiber prüfen

Nach dem Setup sollten die NVIDIA-Treiber korrekt installiert sein:

```bash
nvidia-smi
```

Erwartete Ausgabe:
```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 525.105.17   Driver Version: 525.105.17   CUDA Version: 12.0     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util.  Compute M. |
|===============================+======================+======================|
|   0  NVIDIA L4           Off | 00000000:00:04.0 Off |                    0 |
| N/A   45C    P8    70W /  72W |      0MiB /  23028MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
```

## 4. JAX-CUDA-Unterstützung testen

```bash
# Aktiviere conda-Umgebung
conda activate amp-rbc-md

# Teste JAX-CUDA
python -c "
import jax
import jaxlib
print(f'JAX Version: {jax.__version__}')
print(f'JAXlib Version: {jaxlib.__version__}')
print(f'JAX Devices: {jax.devices()}')
print(f'CUDA verfügbar: {len([d for d in jax.devices() if d.platform == \"gpu\"]) > 0}')
"
```

Erwartete Ausgabe:
```
JAX Version: 0.6.2
JAXlib Version: 0.6.2
JAX Devices: [gpu:0, cpu:0]
CUDA verfügbar: True
```

## 5. ColabFold mit CUDA testen

```bash
# Teste ColabFold
python -c "
import colabfold
print('ColabFold Version:', colabfold.__version__)
print('ColabFold verfügbar:', hasattr(colabfold, 'batch'))
"
```

## 6. Erste Simulation starten

```bash
# Teste mit kurzer Sequenz
amp-rbc-md --seq "AAHHIIGGLFSAGKAIHRLIRRRRR" --n-replica 1 --profile default -j 1
```

## Troubleshooting

### Problem: "Unable to initialize backend 'cuda'"

**Lösung**: JAX-CUDA-Version neu installieren:

```bash
conda activate amp-rbc-md
bash install-jax-cuda.sh
```

### Problem: CUDA-Treiber nicht erkannt

**Lösung**: NVIDIA-Treiber manuell installieren:

```bash
sudo apt-get update
sudo apt-get install -y nvidia-driver-525
sudo reboot
```

### Problem: ColabFold-Abhängigkeitskonflikte

**Lösung**: Umgebung neu erstellen:

```bash
conda env remove -n amp-rbc-md
conda env create -f environment-linux.yml
conda activate amp-rbc-md
bash install-jax-cuda.sh
pip install -e .
```

## Performance-Optimierung

### GPU-Monitoring

```bash
# GPU-Nutzung überwachen
watch -n 1 nvidia-smi

# Detaillierte GPU-Info
nvidia-smi -q
```

### Memory-Optimierung

Für große Batches oder lange Sequenzen:

```bash
# Reduziere Batch-Größe
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

# Starte Simulation mit angepassten Parametern
amp-rbc-md -f batch.fasta --n-replica 2 --profile default -j 2
```

## Kosten-Optimierung

- Verwende Preemptible VMs für nicht-kritische Simulationen
- Nutze Spot-Instances für Batch-Verarbeitung
- Setze automatisches Herunterfahren nach Simulationen

```bash
# Automatisches Herunterfahren nach 2 Stunden Inaktivität
sudo shutdown -h +120
```

## Nächste Schritte

1. **Batch-Simulationen**: Verwende `examples/batch.fasta` für mehrere Peptide
2. **HPC-Deployment**: Nutze Snakemake für Cluster-Ausführung
3. **MLflow-Tracking**: Aktiviere Experiment-Tracking für Reproduzierbarkeit
4. **Docker-Deployment**: Verwende `mydockerhub/amp-rbc-md:latest` für Container-basierte Ausführung 