# amp-rbc-md

Batch-fähige Martini-3-Simulation roter Blutkörperchen-Peptide (≤10) mit Replika-Support, WHAM-Bootstrap-Analyse und konfigurierbarer Tox-Logik.

## Quick-Start (CPU)

```bash
# macOS / Windows (ohne GROMACS):
conda env create -f environment.yml
conda activate amp-rbc-md

# Linux + GPU (inkl. GROMACS):
conda env create -f environment-linux.yml
conda activate amp-rbc-md

amp-rbc-md run_sim --seq "GLSILGKLL" --n-replica 2 --dry-run
```

## Theorie-Hintergrund

Die Insertions-Freie-Energie (ΔG_insert) wird aus Umbrella-Sampling-Fenstern mittels WHAM bestimmt. Ein Bootstrap (N=200) liefert Konfidenzintervalle, die mit heuristischen Schwellen aus `config/judge.yaml` verglichen werden, um toxische Kandidaten zu markieren.

## Slurm-Guide (GPU)

```bash
#!/bin/bash
#SBATCH --gpus=1 --time=24:00:00
module load cuda/12.3  gromacs/2024-gpu

conda activate amp-rbc-md
export GMX_GPU=0
amp-rbc-md run_sim -f examples/batch.fasta --n-replica 3 --profile chol_high --gpu 0
```

## MLflow-Tracking

Setze `MLFLOW_TRACKING_URI`, um Parameter, Metriken und Artefakte (Trajektorie, Pore-Plots, MP4) automatisch zu speichern.

```bash
export MLFLOW_TRACKING_URI=file:$(pwd)/mlruns
```

## Projekt-Layout

Siehe Spezifikation in Issue #1 oder die Ordnerstruktur dieses Repos.
