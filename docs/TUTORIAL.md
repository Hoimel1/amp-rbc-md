# amp-rbc-md – Tutorial

Dieses Tutorial führt dich in ~15 Minuten durch einen vollständigen Workflow:

1. Installation (CPU-Modus)
2. Einzel-Simulation (Dry-Run)
3. Batch + Replika (parallel)
4. Analyse & mlflow
5. HPC-Ausführung (Slurm + Snakemake)

---

## 1  Installation

```bash
conda env create -f environment.yml
conda activate amp-rbc-md
pip install -e .
```

Für GPU + echte GROMACS-Runs verwende `environment-linux.yml` (Linux-Host oder Docker).

---

## 2  Einzel-Simulation (Dry-Run)

```bash
amp-rbc-md run_sim --seq "GLSILGKLL" --n-replica 1 --dry-run
```

Output:
* `results/<hash>/rep0/report.csv` → Metriken (ΔG, CI95, thinning, pore)
* `traj.mp4` → gerendertes Video (Dummy ohne MDAnalysis)

---

## 3  Batch + Replika (parallel)

```bash
amp-rbc-md run_sim -f examples/batch.fasta --n-replica 3 -j 4 --resume
```

* Parallelisiert Replika auf 4 Worker-Threads.
* `--resume` überspringt bereits gelaufene Replika.

---

## 4  Analyse & mlflow

```bash
export MLFLOW_TRACKING_URI=file:$(pwd)/mlruns
amp-rbc-md run_sim --seq "VLAVILVLL" --n-replica 2
mlflow ui
```

Im Browser unter `localhost:5000` findest du Runs inklusive MP4-Artefakten.

---

## 5  HPC (Slurm + Snakemake)

```bash
snakemake --profile workflow/profile/slurm -j 64
```

Das Slurm-Profil erzeugt pro Regel einen `sbatch`. Passe Ressourcen in `workflow/profile/slurm/config.yaml` an.

---

Happy Simulating 🎉 