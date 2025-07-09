from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import mlflow  # type: ignore
from rich.progress import Progress  # type: ignore

from .utils import LOGGER


class GromacsRunner:  # noqa: D101
    def __init__(
        self,
        work_dir: str | Path,
        gpu_id: str | int | None = None,
        replica_id: int = 0,
        dry_run: bool = False,
    ) -> None:
        self.work_dir = Path(work_dir)
        self.gpu_id = str(gpu_id) if gpu_id is not None else os.getenv("GMX_GPU", "")
        self.replica_id = replica_id
        self.dry_run = dry_run

    def _cmd(self, *args: str) -> Sequence[str]:
        base = ["gmx"]
        if self.gpu_id:
            base += ["-gpu_id", self.gpu_id]
        return [*base, *args]

    def run(self, tpr_file: str | Path) -> Path:  # noqa: D401
        """Starte GROMACS-MD-Lauf und logge via mlflow."""
        tpr_file = Path(tpr_file)
        out_prefix = tpr_file.with_suffix("").name
        output_trr = self.work_dir / f"{out_prefix}.trr"

        if self.dry_run:
            LOGGER.info("[Dry] Würde GROMACS mit %s starten", tpr_file)
            return output_trr

        cmd = self._cmd("mdrun", "-s", str(tpr_file), "-deffnm", out_prefix)
        LOGGER.info("Starte GROMACS: %s", " ".join(cmd))
        with Progress() as progress:
            task = progress.add_task("Simuliere", start=False)
            proc = subprocess.Popen(
                cmd,
                cwd=self.work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            progress.start_task(task)
            for line in proc.stdout or []:
                progress.advance(task)
                sys.stdout.write(line)
            return_code = proc.wait()
        if return_code != 0:
            raise RuntimeError("GROMACS lief mit Fehlercode %s", return_code)

        mlflow.log_artifact(str(output_trr))
        return output_trr


__all__ = ["GromacsRunner"]
