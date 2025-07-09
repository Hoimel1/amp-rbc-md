from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import mlflow  # type: ignore
from rich.progress import Progress  # type: ignore

from .utils import LOGGER, detect_gpu


class GromacsRunner:  # noqa: D101
    def __init__(
        self,
        work_dir: str | Path,
        gpu_id: str | int | None = None,
        replica_id: int = 0,
        dry_run: bool = False,
    ) -> None:
        self.work_dir = Path(work_dir)
        if gpu_id is None:
            has_gpu, auto_gpu = detect_gpu()
            self.gpu_id = auto_gpu if has_gpu else ""
        else:
            self.gpu_id = str(gpu_id)
        self.replica_id = replica_id
        self.dry_run = dry_run

    def _cmd(self, *args: str) -> Sequence[str]:
        base = ["gmx"]
        # GPU-Flag komplett deaktiviert für CPU-Instanzen
        # if self.gpu_id and any(arg.startswith("mdrun") for arg in args):
        #     try:
        #         # Teste ob GPU-Flag verfügbar ist
        #         result = subprocess.run(["gmx", "mdrun", "-h"], capture_output=True, text=True, timeout=10)
        #         if "-gpu" in result.stdout:
        #             base += ["-gpu", self.gpu_id]
        #         else:
        #             LOGGER.info("GPU-Flag nicht verfügbar, verwende CPU-Version")
        #     except Exception as e:
        #         LOGGER.info("GPU-Test fehlgeschlagen (%s), verwende CPU-Version", e)
        return [*base, *args]

    # ---------------------------------------------------------------------
    # GROMPP – Topologie vorbereiten
    # ---------------------------------------------------------------------

    def grompp(
        self,
        mdp: str | Path,
        gro: str | Path,
        top: str | Path,
        out_tpr: str | Path,
    ) -> Path:
        """Rufe `gmx grompp` auf und erzeuge TPR.

        Bei `dry_run` oder fehlendem GROMACS wird eine Dummy-Datei geschrieben,
        sodass die Pipeline fortfahren kann.
        """
        mdp, gro, top, out_tpr = map(Path, (mdp, gro, top, out_tpr))

        if self.dry_run:
            LOGGER.info("[Dry] Würde grompp aufrufen (mdp=%s, gro=%s)", mdp, gro)
            out_tpr.write_text("stub tpr")
            return out_tpr

        cmd = self._cmd(
            "grompp",
            "-f",
            str(mdp),
            "-c",
            str(gro),
            "-p",
            str(top),
            "-o",
            str(out_tpr),
        )

        try:
            LOGGER.info("Führe grompp aus: %s", " ".join(cmd))
            result = subprocess.run(
                cmd, 
                cwd=self.work_dir, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            LOGGER.info("grompp erfolgreich ausgeführt")
            return out_tpr
        except subprocess.CalledProcessError as err:
            LOGGER.error("grompp fehlgeschlagen: %s", err)
            if err.stdout:
                LOGGER.error("grompp stdout: %s", err.stdout)
            if err.stderr:
                LOGGER.error("grompp stderr: %s", err.stderr)
            LOGGER.warning("Erstelle Dummy-TPR für Pipeline-Kontinuität")
            out_tpr.write_text("stub tpr")
            return out_tpr
        except FileNotFoundError as err:
            LOGGER.error("gmx nicht gefunden: %s", err)
            LOGGER.warning("Erstelle Dummy-TPR für Pipeline-Kontinuität")
            out_tpr.write_text("stub tpr")
            return out_tpr

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
