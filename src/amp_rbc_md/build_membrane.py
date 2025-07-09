from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Sequence

import yaml

from .utils import LOGGER, ensure_dir

CONFIG = yaml.safe_load(Path("config/system.yaml").read_text())


def build(profile: str, peptide_gro: str | Path, out_dir: str | Path) -> Path:
    """Baue eine Lipidmembran mithilfe von insane.py oder Dummy-Version."""
    out_dir = ensure_dir(out_dir)
    system_gro = Path(out_dir) / "membrane.gro"
    lipids = CONFIG.get(profile, CONFIG["default"])  # type: ignore[index]

    cmd: Sequence[str] = [
        "insane",
        "-f",
        str(peptide_gro),
        "-o",
        str(system_gro),
        "-l", "POPC",
        "-l", "CHOL",
        "-u", "W",
        "-sol", "W",
        "-x", "10.0",
        "-y", "10.0",
        "-z", "10.0",
    ]

    try:
        LOGGER.info("Baue Membran mit insane.py: %s", " ".join(cmd))
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        LOGGER.info("insane erfolgreich ausgeführt")
    except subprocess.CalledProcessError as e:
        LOGGER.error("insane fehlgeschlagen: %s", e)
        if e.stderr:
            LOGGER.error("insane stderr: %s", e.stderr)
        if e.stdout:
            LOGGER.error("insane stdout: %s", e.stdout)
        raise RuntimeError(f"insane konnte nicht ausgeführt werden: {e}")
    except FileNotFoundError as e:
        LOGGER.error("insane nicht gefunden: %s", e)
        raise RuntimeError(f"insane konnte nicht ausgeführt werden: {e}")
    return system_gro


__all__ = ["build"]
