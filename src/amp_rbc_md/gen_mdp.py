from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml

from .utils import LOGGER, ensure_dir

TEMPLATE_DIR = Path("mdp_templates")
CONFIG = yaml.safe_load(Path("config/md.yaml").read_text())


PROD_LENGTH_MAP = {
    "short": 500_000_000,  # 500 ns bei 1 ps steps? Wait our dt is 20 fs, 500 ns = 25e3 steps? Let's compute.
}


def _infer_nsteps(length_aa: int) -> int:
    """Gib nsteps für Produktion basierend auf Länge zurück."""
    # dt=20 fs => 50,000 steps per ns
    steps_per_ns = int(1e6 / 20)  # 50,000
    if length_aa < 20:
        return steps_per_ns * 500
    if length_aa <= 60:
        return steps_per_ns * 1000
    return steps_per_ns * 2000


def generate(
    phase: Literal["minim", "nvt", "npt", "prod"],
    out_dir: str | Path,
    *,
    length_aa: int | None = None,
    override_nsteps: int | None = None,
) -> Path:
    """Erzeuge eine MDP-Datei für die gewünschte Phase."""
    out_dir = ensure_dir(out_dir)
    tpl_path = TEMPLATE_DIR / f"{phase}.mdp"
    mdp_out = Path(out_dir) / f"{phase}.mdp"
    content = tpl_path.read_text()

    if phase == "prod":
        nsteps = (
            override_nsteps
            if override_nsteps is not None
            else _infer_nsteps(length_aa or 0)
        )
        # Ersetze placeholder nsteps in Template
        content = content.replace(
            "nsteps           = 0", f"nsteps           = {nsteps}"
        )
        LOGGER.info("nsteps für Prod auf %d gesetzt", nsteps)

    mdp_out.write_text(content)
    return mdp_out


__all__ = ["generate"]
