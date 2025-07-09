from __future__ import annotations

import random
from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml

from .utils import LOGGER, ensure_dir

JUDGE_CFG = yaml.safe_load(Path("config/judge.yaml").read_text())


class BootstrapResult(NamedTuple):
    """ΔG und 95%-CI."""

    delta_g: float
    ci95: float


def bootstrap_delta_g(values: np.ndarray, n_boot: int = 200) -> BootstrapResult:
    """Bootstrappe ΔG über gegebene Werte."""
    boots = [np.mean(np.random.choice(values, size=len(values))) for _ in range(n_boot)]
    ci = np.percentile(boots, [2.5, 97.5])
    return BootstrapResult(float(np.mean(values)), float(ci[1] - ci[0]))


def analyse_trajectory(traj: str | Path, out_dir: str | Path) -> dict[str, float]:
    """Dummy-Analyse zur CI: generiert zufällige ΔG & thinning."""
    out_dir = ensure_dir(out_dir)
    random.seed(0)
    delta_g_values = np.random.normal(loc=-4.0, scale=1.0, size=16)
    dg, ci = bootstrap_delta_g(delta_g_values)
    thinning = abs(np.random.normal(loc=0.3, scale=0.05))

    # Speichere CSV
    csv_path = Path(out_dir) / "report.csv"
    csv_path.write_text(
        "delta_g,ci95,thinning\n{:.3f},{:.3f},{:.3f}\n".format(dg, ci, thinning)
    )
    LOGGER.info("Analyse gespeicher in %s", csv_path)
    return {"delta_g": dg, "ci95": ci, "thinning": thinning}


__all__ = [
    "bootstrap_delta_g",
    "analyse_trajectory",
]
