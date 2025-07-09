from __future__ import annotations

import random
from pathlib import Path
from typing import Tuple

from .utils import LOGGER

try:
    import MDAnalysis as mda  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    mda = None  # type: ignore


DEFAULT_Z_CUTOFF = 0.5  # nm
DEFAULT_MIN_DURATION_NS = 5.0
DT_FS = 20.0  # Zeitschritt aus mdp


def detect_pore(
    traj: str | Path,
    gro: str | Path | None = None,
    *,
    z_cutoff_nm: float = DEFAULT_Z_CUTOFF,
    min_duration_ns: float = DEFAULT_MIN_DURATION_NS,
) -> bool:
    """Detektiere transienten Porenwassereintritt.

    Heuristik: ≥1 Wasser innerhalb |z| ≤ z_cutoff für mindestens min_duration.
    """
    traj = Path(traj)

    if mda is None:
        raise RuntimeError("MDAnalysis ist nicht verfügbar. Bitte installieren Sie es mit: pip install MDAnalysis")

    if not traj.exists():
        raise FileNotFoundError(f"Trajektorie {traj} existiert nicht")

    if gro is None:
        gro_guess = traj.with_suffix(".gro")
        if gro_guess.exists():
            gro = gro_guess
        else:
            raise FileNotFoundError("Gro-Datei für Pore-Analyse nicht gefunden")

    u = mda.Universe(str(gro), str(traj))  # type: ignore[arg-type]
    water = u.select_atoms("resname W* and name OW")

    # Analyse entlang der Trajektorie
    frames_inside = 0
    min_frames = int(min_duration_ns * 1e6 / DT_FS)
    for ts in u.trajectory:  # pragma: no cover
        if any(abs(atom.position[2] / 10) <= z_cutoff_nm for atom in water):  # Å → nm
            frames_inside += 1
            if frames_inside >= min_frames:
                return True
        else:
            frames_inside = 0
    return False


__all__ = ["detect_pore"] 