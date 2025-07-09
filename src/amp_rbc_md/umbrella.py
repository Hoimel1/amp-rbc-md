from __future__ import annotations

from pathlib import Path
from typing import List

from .utils import LOGGER, ensure_dir

PULL_TEMPLATE = Path("mdp_templates/pull.mdp").read_text()
DEFAULT_K = 1000  # kJ mol^-1 nm^-2


def generate_windows(out_dir: str | Path, n_windows: int = 16, spacing_nm: float = 0.1, k: int = DEFAULT_K) -> List[Path]:
    """Erzeuge Pull-MDP-Dateien für Umbrella-Sampling.

    Windows werden symmetrisch um z0=0 verteilt: -((n-1)/2)*spacing … +
    """
    out_dir = ensure_dir(out_dir)
    center_index = (n_windows - 1) / 2
    mdp_paths: list[Path] = []
    for idx in range(n_windows):
        offset = (idx - center_index) * spacing_nm
        mdp_content = PULL_TEMPLATE.format(k=k, z0=f"{offset:.2f}")
        path = out_dir / f"pull_{idx:02d}.mdp"
        path.write_text(mdp_content)
        mdp_paths.append(path)
    LOGGER.info("%d Umbrella-Fenster erzeugt", n_windows)
    return mdp_paths

__all__ = ["generate_windows"] 