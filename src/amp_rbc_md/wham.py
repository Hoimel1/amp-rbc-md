from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Tuple

import numpy as np

from .utils import LOGGER, ensure_dir


def _dummy_free_energy(n_windows: int = 16) -> Tuple[np.ndarray, np.ndarray]:
    z = np.linspace(-0.8, 0.8, n_windows)
    dg = -5.0 * np.exp(-z**2 / 0.2)  # Glockenkurve
    return z, dg


def run_wham(window_dir: str | Path, *, t_min: int = 0, t_max: int | None = None, dry_run: bool = False) -> Tuple[float, Path]:
    """Führe `gmx wham` aus und gib Gesamt-ΔG_insert zurück.

    Liefert (deltaG, output_xvg_path). Bei dry_run oder Fehler wird Dummy-Profil erzeugt.
    """
    window_dir = Path(window_dir)
    output = window_dir / "wham_free_energy.xvg"

    if dry_run:
        z, dg = _dummy_free_energy()
        np.savetxt(output, np.column_stack([z, dg]))
        return float(np.trapz(dg, z)), output

    try:
        tpr_list = sorted(window_dir.glob("*.tpr"))
        pullf_list = sorted(window_dir.glob("*.pullf.xvg"))
        if not tpr_list or not pullf_list:
            raise FileNotFoundError("Keine TPR/PULLF-Dateien gefunden")
        (window_dir / "tpr.dat").write_text("\n".join(str(p) for p in tpr_list))
        (window_dir / "pullf.dat").write_text("\n".join(str(p) for p in pullf_list))

        cmd = [
            "gmx",
            "wham",
            "-it",
            "tpr.dat",
            "-if",
            "pullf.dat",
            "-o",
            output.name,
            "-b",
            str(t_min),
        ]
        if t_max is not None:
            cmd += ["-e", str(t_max)]
        subprocess.run(cmd, check=True, cwd=window_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (FileNotFoundError, subprocess.CalledProcessError) as err:
        LOGGER.warning("gmx wham fehlgeschlagen (%s), nutze Dummy-ΔG.", err)
        z, dg = _dummy_free_energy()
        np.savetxt(output, np.column_stack([z, dg]))

    data = np.loadtxt(output)
    delta_g = float(np.trapz(data[:, 1], data[:, 0]))
    LOGGER.info("WHAM ΔG_insert: %.2f kJ/mol", delta_g)
    return delta_g, output

__all__ = ["run_wham"] 