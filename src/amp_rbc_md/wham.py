from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Tuple

import numpy as np

from .utils import LOGGER, ensure_dir


def run_wham(window_dir: str | Path, *, t_min: int = 0, t_max: int | None = None, dry_run: bool = False) -> Tuple[float, Path]:
    """Führe `gmx wham` aus und gib Gesamt-ΔG_insert zurück.

    Liefert (deltaG, output_xvg_path). 
    Diese Funktion führt echte WHAM-Analysen durch und erzeugt keine Dummy-Daten.
    """
    window_dir = Path(window_dir)
    output = window_dir / "wham_free_energy.xvg"

    if dry_run:
        LOGGER.info("WHAM Dry-Run: Überspringe WHAM-Analyse")
        # Erstelle leere Ausgabedatei für Dry-Run
        output.write_text("# WHAM Dry-Run - Keine echten Daten\n")
        return 0.0, output

    # Prüfe ob GROMACS verfügbar ist
    try:
        subprocess.run(["gmx", "--version"], check=True, capture_output=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        raise RuntimeError("GROMACS ist nicht verfügbar. Bitte installieren Sie GROMACS.")

    # Suche nach WHAM-Eingabedateien
    tpr_list = sorted(window_dir.glob("*.tpr"))
    pullf_list = sorted(window_dir.glob("*.pullf.xvg"))
    
    if not tpr_list:
        raise FileNotFoundError(f"Keine TPR-Dateien in {window_dir} gefunden")
    if not pullf_list:
        raise FileNotFoundError(f"Keine PULLF-Dateien in {window_dir} gefunden")
    
    LOGGER.info("WHAM-Analyse: %d TPR-Dateien, %d PULLF-Dateien gefunden", len(tpr_list), len(pullf_list))
    
    # Erstelle Eingabedateien für WHAM
    tpr_dat = window_dir / "tpr.dat"
    pullf_dat = window_dir / "pullf.dat"
    
    tpr_dat.write_text("\n".join(str(p) for p in tpr_list))
    pullf_dat.write_text("\n".join(str(p) for p in pullf_list))
    
    # Führe WHAM aus
    cmd = [
        "gmx",
        "wham",
        "-it", "tpr.dat",
        "-if", "pullf.dat",
        "-o", output.name,
        "-b", str(t_min),
    ]
    
    if t_max is not None:
        cmd += ["-e", str(t_max)]
    
    LOGGER.info("Führe WHAM aus: %s", " ".join(cmd))
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            cwd=window_dir, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        LOGGER.info("WHAM erfolgreich ausgeführt")
        if result.stdout:
            LOGGER.debug("WHAM stdout: %s", result.stdout)
    except subprocess.CalledProcessError as e:
        LOGGER.error("WHAM fehlgeschlagen: %s", e)
        if e.stderr:
            LOGGER.error("WHAM stderr: %s", e.stderr)
        if e.stdout:
            LOGGER.error("WHAM stdout: %s", e.stdout)
        raise RuntimeError(f"WHAM-Analyse fehlgeschlagen: {e}")
    
    # Lade und analysiere WHAM-Ergebnisse
    if not output.exists():
        raise FileNotFoundError(f"WHAM-Ausgabedatei {output} wurde nicht erstellt")
    
    try:
        data = np.loadtxt(output)
        if len(data) == 0:
            raise ValueError("WHAM-Ausgabedatei ist leer")
        
        # Berechne ΔG als Integral über die Freie-Energie-Kurve
        z_coords = data[:, 0]  # Z-Koordinaten
        free_energy = data[:, 1]  # Freie Energie
        
        # Entferne NaN-Werte
        valid_mask = ~(np.isnan(z_coords) | np.isnan(free_energy))
        if not np.any(valid_mask):
            raise ValueError("Keine gültigen Daten in WHAM-Ausgabe")
        
        z_coords = z_coords[valid_mask]
        free_energy = free_energy[valid_mask]
        
        # Berechne ΔG als Differenz zwischen Minimum und Maximum
        delta_g = np.max(free_energy) - np.min(free_energy)
        
        LOGGER.info("WHAM ΔG_insert: %.2f kJ/mol", delta_g)
        return delta_g, output
        
    except Exception as e:
        LOGGER.error("Fehler beim Laden der WHAM-Daten: %s", e)
        raise RuntimeError(f"WHAM-Daten konnten nicht verarbeitet werden: {e}")


__all__ = ["run_wham"] 