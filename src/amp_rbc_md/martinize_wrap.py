from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Sequence

from .utils import LOGGER, ensure_dir


def martinize(
    pdb: str | Path, out_dir: str | Path, elastic: int = 0
) -> tuple[Path, Path]:
    """Rufe martinize2 auf, um CG-Topologie & -Struktur zu erstellen.

    Für die Testumgebung wird, falls martinize2 nicht verfügbar ist, ein Dummy-Top
    und .gro erzeugt.
    """
    out_dir = ensure_dir(out_dir)
    top_path = Path(out_dir) / "system.top"
    gro_path = Path(out_dir) / "system.gro"

    cmd: Sequence[str] = [
        "martinize2",
        "-f",
        str(pdb),
        "-o",
        str(top_path),
        "-x",
        str(gro_path),
        "-ff", "martini3001",
        "-elastic" if elastic else "",
    ]
    try:
        LOGGER.info("Führe martinize2 aus: %s", " ".join(cmd))
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        LOGGER.info("martinize2 erfolgreich ausgeführt")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        LOGGER.error("martinize2 fehlgeschlagen: %s", e)
        raise RuntimeError(f"martinize2 konnte nicht ausgeführt werden: {e}")

    return gro_path, top_path


__all__: list[str] = ["martinize"]
