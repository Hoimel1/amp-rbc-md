"""Kommandozeilen-Einstiegspunkte.

Der Setup-Eintrag `amp-rbc-md = amp_rbc_md.cli:main` verweist hierher.
`main()` ruft lediglich den `run_sim`-Click-Command aus `cli/run_sim.py` auf.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

# Dynamisch sicherstellen, dass das Verzeichnis `src/cli` importierbar ist
_cli_path = Path(__file__).resolve().parent.parent / "cli"
if str(_cli_path) not in sys.path:
    sys.path.append(str(_cli_path))


def main() -> None:  # noqa: D401
    """Delegiert an den Click-Command `run_sim`."""
    run_sim_mod = importlib.import_module("cli.run_sim")
    run_sim = getattr(run_sim_mod, "run_sim")
    run_sim.main(standalone_mode=True)
