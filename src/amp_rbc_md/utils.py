from __future__ import annotations

import logging
import os
import random
from pathlib import Path
from typing import Any

import numpy as np
from rich.logging import RichHandler  # type: ignore
import subprocess

_LOGGER = None


def get_logger(name: str | None = None) -> logging.Logger:  # noqa: D401
    """Hol einen konfigurierten Logger mit Rich-Handler."""
    global _LOGGER
    if _LOGGER is None:
        logging.basicConfig(
            level=os.getenv("LOGLEVEL", "INFO"),
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(rich_tracebacks=True)],
        )
    return logging.getLogger(name) if name else logging.getLogger()


LOGGER = get_logger(__name__)


def set_seed(seed: int = 42) -> None:
    """Setze deterministische Seeds für random, numpy und (falls verfügbar) torch."""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch  # type: ignore

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ModuleNotFoundError:
        # torch ist optional
        pass


def ensure_dir(path: str | Path) -> Path:
    """Stelle sicher, dass ein Verzeichnis existiert und gib Pfad zurück."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def hash_sequence(sequence: str) -> str:
    """Liefere einen kurzen SHA1-Hash für eine Peptidsequenz."""
    import hashlib

    return hashlib.sha1(sequence.encode()).hexdigest()[:10]


# ---------------------------------------------------------------------------
# Hardware Detection
# ---------------------------------------------------------------------------


def detect_gpu() -> tuple[bool, str | None]:
    """Ermittle, ob eine NVIDIA-GPU vorhanden ist und liefere erste GPU-ID.

    Strategien:
    1. `nvidia-smi --query-gpu=index --format=csv,noheader`
    2. Fallback: Prüfe, ob GROMACS als GPU-Build kompiliert wurde (`gmx --version`).
    Liefert `(has_gpu, gpu_id)` – bei CPU-only `(False, None)`.
    """

    try:
        res = subprocess.run(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        gpu_ids = [line.strip() for line in res.stdout.splitlines() if line.strip().isdigit()]
        if gpu_ids:
            LOGGER.info("GPU erkannt via nvidia-smi: %s", gpu_ids[0])
            return True, gpu_ids[0]
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass

    try:
        res = subprocess.run(["gmx", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if "GPU support" in res.stdout:
            # Heuristik: setze auf 0
            LOGGER.info("GPU-Build von GROMACS erkannt, nehme Device 0")
            return True, "0"
    except FileNotFoundError:
        pass

    LOGGER.info("Keine GPU erkannt, wechsle auf CPU-Modus")
    return False, None


__all__: list[str] = [
    "get_logger",
    "LOGGER",
    "set_seed",
    "ensure_dir",
    "hash_sequence",
]
__all__.append("detect_gpu")
