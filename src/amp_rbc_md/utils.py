from __future__ import annotations

import logging
import os
import random
from pathlib import Path
from typing import Any

import numpy as np
from rich.logging import RichHandler  # type: ignore

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


__all__: list[str] = [
    "get_logger",
    "LOGGER",
    "set_seed",
    "ensure_dir",
    "hash_sequence",
]
