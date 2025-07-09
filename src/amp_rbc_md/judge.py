from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import yaml

from .utils import LOGGER

RULES = yaml.safe_load(Path("config/judge.yaml").read_text())

LABELS = RULES["labels"]


# noinspection PyPyDialects
class Verdict(NamedTuple):
    """(Label, Confidence)."""

    label: str
    confidence: float


def judge_metrics(delta_g: float, ci95: float, thinning: float) -> Verdict:
    """Bewerte Ergebnis gemäss YAML-Regeln."""
    score = 0.0
    if delta_g < RULES["delta_g"]["threshold_kcal_per_mol"]:
        score += 0.4
    if ci95 < RULES["delta_g"]["ci_width_kcal_per_mol"]:
        score += 0.2
    if thinning < RULES["thinning"]["threshold_nm"]:
        score += 0.4

    label = LABELS["pass"] if score >= 0.7 else LABELS["fail"]
    LOGGER.info("Bewertung: %s (%.2f)", label, score)
    return Verdict(label=label, confidence=score)


__all__ = ["judge_metrics", "Verdict"]
