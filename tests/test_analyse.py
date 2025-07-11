import numpy as np
import pytest
from pathlib import Path

from amp_rbc_md.analyse import bootstrap_delta_g, analyse_membrane_thinning


def test_bootstrap_delta_g():
    """Teste Bootstrap-Funktion mit echten Daten."""
    # Erstelle Test-Daten
    test_values = np.array([-3.5, -4.2, -3.8, -4.0, -3.9])

    result = bootstrap_delta_g(test_values, n_boot=100)

    assert isinstance(result.delta_g, float)
    assert isinstance(result.ci95, float)
    assert result.ci95 > 0  # CI sollte positiv sein
    assert -5.0 < result.delta_g < -2.0  # Realistischer Bereich


def test_analyse_membrane_thinning_missing_mda():
    """Teste Thinning-Analyse ohne MDAnalysis."""
    with pytest.raises(
        RuntimeError, match="MDAnalysis ist für Thinning-Analyse erforderlich"
    ):
        analyse_membrane_thinning(Path("nonexistent.trr"))


def test_analyse_membrane_thinning_missing_file():
    """Teste Thinning-Analyse mit fehlender Trajektorie."""
    # Mock MDAnalysis import
    import sys
    from unittest.mock import MagicMock

    mock_mda = MagicMock()
    sys.modules["MDAnalysis"] = mock_mda

    with pytest.raises(FileNotFoundError):
        analyse_membrane_thinning(Path("nonexistent.trr"))


def test_analyse_trajectory_missing_file(tmp_path):
    """Teste Analyse mit fehlender Trajektorie."""
    from amp_rbc_md.analyse import analyse_trajectory

    with pytest.raises(FileNotFoundError):
        analyse_trajectory(Path("nonexistent.trr"), tmp_path)
