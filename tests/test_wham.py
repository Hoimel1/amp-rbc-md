import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from amp_rbc_md.wham import run_wham


def test_run_wham_dry_run(tmp_path):
    """Teste WHAM Dry-Run."""
    delta_g, xvg = run_wham(tmp_path, dry_run=True)

    assert delta_g == 0.0
    assert xvg.exists()
    assert "WHAM Dry-Run" in xvg.read_text()


def test_run_wham_missing_gmx(tmp_path):
    """Teste WHAM ohne GROMACS."""
    with patch("subprocess.run") as mock_run:
        mock_run.side_effect = FileNotFoundError("gmx not found")

        with pytest.raises(RuntimeError, match="GROMACS ist nicht verfügbar"):
            run_wham(tmp_path, dry_run=False)


def test_run_wham_missing_files(tmp_path):
    """Teste WHAM ohne Eingabedateien."""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)

        with pytest.raises(FileNotFoundError, match="Keine TPR-Dateien"):
            run_wham(tmp_path, dry_run=False)


def test_run_wham_success(tmp_path):
    """Teste erfolgreiche WHAM-Analyse."""
    # Erstelle Mock-Dateien
    tpr_file = tmp_path / "test.tpr"
    pullf_file = tmp_path / "test.pullf.xvg"
    tpr_file.write_text("dummy tpr")
    pullf_file.write_text("dummy pullf")

    # Mock WHAM-Ausgabe
    wham_output = tmp_path / "wham_free_energy.xvg"
    wham_data = """# Z-Koordinate    Freie-Energie
-0.8    -2.5
-0.4    -4.0
0.0     -5.0
0.4     -4.0
0.8     -2.5"""
    wham_output.write_text(wham_data)

    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0, stdout="WHAM success")

        delta_g, output = run_wham(tmp_path, dry_run=False)

        assert delta_g > 0  # ΔG sollte positiv sein
        assert output.exists()


def test_run_wham_empty_output(tmp_path):
    """Teste WHAM mit leerer Ausgabe."""
    # Erstelle Mock-Dateien
    tpr_file = tmp_path / "test.tpr"
    pullf_file = tmp_path / "test.pullf.xvg"
    tpr_file.write_text("dummy tpr")
    pullf_file.write_text("dummy pullf")

    # Leere WHAM-Ausgabe
    wham_output = tmp_path / "wham_free_energy.xvg"
    wham_output.write_text("")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)

        with pytest.raises(ValueError, match="WHAM-Ausgabedatei ist leer"):
            run_wham(tmp_path, dry_run=False)
