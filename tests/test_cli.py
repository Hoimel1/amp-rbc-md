from __future__ import annotations

import sys
from pathlib import Path

from click.testing import CliRunner  # type: ignore

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from cli.run_sim import run_sim  # type: ignore  # noqa: E402
import subprocess


def _fake_run(self, tpr_file):  # noqa: D401
    return tpr_file.with_suffix(".trr")


def test_cli_dry(monkeypatch):
    monkeypatch.setattr("amp_rbc_md.gmx_runner.GromacsRunner.run", _fake_run)
    runner = CliRunner()
    result = runner.invoke(
        run_sim, ["--seq", "ACDEFGHIK", "--n-replica", "1", "--dry-run"]
    )
    assert result.exit_code == 0


def test_cli_help() -> None:
    runner = CliRunner()
    result = runner.invoke(run_sim, ["--help"])
    assert result.exit_code == 0
    assert "Peptidsequenz" in result.output
