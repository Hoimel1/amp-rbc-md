import subprocess
from pathlib import Path

import pytest  # type: ignore

from amp_rbc_md.build_membrane import build  # type: ignore
from amp_rbc_md.martinize_wrap import martinize  # type: ignore


@pytest.fixture(autouse=True)
def no_subprocess(monkeypatch):
    """Verhindere echte Subprocess-Aufrufe."""

    def _raise(*_args, **_kwargs):
        raise FileNotFoundError

    monkeypatch.setattr(subprocess, "run", _raise)


def test_build_membrane_dummy(tmp_path: Path) -> None:
    peptide_gro = tmp_path / "peptide.gro"
    peptide_gro.write_text("dummy")
    out = build("default", peptide_gro, tmp_path)
    assert out.exists()


def test_martinize_dummy(tmp_path: Path) -> None:
    pdb = tmp_path / "model.pdb"
    pdb.write_text("HEADER")
    gro, top = martinize(pdb, tmp_path)
    assert gro.exists() and top.exists() 