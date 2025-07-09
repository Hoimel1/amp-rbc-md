from pathlib import Path

from amp_rbc_md.gmx_runner import GromacsRunner  # type: ignore


def test_gmx_runner_dry(tmp_path: Path) -> None:
    runner = GromacsRunner(work_dir=tmp_path, dry_run=True)
    tpr = tmp_path / "system.tpr"
    tpr.write_text("stub")

    trr_path = runner.run(tpr)
    # Pfad wird zurückgegeben, Datei aber nicht geschrieben.
    assert trr_path.name.endswith(".trr")
    assert not trr_path.exists() 