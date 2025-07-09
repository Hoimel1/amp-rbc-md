from pathlib import Path

from amp_rbc_md.wham import run_wham  # type: ignore


def test_run_wham_dummy(tmp_path: Path):
    delta_g, xvg = run_wham(tmp_path, dry_run=True)
    assert xvg.exists()
    assert delta_g < 0 