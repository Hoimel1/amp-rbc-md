from amp_rbc_md.pore import detect_pore  # type: ignore


def test_detect_pore_dummy(tmp_path):
    traj = tmp_path / "traj.trr"
    traj.write_text("dummy")
    result = detect_pore(traj)
    assert isinstance(result, bool)
