from amp_rbc_md.render import traj_to_mp4  # type: ignore


def test_traj_to_mp4_dummy(tmp_path):
    trr = tmp_path / "traj.trr"
    trr.write_text("stub")
    mp4 = traj_to_mp4(trr, gro=None, out_mp4=tmp_path / "out.mp4")
    assert mp4.exists()
