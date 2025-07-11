from amp_rbc_md.umbrella import generate_windows  # type: ignore


def test_generate_windows(tmp_path):
    paths = generate_windows(tmp_path, n_windows=4, spacing_nm=0.1, k=500)
    assert len(paths) == 4
    assert "pull_00.mdp" in {p.name for p in paths}
    content = paths[0].read_text()
    assert "500" in content and "pull" in content
