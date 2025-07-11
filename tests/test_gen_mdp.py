from pathlib import Path

import pytest  # type: ignore

from amp_rbc_md.gen_mdp import _infer_nsteps, generate  # type: ignore


@pytest.mark.parametrize(
    "length,expected_ns",
    [
        (10, 500),
        (20, 1000),
        (59, 1000),
        (61, 2000),
    ],
)
def test_infer_nsteps(length: int, expected_ns: int) -> None:
    steps_per_ns = int(1e6 / 20)
    assert _infer_nsteps(length) == steps_per_ns * expected_ns


def test_generate_prod(tmp_path: Path) -> None:
    mdp_path = generate("prod", tmp_path, length_aa=10)
    assert "nsteps" in mdp_path.read_text()
