import numpy as np

from amp_rbc_md.analyse import bootstrap_delta_g  # type: ignore


def test_bootstrap_delta_g() -> None:
    values = np.array([-5.0, -4.5, -3.8, -4.1])
    dg, ci = bootstrap_delta_g(values, n_boot=50)
    assert dg < 0
    assert ci > 0 