import sys
from pathlib import Path

# Ermögliche Direktlauf ohne Installation
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from amp_rbc_md.utils import hash_sequence, set_seed  # type: ignore # noqa: E402


def test_hash_sequence_deterministic() -> None:
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert hash_sequence(seq) == hash_sequence(seq)


def test_set_seed_reproducibility() -> None:
    import random
    import numpy as np

    set_seed(42)
    a = random.random()
    na = np.random.rand()

    set_seed(42)
    b = random.random()
    nb = np.random.rand()

    assert a == b
    assert na == nb
