import subprocess
from pathlib import Path

import pytest  # type: ignore

from amp_rbc_md.utils import detect_gpu  # type: ignore
from amp_rbc_md.gmx_runner import GromacsRunner  # type: ignore


@pytest.mark.parametrize("nvidia_available", [True, False])
def test_detect_gpu(monkeypatch, nvidia_available: bool) -> None:
    def fake_run(cmd, **kwargs):  # noqa: D401
        if cmd[0] == "nvidia-smi":
            if nvidia_available:
                class Res:  # noqa: D401
                    stdout = "0\n"
                    returncode = 0
                return Res()
            raise FileNotFoundError
        raise FileNotFoundError
    monkeypatch.setattr(subprocess, "run", fake_run)
    has_gpu, gpu_id = detect_gpu()
    assert has_gpu == nvidia_available
    assert (gpu_id == "0") == nvidia_available


def test_runner_auto_cpu(monkeypatch, tmp_path: Path) -> None:
    # Simuliere kein GPU und kein nvidia-smi
    def fake_run(*_a, **_kw):
        raise FileNotFoundError
    monkeypatch.setattr(subprocess, "run", fake_run)
    runner = GromacsRunner(tmp_path, gpu_id=None, dry_run=True)
    assert runner.gpu_id == "" 