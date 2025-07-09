from amp_rbc_md.judge import judge_metrics, LABELS  # type: ignore


def test_judge_pass() -> None:
    label, conf = judge_metrics(delta_g=-5.0, ci95=0.5, thinning=0.1)
    assert label == LABELS["pass"]
    assert conf >= 0.7


def test_judge_fail() -> None:
    label, conf = judge_metrics(delta_g=2.0, ci95=2.0, thinning=1.0)
    assert label == LABELS["fail"]
    assert conf < 0.7 