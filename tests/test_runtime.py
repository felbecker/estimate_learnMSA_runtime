import re
import subprocess
import sys
from pathlib import Path

import pytest

from estimate_learnmsa_runtime.runtime import (
    estimate_runtime, estimate_runtime_seconds, format_runtime, round_runtime
)
from estimate_learnmsa_runtime.fasta_stats import FastaStats, read_fasta_stats


ROOT = Path(__file__).parents[1]
DIR = ROOT / "tests" / "data"


@pytest.mark.parametrize(
    "seconds, minutes",
    [
        (0, 15), (1, 15), (15 * 60, 15), (15 * 60 + 1, 30), (59 * 60, 60),
        (60 * 60, 60), (61 * 60, 120), (10 * 3600 + 1, 11 * 60),
    ],
)
def test_round_runtime(seconds: float, minutes: int) -> None:
    assert round_runtime(seconds) == minutes


def test_format_runtime() -> None:
    assert format_runtime(15) == "00:15:00"
    assert format_runtime(120) == "02:00:00"
    assert format_runtime(130 * 60) == "130:00:00"


def _stats(n: int, avg: float, median: float, max_len: int) -> FastaStats:
    return FastaStats(
        num_seqs=n, max_len=max_len, avg_len=avg, median_len=median,
        total_residues=int(n * avg),
    )


def test_estimate_is_monotonic() -> None:
    base = estimate_runtime_seconds(_stats(1000, 100, 100, 200))
    assert estimate_runtime_seconds(_stats(10**6, 100, 100, 200)) > base
    assert estimate_runtime_seconds(_stats(1000, 300, 100, 600)) > base
    assert estimate_runtime_seconds(_stats(1000, 100, 300, 600)) > base


def test_estimate_ignores_max_len() -> None:
    assert estimate_runtime_seconds(_stats(1000, 100, 100, 200)) == \
        estimate_runtime_seconds(_stats(1000, 100, 100, 5000))


def test_scale() -> None:
    stats = read_fasta_stats(f"{DIR}/egf.fasta")
    assert estimate_runtime(stats) == 15
    assert estimate_runtime(stats, scale=100) > 60
    for scale in [0, -1]:
        with pytest.raises(ValueError):
            estimate_runtime(stats, scale=scale)


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "estimate_learnmsa_runtime", *args],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def test_runtime_cli_rejects_non_positive_scale() -> None:
    result = _run_cli(f"{DIR}/egf.fasta", "--scale", "0")
    assert result.returncode != 0
    assert "scale must be positive" in result.stderr


@pytest.mark.parametrize("extra_args", [[], ["--scale", "2.5"]])
def test_runtime_cli(extra_args: list[str]) -> None:
    result = _run_cli(f"{DIR}/egf.fasta", *extra_args)
    assert result.returncode == 0, result.stderr
    assert re.fullmatch(r"\d{2,}:\d{2}:00\n", result.stdout)
