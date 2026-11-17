import subprocess
import sys
from pathlib import Path

import pytest

from estimate_learnmsa_runtime.fasta_stats import FastaStats, read_fasta_stats
from estimate_learnmsa_runtime.memory import estimate_memory


ROOT = Path(__file__).parents[1]
DIR = ROOT / "tests" / "data"


def _residues(total: int) -> FastaStats:
    return FastaStats(
        num_seqs=1, max_len=1, avg_len=1.0, median_len=1.0,
        total_residues=total,
    )


@pytest.mark.parametrize(
    "residues, struct, language_model, gib",
    [
        (10**6, False, False, 8),           # small inputs get the minimum
        (10**6, True, False, 8),
        (151_805_742, True, False, 24),     # p450, peak PSS 17.8 GiB
        (136_407_952, True, False, 24),     # HATPase_c, peak PSS 11.3 GiB
        (497_946_013, False, False, 24),    # ABC_tran
        (497_946_013, True, False, 64),     # ABC_tran, 37 GiB, 32 GiB failed
        (497_946_013, False, True, 64),     # ProtT5 has the ProstT5 budget
        (497_946_013, True, True, 100),
    ],
)
def test_estimate_memory(
    residues: int, struct: bool, language_model: bool, gib: int
) -> None:
    assert estimate_memory(
        _residues(residues), struct=struct, language_model=language_model
    ) == gib


def test_estimate_memory_is_monotonic() -> None:
    small, large = _residues(10**8), _residues(2 * 10**8)
    assert estimate_memory(large) > estimate_memory(small)
    assert estimate_memory(small, struct=True) > estimate_memory(small)
    assert estimate_memory(small, struct=True, language_model=True) > \
        estimate_memory(small, struct=True)


def test_estimate_memory_from_file() -> None:
    assert estimate_memory(read_fasta_stats(f"{DIR}/egf.fasta"), struct=True) == 8


@pytest.mark.parametrize(
    "extra_args, expected", [([], "8G"), (["--struct", "--language_model"], "8G")]
)
def test_memory_cli(extra_args: list[str], expected: str) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "estimate_learnmsa_runtime.memory",
         f"{DIR}/egf.fasta", *extra_args],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == f"{expected}\n"
