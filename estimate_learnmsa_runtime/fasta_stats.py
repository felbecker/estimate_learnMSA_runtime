from collections import Counter
from dataclasses import dataclass
from pathlib import Path


# Whitespace and control characters, as bytes.translate delete tables
_WHITESPACE = bytes(range(33))
_NON_RESIDUES = _WHITESPACE + b"-."


@dataclass(frozen=True)
class FastaStats:
    """Summary statistics of the sequences in a FASTA file."""

    num_seqs: int
    max_len: int
    avg_len: float
    median_len: float
    total_residues: int


def read_fasta_stats(
    path: str | Path, chunk_size: int = 16 * 1024 * 1024
) -> FastaStats:
    """Computes the number of sequences and their maximum, average and median
    lengths in a single streaming pass over a FASTA file, without parsing
    records.

    Records are separated with bytes.split on "\\n>" and residues are counted
    with bytes.translate, so the work in Python is one step per sequence,
    however short the lines are. No numpy is needed.

    Gaps (- and .) and whitespace do not count towards sequence lengths, so
    aligned FASTA/A2M files yield the lengths of the unaligned sequences.

    Args:
        path: Path to the FASTA file.
        chunk_size: Number of bytes read at once. Memory usage is a small
            multiple of this.

    Returns:
        The statistics of the file.

    Raises:
        ValueError: If the file is empty or not in FASTA format.
    """
    num_headers = 0
    # Residues of the sequence that is still open
    open_len = 0
    total = 0
    # Histogram of the lengths of the complete sequences for the median
    hist: Counter[int] = Counter()
    in_header = False
    # Whether the previous chunk ended with a line break (or there was none)
    line_start = True

    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            pos = 0
            if in_header:
                newline = chunk.find(b"\n")
                if newline < 0:
                    continue
                in_header, pos = False, newline + 1
            # A leading line break turns a header at the start into "\n>"
            if pos or line_start:
                data = b"\n" + chunk[pos:]
            else:
                data = chunk
            line_start = chunk.endswith(b"\n")
            # The text of the open sequence, then one piece per header
            pieces = data.split(b"\n>")
            if num_headers == 0:
                if pieces[0].translate(None, _WHITESPACE):
                    raise ValueError(f"{path} is not a FASTA file.")
            else:
                open_len += len(pieces[0].translate(None, _NON_RESIDUES))
            if len(pieces) == 1:
                continue
            if num_headers:
                hist[open_len] += 1
                total += open_len
            num_headers += len(pieces) - 1
            # All but the last record are complete
            for record in pieces[1:-1]:
                newline = record.find(b"\n")
                n = (
                    0 if newline < 0
                    else len(record[newline:].translate(None, _NON_RESIDUES))
                )
                hist[n] += 1
                total += n
            last = pieces[-1]
            newline = last.find(b"\n")
            in_header = newline < 0
            open_len = (
                0 if in_header
                else len(last[newline:].translate(None, _NON_RESIDUES))
            )

    if num_headers == 0:
        raise ValueError(f"{path} is empty or not a FASTA file.")
    hist[open_len] += 1
    total += open_len

    return FastaStats(
        num_seqs=num_headers,
        max_len=max(hist),
        avg_len=total / num_headers,
        median_len=_histogram_median(hist),
        total_residues=total,
    )


def _histogram_median(hist: Counter[int]) -> float:
    """Median of the values counted in hist, with the same convention as
    statistics.median (mean of the two middle values for an even count)."""
    n = hist.total()
    lower_rank, upper_rank = (n - 1) // 2, n // 2
    seen = 0
    lower = None
    for value in sorted(hist):
        seen += hist[value]
        if lower is None and seen > lower_rank:
            lower = value
        if seen > upper_rank:
            return (lower + value) / 2
    raise ValueError("Empty histogram.")
