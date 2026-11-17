"""Conservative memory estimate for aligning a FASTA file with learnMSA, e.g.
to configure a cluster job. Run standalone with

    estimate-learnmsa-memory proteins.fasta [--struct] [--language_model]

which prints the estimate in GiB with Slurm's unit suffix, for example:

    sbatch --mem=$(estimate-learnmsa-memory proteins.fasta --struct) job.sh
"""
import argparse
import math

from .fasta_stats import FastaStats, read_fasta_stats


# Fitted on the peak proportional set size (PSS, every physical page counted
# once, as a cgroup memory limit does) of the process tree of 36 runs with
# default settings on a GPU workstation (2026-11, learnMSA 3.0.0, extHomFam
# families with up to 1.3M sequences and 152M residues), plus one older
# --struct run of ABC_tran (3M sequences, 498M residues, 37 GiB). Memory grows
# linearly with the total number of residues R:
# - default: 2.1 GiB + 29 B * R (largest deviation from the fit 0.3 GiB).
# - --struct (3Di states predicted with ProstT5): 3.3 GiB + 74 B * R. The 3Di
#   logits of every residue are kept in memory (float16 * 20 = 40 B).
# In those runs nearly every sequence was its own representative. With many
# near-identical sequences, learnMSA predicts only the representatives and
# copies their logits to the other sequences, which briefly holds the cache
# about twice (80 B * R plus 16 B index per copied residue). Each language
# model is therefore budgeted 64 B * R on top of the default 32 B * R.
# --use_language_model (ProtT5) was not measured. It takes the same path with
# a smaller cache (16 dimensions, 32 B) and gets the same budget.
# The fit is scaled by a safety factor. The estimate covers every run used for
# fitting by at least 1.23x (p450 with --struct), most by 1.5x to 2.3x.
_BASE_GIB = 4.0
_BYTES_PER_RESIDUE = 32
_BYTES_PER_RESIDUE_LANGUAGE_MODEL = 64
_SAFETY = 1.25
# rounding step and minimum of the estimate in GiB
_STEP_GIB = 4
_MIN_GIB = 8


def estimate_memory(
    stats: FastaStats, struct: bool = False, language_model: bool = False
) -> int:
    """Conservative estimate of the peak memory in GiB needed to align the
    sequences described by stats with learnMSA, e.g. to configure a cluster
    job. Get stats with read_fasta_stats, once for all estimates of a file.

    Args:
        stats: Statistics of the input file.
        struct: Whether learnMSA predicts 3Di states with ProstT5 (--struct
            without a structure file).
        language_model: Whether learnMSA uses ProtT5 embeddings
            (--use_language_model).

    Returns:
        The estimated memory in GiB, rounded up to 4 GiB steps, at least 8.
    """
    per_residue = (
        _BYTES_PER_RESIDUE
        + _BYTES_PER_RESIDUE_LANGUAGE_MODEL * (int(struct) + int(language_model))
    )
    gib = _SAFETY * (_BASE_GIB + stats.total_residues * per_residue / 2**30)
    return max(_MIN_GIB, math.ceil(gib / _STEP_GIB) * _STEP_GIB)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="estimate-learnmsa-memory",
        description="Print a conservative memory estimate for aligning a "
            "FASTA file with learnMSA using default settings, in GiB with "
            "the suffix G accepted by Slurm's --mem.",
    )
    parser.add_argument("input_file", help="FASTA file to align.")
    parser.add_argument(
        "--struct",
        action="store_true",
        help="learnMSA predicts 3Di states with ProstT5 (--struct).",
    )
    parser.add_argument(
        "--language_model",
        action="store_true",
        help="learnMSA uses ProtT5 embeddings (--use_language_model).",
    )
    args = parser.parse_args(argv)
    try:
        gib = estimate_memory(
            read_fasta_stats(args.input_file), args.struct, args.language_model
        )
    except (ValueError, OSError) as e:
        parser.error(str(e))
    print(f"{gib}G")


if __name__ == "__main__":
    main()
