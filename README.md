# estimate-learnmsa-runtime

Prints a conservative wall time estimate for aligning a protein FASTA file
with [learnMSA](https://github.com/Gaius-Augustus/learnMSA), e.g. to set the
time limit of a cluster job. It needs neither learnMSA nor any other
dependency.

## Installation

```
pip install .
```

## Usage

```
estimate-learnmsa-runtime proteins.fasta [--scale SCALE]
```

or `python -m estimate_learnmsa_runtime proteins.fasta`. The estimate is
printed as `HH:MM:SS`, rounded up to 15 minute steps below one hour and to
full hours above. This format is accepted by Slurm's `--time` and
Snakemake's `runtime` resource:

```
sbatch --time=$(estimate-learnmsa-runtime proteins.fasta) job.sh
```

From Python:

```python
from estimate_learnmsa_runtime import (
    estimate_memory, estimate_runtime, format_runtime, read_fasta_stats
)

stats = read_fasta_stats("proteins.fasta")  # one pass over the file
print(format_runtime(estimate_runtime(stats, scale=1.0)))
print(estimate_memory(stats, struct=True), "GiB")
```

## How it works

The input file is scanned once (vectorized, without parsing records) for the
number of sequences N and their mean length L and median length M. Gaps
(`-`, `.`) are not counted, so aligned FASTA/A2M input works as well. The
runtime is modelled as

    overhead + a * N + b * sqrt(N) * L * M + c * N * L * M

covering per-sequence costs, training and the prediction passes (model
selection, decoding). The coefficients were fitted on 104 runs of learnMSA
with default settings on a GPU workstation (HomFam and extHomFam families
with up to 3M sequences) and multiplied by a safety factor of 2.5 that covers
all of them, also in leave-one-out fits.

## Memory

```
estimate-learnmsa-memory proteins.fasta [--struct] [--language_model]
```

prints a conservative memory estimate in GiB with Slurm's unit suffix:

```
sbatch --mem=$(estimate-learnmsa-memory proteins.fasta --struct) job.sh
```

From Python, `estimate_memory(stats, struct=..., language_model=...)` returns
the GiB as an integer. Memory grows
linearly with the total number of residues R:

    1.25 * (4 GiB + R * (32 B + 64 B per language model))

rounded up to 4 GiB steps, at least 8 GiB. `--struct` (3Di states predicted
with ProstT5) and `--language_model` (ProtT5 embeddings) each count as one
language model. The coefficients cover the peak proportional set size (PSS)
of 36 extHomFam runs with default settings on a GPU workstation, with and
without `--struct`, by at least 1.23x. The language model budget also covers
the brief copy of the per-residue cache that learnMSA makes when many
sequences are near-identical, which those runs did not exercise. ProtT5 was
not measured and gets the ProstT5 budget.

## Scale

Slower machines, e.g. with few CPU cores, a slow shared filesystem or a
smaller GPU, or non-default learnMSA settings (e.g. `--use_language_model`)
may need more time. `--scale` multiplies the estimate before rounding.

## Tests

```
pip install .[test]
pytest
```
