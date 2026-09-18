# estimate-learnmsa-runtime

Prints a conservative wall time estimate for aligning a protein FASTA file
with [learnMSA](https://github.com/Gaius-Augustus/learnMSA), e.g. to set the
time limit of a cluster job. It does not need learnMSA, only numpy.

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
from estimate_learnmsa_runtime import estimate_runtime, format_runtime

minutes = estimate_runtime("proteins.fasta", scale=1.0)
print(format_runtime(minutes))
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

## Scale

Slower machines, e.g. with few CPU cores, a slow shared filesystem or a
smaller GPU, or non-default learnMSA settings (e.g. `--use_language_model`)
may need more time. `--scale` multiplies the estimate before rounding.

## Tests

```
pip install .[test]
pytest
```
