# Klebsiella pneumoniae AMR prediction

**Live app:** <https://nasirnesirli.com/portfolio/klebsiella-amr/app>
(username `demo`, password `klebsiella2026`)

This project predicts antibiotic resistance in *Klebsiella pneumoniae* from
sequencing data. It downloads bacterial genomes, finds resistance genes in
them, and trains machine learning models to predict which antibiotics will
work.

The whole pipeline runs with one command.

## What it does

1. Reads sample metadata and splits it into a training set and a test set.
2. Downloads the sequencing reads for each sample from ENA.
3. Cleans the reads, checks their quality, and assembles a genome.
4. Finds antibiotic resistance genes in each genome with AMRFinderPlus.
5. Builds a table of which genes each sample has.
6. Trains four models: XGBoost, LightGBM, a small neural network, and
   DNABERT-2 (optional).
7. Writes a quality report and a summary of the results.

The pipeline deletes large temporary files as soon as it no longer needs them,
so it can run on a laptop.

## Requirements

- [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba)
- GNU Make. On macOS, install it with `brew install make` and use `gmake`.
- About 10 GB of disk space for the reference databases.

## Setup

Create the two conda environments:

```bash
make setup
```

## How to run

Start with a small test run of 5 samples per split:

```bash
gmake metadata MAX_SAMPLES=5
gmake -j 3 all BATCH_SIZE=3
```

To use all samples, set `max_samples: -1` in `config.yaml` and run:

```bash
caffeinate -i gmake -j 3 all BATCH_SIZE=3
```

A full run takes several days. You can stop it at any time. If you run the
same command again, it continues from where it stopped.

To also train the optional DNABERT-2 model, run `gmake dnabert` **before**
`gmake all`. It takes about 14 hours.

## Configuration

Edit `config.yaml`. The most useful settings are:

| Setting | What it does |
|---|---|
| `antibiotics` | Which antibiotics to model |
| `splits.mode` | `temporal` (train on old samples, test on new) or `random` |
| `max_samples` | Limit the number of samples; `-1` means use all |
| `resources` | Threads and memory for each tool |

**Important:** set `resources` to match your computer. The pipeline runs
`BATCH_SIZE` samples at the same time, and each one uses the threads you set.
On a 12-core laptop with 18 GB of RAM, use `BATCH_SIZE=3` with 3 threads per
tool and `spades_memory: 4`. Higher values will slow the run down or fill the
memory.

## Results

From a run with 1184 samples. The numbers are ROC AUC on the test set (1.0 is
perfect, 0.5 is random guessing):

| Antibiotic | XGBoost | LightGBM | Neural net | DNABERT-2 |
|---|---|---|---|---|
| Amikacin | 1.00 | 1.00 | 1.00 | 0.95 |
| Ciprofloxacin | 0.85 | 0.87 | 0.85 | 0.84 |
| Ceftazidime | 0.62 | 0.62 | 0.68 | 0.54 |
| Meropenem | 0.87 | 0.88 | 0.87 | 0.78 |

Three things to know when you read these numbers:

- **The models find the right genes.** For amikacin they rank `armA` and
  `aac(6')-Ib` highest. For meropenem they rank `blaKPC-3` and `ompK36`
  highest. These are the genes that really cause resistance.
- **Ceftazidime looks bad, but the data is the problem.** With
  `splits.mode: random`, the same models reach 0.94. The training years and
  the test years contain very different samples, so the model cannot transfer
  between them.
- **DNABERT-2 is not worth its cost here.** It is as good as the other models
  on two antibiotics and worse on two, but it needs about 14 hours instead of
  a few seconds.

## Web app

**<https://nasirnesirli.com/portfolio/klebsiella-amr/app>** — username `demo`,
password `klebsiella2026`.

The app takes an assembled genome and shows the predicted resistance pattern.

First set a password. Copy the example file and change the value:

```bash
cp .streamlit/secrets.toml.example .streamlit/secrets.toml
```

Then build the model files for the app and start it:

```bash
make app-artifacts
make app
```

The app lives in the `app/` folder. It can also run in Docker, which is how it
is deployed:

```bash
docker build -t klebsiella-amr-app app/
docker run --rm -p 8501:8501 \
  -e AMR_APP_USERNAME=demo -e AMR_APP_PASSWORD=klebsiella2026 \
  klebsiella-amr-app
```

See `app/README.md` for the Coolify settings.

The app first checks that the file really is a *Klebsiella pneumoniae* genome.
It rejects other species, sequencing reads, and protein files, because the
models only know *Klebsiella* and would give a confident but wrong answer for
anything else.

For each antibiotic the app shows two answers:

- **AMRFinder rule** — is a resistance gene of this class present?
- **Models** — XGBoost and LightGBM predictions.

The two do not always agree, and the app marks it when they differ. A
disagreement means you should look at the genes yourself.

Neither answer is better everywhere. The rule is better for meropenem. The
models are much better for ciprofloxacin, because almost every *Klebsiella*
carries a quinolone gene, so the rule says "resistant" for nearly every sample.

**Security:** the password is one shared word for a small group. It is not real
user management. Run the app over HTTPS, and put a proper login in front of it
if more people need access.

## Limits

- The pipeline runs Kraken2 to check species, but it does not remove samples
  that fail this check.
- `build_sequences.py` uses the wrong key for some genes, so DNABERT-2 sees
  extra genes that should not be there.
- Each model chooses its decision threshold from the training data. This helps
  a lot when a model is badly calibrated, but it can be slightly too
  optimistic.

## More information

`NOTES.md` explains the pipeline in more detail: how the Makefile dependencies
work, how the pipeline handles network errors and failed samples, and the full
results with their caveats. Read it before you change the Makefile.

## License

MIT
