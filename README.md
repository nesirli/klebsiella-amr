# Klebsiella pneumoniae AMR prediction

**Live app:** <https://klebsiella-amr.nasirnesirli.com>
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

**<https://klebsiella-amr.nasirnesirli.com>** — username `demo`,
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
docker build -t klebsiella-amr-app -f app/Dockerfile .
docker run --rm -p 8501:8501 \
  -e AMR_APP_USERNAME=demo -e AMR_APP_PASSWORD=klebsiella2026 \
  klebsiella-amr-app
```

See `app/README.md` for the Dokploy settings.

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

## Limits and known issues

- The pipeline runs Kraken2 on every sample, but nothing uses the result, so a
  contaminated isolate can still reach the training set. (The web app *does*
  check the species — see below.)
- `build_sequences.py` keys genes on the third field of the AMRFinder FASTA
  header, which its docstring claims matches `build_features.py`. It does not:
  378 of 686 keys are coordinate-stamped frameshift entries like
  `cirA_@del_240_241_24418_24421_0_STOP`, unique to one isolate, where the
  tabular matrix has 451 clean `Element symbol` values. DNABERT-2 therefore
  sees hundreds of near-singleton pseudo-genes. Fixing it needs re-tuning.
- `make report` deletes `data/amr/`, so adding DNABERT-2 afterwards means
  re-running AMRFinder over the kept assemblies (a couple of hours, no
  re-assembly). Run `make dnabert` before `make report`.
- Each model chooses its decision threshold from the training data. This helps
  a lot when a model is badly calibrated, but it can be slightly too
  optimistic.

---

The rest of this file is for anyone changing the pipeline.

## Project structure

```text
.
├── Makefile                  # pipeline orchestration
├── config.yaml               # user-editable configuration
├── modules/                  # Makefile includes, one per pipeline stage
│   ├── config.mk             # config.yaml -> make variables
│   ├── metadata.mk           # train/test split
│   ├── pipeline.mk           # per-sample rules: download -> trim -> assemble -> AMR
│   ├── features.mk           # gene matrices
│   ├── models.mk             # tuning + training rules
│   └── batch.mk              # batch driver, MultiQC, report
├── envs/
│   ├── env-ml.yml            # amr: Python + ML dependencies
│   └── env-bio.yml           # bioinfo: fastp, spades, kraken2, quast, amrfinderplus, multiqc, seqtk
├── metadata.csv              # input metadata (semicolon-delimited)
├── scripts/
│   ├── common.py             # shared data loading, metrics, artifacts, Optuna driver
│   ├── mlp.py                # the MLP shared by train_nn.py and tune_nn.py
│   ├── metadata.py
│   ├── download_reads.py
│   ├── build_features.py
│   ├── build_sequences.py
│   ├── export_config.py
│   ├── export_manifest.py    # freezes feature layout + thresholds for the app
│   ├── summarize.py
│   ├── train_{xgboost,lightgbm,nn,dnabert}.py
│   └── tune_{xgboost,lightgbm,nn,dnabert}.py
├── app/                      # the deployed Streamlit service
│   ├── app.py                # UI: login, upload, species gate, results
│   ├── predict.py            # AMRFinder -> features -> rule + model calls
│   ├── artifacts/            # manifest, served models, species sketch
│   ├── environment.yml       # micromamba spec for the image
│   ├── Dockerfile
│   └── README.md             # deployment settings
├── tests/
├── reference/                # downloaded reference databases
├── data/                     # transient per-sample files (auto-deleted)
└── results/                  # features, models, reports (kept)
```

Every model script shares one contract, implemented in `scripts/common.py`: read a
gene matrix, fit one binary R-vs-S classifier for one antibiotic, and emit the same
six artifacts (model, params, metrics, predictions, importance CSV, importance
plot). Only the estimator and its importance method differ between models. An
antibiotic without enough labelled data to fit still emits all six, marked
`"skipped"`, so a partial dataset cannot fail the run.

Tuners hand hyperparameters to trainers as
`{"hyperparameters": {...}, "tuning": {...}}` — search bookkeeping is kept out of
the block that gets splatted into the estimator.

DNABERT-2 needs its own tuner (`tune_dnabert.py`) rather than the shared
`common.run_tuning`, because it consumes nucleotide strings where the others
consume a float matrix. It reuses the same fold logic and writes the same
handoff file. Its search is deliberately small — 5 trials × 2 folds on a capped
subsample — because every fold refits 117M parameters, so a trial costs minutes
rather than milliseconds.

### Decision thresholds

Each model picks its threshold by maximising balanced accuracy on the
**training** split (`common.choose_threshold`) rather than assuming 0.5, and
records it in the metrics. A fixed 0.5 assumes calibrated probabilities over a
balanced prior, and neither holds here: DNABERT-2 ranked amikacin at 0.980 ROC
AUC yet called all 178 test isolates susceptible, because 74-vs-228 training
labels kept every probability under 0.5.

Caveat worth knowing: the threshold is fitted on predictions the model has
already partly memorised, so it is biased toward the training distribution. For
badly calibrated models this is a large net win; for well-calibrated ones it
can cost a little (XGBoost amikacin balanced accuracy 0.985 → 0.956). Choosing
it on cross-validated out-of-fold predictions instead would remove the bias and
is the obvious next improvement.

## Choosing a split

`splits.mode` decides how train and test are separated, and it changes what a
score means:

- **`temporal`** (default) trains on isolates collected up to `train_cutoff`
  and tests on `test_years`. This is the deployment question — does a model
  fitted on old isolates work on new ones? — and it is the harder, more honest
  number.
- **`random`** holds out a stratified `test_fraction` at random, so both sides
  share a distribution. Use it as a *control*.

Running both is how you tell "the features carry no signal" apart from "the
signal does not survive a shift between eras". On this dataset ceftazidime
scores 0.62–0.68 test ROC AUC under the temporal split and 0.94 under the
random one, from the same features and the same models: its training years are
47% resistant and its test years 85%, because the split is confounded with
which studies deposited isolates when. Amikacin, by contrast, is ~1.0 either
way. A drug that collapses only under the temporal split has a population
problem, not a feature problem.

## Interpretability

Each model emits a `gene,importance` CSV (most important first) and a plot. The
column is the same across models so they can be compared side by side; what fills
it differs:

- **XGBoost / LightGBM**: mean |SHAP| over the training rows, plotted as a SHAP beeswarm
- **MLP**: permutation importance — drop in test ROC AUC when a gene column is shuffled
- **DNABERT-2**: occlusion importance — drop in test P(resistant) when a gene is removed from the pooled embedding

The `report` target aggregates per-model metrics and the top genes per
model/antibiotic into `results/report/summary.json`. Models that were not run
(DNABERT-2 is optional) are simply absent from it.

## Make semantics this pipeline depends on

Non-obvious rules that took a full-scale run to find. Changing them will
silently re-run days of work, or silently delete it.

- **`.SECONDARY:`** (in `Makefile`). GNU make 4.4 deletes files it considers
  intermediate when it exits — here that meant every assembly, AMR table,
  fastp JSON and Kraken report vanishing mid-run, while their grouped-target
  siblings survived. `.SECONDARY` with no prerequisites keeps all of them.
  `.DELETE_ON_ERROR` still removes the target of a recipe that failed.
- **Transient FASTQs are order-only.** Reads, trimmed and downsampled files are
  deleted per sample by design. As *normal* prerequisites, make plans to
  recreate the deleted files, treats that plan as brand-new, and so declares
  every kept artifact stale — sending the analyze phase back to re-download and
  re-assemble every finished sample. After the `|`, they are still built when
  something genuinely needs them. The kept files carry the real timestamps.
- **`download_reads.py` is order-only too.** Editing the downloader must not
  invalidate reads already on disk. To genuinely refetch a sample, delete its
  reads.
- **MultiQC depends on the `cleaned` markers, not on `qc/.done`.** Reaching for
  a fastp JSON drags in the grouped rule that produced it, and that group also
  lists the trimmed FASTQs which cleanup has deleted — an incomplete group, so
  make re-runs fastp for every sample. A `cleaned` marker proves QC, Kraken2
  and QUAST all finished and belongs to no group.
- **A capped run must sample, not truncate.** `MAX_SAMPLES` uses a stratified
  draw; NCBI metadata arrives grouped by submitting study, so `.head(n)` once
  produced a split that was 18 R and 1 S for a drug whose pool is 51% R.

## Failure handling

Transient problems must not be mistaken for bad data, because the pipeline's
response to bad data is permanent: a failed sample is *retired* and drops out
of the run.

- `download_reads.py` separates "cannot reach ENA" (DNS failures, resets,
  timeouts, HTTP 403/429/5xx) from "this accession is bad" (404, malformed
  filereport). Network failures do not consume the retry budget; they back off
  and wait up to 30 minutes. A brief outage once retired 57 healthy samples.
- A truncated-but-advancing transfer counts as progress, not failure. ENA
  closes long transfers early, and one 186 MB file needed several resumed
  passes; charging each pass against the retry budget retired a sample that was
  downloading correctly.
- Downloads are verified by decompressing them, not by size. One FASTQ arrived
  at exactly its advertised `Content-Length` and still ended mid-stream — only
  fastp noticed, three steps later, and the sample was retired for it.
- Retired samples are listed in `results/metadata/retired.txt`. **Re-running
  the same command recovers them**, because finished samples are skipped via
  their `data/cleaned/*_cleaned` markers. Check that file when a run ends.

## Serving

`app/` is a self-contained Streamlit service: upload an assembled genome, get a
resistance profile. `make app-artifacts` then `make app`, or build the
Dockerfile in that folder. `app/README.md` covers deployment; the notes below
are about why it is built the way it is.

**Why it shows a rule and a model side by side.** Neither wins everywhere.
Measured on the same test set, balanced accuracy:

| antibiotic | rule | xgboost | lightgbm |
|---|---|---|---|
| amikacin | 0.703 | 0.956 | 0.959 |
| ciprofloxacin | 0.508 | 0.761 | 0.741 |
| ceftazidime | 0.509 | 0.561 | 0.531 |
| meropenem | **0.928** | 0.821 | 0.856 |

The rule is "resistant if AMRFinder reports any gene of that drug class". It
wins for meropenem, because carbapenemase carriage is a clean determinant. It
is useless for ciprofloxacin and ceftazidime — specificity 0.02 and 0.04 —
because *K. pneumoniae* intrinsically carries quinolone efflux genes
(`oqxA`/`oqxB`) and a chromosomal beta-lactamase, so "has a gene of this class"
is true of nearly every isolate. The models learn to ignore those and weigh the
QRDR mutations that matter. AMRFinder runs either way: it is the models' input,
not an alternative to them.

**The manifest is the correctness boundary.** `results/models/manifest.json`
freezes the gene columns *in training order*, plus each model's fitted
threshold. Serving must index features by the manifest, never by the uploaded
sample — a misaligned vector produces a confident number rather than an error.
Genes the sample has that training lacked are dropped and reported; a large
count there means the AMRFinder database has moved since training.
`test_predict.py` asserts a served prediction reproduces the training
prediction, which is what catches an alignment regression.

Model paths inside the manifest are stored **relative to the manifest file**,
so the same document serves the repo tree (`results/models/`) and the image
(`/app/artifacts/`) with no rewriting; `predict.py` resolves them against the
manifest's own directory.

Build the manifest **before** `make report`: it needs the feature matrices,
which that target deletes.

The image cannot use the `python:slim` + `uv` pattern of the other portfolio
services, because AMRFinderPlus is a bioconda package that shells out to BLAST+
and HMMER. micromamba installs both stacks from one spec. `predict.py` calls
`amrfinder` directly when it is on PATH (the container) and falls back to
`conda run -n bioinfo` otherwise (a local checkout).

**Species gate.** The models were fitted only on *K. pneumoniae*, so the app
refuses anything else rather than scoring it. Three checks run before
AMRFinder, cheapest first:

1. Is it nucleotide FASTA at all? Catches uploaded FASTQ reads and protein
   FASTA, which would otherwise fail obscurely several tools later.
2. Is the assembly 3-9 Mb? Catches single contigs, read files and eukaryotic
   sequence without running anything expensive.
3. Mash distance to `artifacts/kpneumoniae_ref.msh`, a 1000-hash sketch of 20
   assemblies from this dataset (164 KB).

The 0.05 cutoff was calibrated, not guessed: 60 held-out *K. pneumoniae*
assemblies scored at most **0.028** (median 0.009), while *E. coli* K-12 --
the hardest realistic confusion, same family -- scored **0.075**. The cutoff
sits in that gap and is close to the conventional 95% ANI species boundary.
Members of the *K. pneumoniae* complex (*K. variicola*, *K. quasipneumoniae*)
fall under it and are accepted, which is intended: they share the AMR genetics.

It **fails closed**. A missing sketch or a missing `mash` binary rejects the
upload rather than skipping the check, so a broken image cannot silently start
accepting any organism. `docker-entrypoint.sh` also asserts the sketch exists
at start-up.

This is the gate the pipeline itself still lacks -- Kraken2 runs over every
training sample but nothing consumes its verdict, so a contaminated isolate can
still reach the training set.

**Not served:** the MLP (adds little over the trees) and DNABERT-2 (450 MB per
antibiotic, slow, and worse on every drug).

**The password is a gate, not authentication** — one shared secret, no
accounts, no lockout, no audit trail. It is compared with `hmac.compare_digest`
and read from `.streamlit/secrets.toml` (gitignored) or `$AMR_APP_PASSWORD`.
Serve over HTTPS or it crosses the wire in clear text, and put a real proxy in
front before widening access. Uploads go to a temporary directory that is
deleted when the analysis finishes, and are never logged.

## License

MIT — see [LICENSE](LICENSE).
