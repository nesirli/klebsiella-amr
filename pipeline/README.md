# pipeline

Nextflow (DSL2) rewrite of the kleb-amr-ml pipeline. Simpler than the
original Snakemake version: one workflow engine, one language (Python)
for process logic, AWS Batch as the target executor for scale-out.

From an NCBI Pathogen Detection metadata export it: splits samples into
train/test by collection year, downloads reads (ENA), QC-trims (fastp),
classifies (kraken2), downsamples to a target coverage (seqtk),
assembles (SPAdes), QCs the assembly (quast), detects AMR genes
(amrfinder), calls SNPs vs. a reference (snippy), aggregates reports
(MultiQC), then builds an AMR-gene feature matrix and trains one
XGBoost resistance classifier per antibiotic.

## Run locally

```bash
nextflow run main.nf                    # all samples
nextflow run main.nf --max_samples 4    # 4 per split (train + test) for a quick run
```

Requires Docker (colima or Docker Desktop). The kraken2 (~6GB),
amrfinder (~250MB) and reference-genome downloads are cached under
`reference/` via `storeDir`, so they're fetched once.

## Run on AWS Batch

```bash
nextflow run main.nf -profile awsbatch \
    --aws_queue <batch-queue-name> \
    --aws_region <region> \
    --aws_workdir s3://<bucket>/work
```

Infrastructure (compute environment, job queue, IAM, S3 bucket) is
provisioned by the Terraform config in `terraform/`.

## Layout

- `main.nf` — entrypoint workflow (stage wiring + `output {}` publishing)
- `nextflow.config` — params, `standard`/`awsbatch` profiles, resource pool
- `modules/` — one `.nf` file per pipeline stage
- `bin/` — self-contained Python scripts (PEP 723 inline deps) that
  modules call; Nextflow auto-stages `bin/` onto every task's PATH
- `terraform/` — AWS Batch infrastructure as code
