# pipeline

Nextflow (DSL2) rewrite of the kleb-amr-ml pipeline. Simpler than the
original Snakemake version: one workflow engine, one language (Python)
for process logic, AWS Batch as the target executor for scale-out.

## Run locally

```bash
nextflow run main.nf
```

## Run on AWS Batch

```bash
nextflow run main.nf -profile awsbatch \
    --aws_queue <batch-queue-name> \
    --aws_region <region> \
    --aws_workdir s3://<bucket>/work
```

## Layout

- `main.nf` — entrypoint workflow
- `modules/` — one `.nf` file per pipeline stage
- `scripts/` — Python scripts called by modules (the actual logic)
