# Klebsiella AMR prediction — app

Streamlit app: upload an assembled *Klebsiella pneumoniae* genome, get a
predicted resistance pattern. Deployed to Coolify from this folder's
Dockerfile.

## Layout

```
app/
├── app.py                 Streamlit UI (password gate, upload, results)
├── predict.py             AMRFinder -> feature vector -> rule + model calls
├── artifacts/             baked into the image, produced by `make app-artifacts`
│   ├── manifest.json      gene columns in training order + fitted thresholds
│   ├── xgboost/*.json
│   └── lightgbm/*.txt
├── environment.yml        micromamba spec (bioinformatics + Python stack)
├── Dockerfile
└── docker-entrypoint.sh   start-up checks
```

## Why micromamba and not python-slim + uv

AMRFinderPlus is a bioconda package that shells out to BLAST+ and HMMER. It is
not on PyPI, so the usual `python:3.12-slim` + `uv sync` image cannot build it.
micromamba installs the bioinformatics tools and the Python packages from one
lockable spec. The image is correspondingly large (~2 GB, mostly BLAST and the
AMRFinder database).

The ~250 MB AMRFinder database is baked in at build time rather than fetched on
start-up: otherwise every cold start waits on NCBI, and an NCBI outage takes
the service down.

## Build and run locally

```bash
make app-artifacts                 # from the repo root, refreshes artifacts/
docker build -t klebsiella-amr-app app/
docker run --rm -p 8501:8501 \
  -e AMR_APP_USERNAME=demo -e AMR_APP_PASSWORD=klebsiella2026 \
  klebsiella-amr-app
```

Then open <http://localhost:8501>.

To run without Docker (uses the repo's conda environments):

```bash
make app
```

## Coolify deployment

Live at <https://nasirnesirli.com/portfolio/klebsiella-amr/app>

| Setting | Value |
|---|---|
| Build pack | Dockerfile |
| Base directory | `/app` |
| Dockerfile location | `/app/Dockerfile` |
| Port | `8501` |

Environment variables:

| Variable | Value | Notes |
|---|---|---|
| `AMR_APP_USERNAME` | `demo` | Defaults to `demo` if unset. |
| `AMR_APP_PASSWORD` | `klebsiella2026` | Required. Set it here, never in the image. |
| `BASE_URL_PATH` | `/portfolio/klebsiella-amr/app` | Must match the proxy path exactly, no trailing slash. |

`BASE_URL_PATH` is passed to `--server.baseUrlPath`, which is what makes
Streamlit emit correctly prefixed asset **and websocket** URLs behind the
subpath. Without it the page loads but the websocket 404s and the app hangs on
"Please wait…".

The image sets `--server.enableCORS=false` and
`--server.enableXsrfProtection=false`, which Streamlit needs behind a reverse
proxy that terminates TLS.

Coolify's healthcheck can use the built-in `HEALTHCHECK`, which already accounts
for `BASE_URL_PATH` (`{base}/_stcore/health`).

## Updating the models

Artifacts are committed so the image builds without access to `results/`. After
retraining:

```bash
make app-artifacts
git add app/artifacts && git commit -m "Refresh served models"
```

Then redeploy. `make app-artifacts` must run before `make report`, which
deletes the feature matrices the manifest is generated from.

## Species gate

Uploads are checked before anything expensive runs: nucleotide FASTA, then a
3-9 Mb size band, then Mash distance to `artifacts/kpneumoniae_ref.msh`. The
0.05 cutoff was calibrated against this dataset -- 60 held-out *K. pneumoniae*
scored at most 0.028, *E. coli* K-12 scored 0.075.

The check fails closed: a missing sketch or missing `mash` rejects the upload
instead of skipping verification.

## Notes

- Research use only. Not a diagnostic tool.
- The password is a single shared secret with no accounts or lockout. Coolify
  terminates TLS in front of it; do not expose the container port directly.
- Uploads are written to a temporary directory and deleted when the analysis
  finishes. They are not logged.
- Predictions for ceftazidime are flagged in the UI as unreliable on recent
  isolates — see the main `README.md`.
