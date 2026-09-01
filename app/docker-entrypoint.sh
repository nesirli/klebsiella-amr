#!/bin/sh
set -e

# Fail loudly at start-up rather than when the first user uploads a genome.
# Both of these are baked in at build time; a miss means the image is broken.
[ -f /app/artifacts/manifest.json ] || {
  echo "FATAL: /app/artifacts/manifest.json missing -- run 'make app-artifacts' before building" >&2
  exit 1
}
[ -f /app/artifacts/kpneumoniae_ref.msh ] || {
  echo "FATAL: species reference sketch missing -- the app would accept any organism" >&2
  exit 1
}
[ -d "${AMRFINDER_DB}/latest" ] || {
  echo "FATAL: AMRFinderPlus database missing at ${AMRFINDER_DB}/latest" >&2
  exit 1
}

exec "$@"
