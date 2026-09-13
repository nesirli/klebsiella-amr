# Klebsiella AMR prediction app.
#
# Unlike the other portfolio services this cannot use python-slim + uv:
# AMRFinderPlus is a bioconda package that shells out to BLAST+ and HMMER, so
# the whole stack is installed with micromamba instead.
#
# Built with the repository root as context:  docker build -f app/Dockerfile .
FROM mambaorg/micromamba:1.5.8 AS base

COPY --chown=$MAMBA_USER:$MAMBA_USER app/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Bake the ~250 MB AMRFinderPlus database into the image. It could be fetched
# at start-up instead, but then every cold start waits on NCBI and an outage
# takes the service down.
ENV AMRFINDER_DB=/opt/amrfinder_db

# The micromamba image runs as an unprivileged user, which cannot create
# directories under /opt; make the target writable before populating it.
USER root
RUN mkdir -p ${AMRFINDER_DB} && chown -R $MAMBA_USER:$MAMBA_USER ${AMRFINDER_DB}
USER $MAMBA_USER

RUN amrfinder_update -d ${AMRFINDER_DB}

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER app/app.py app/predict.py ./
COPY --chown=$MAMBA_USER:$MAMBA_USER app/artifacts ./artifacts
COPY --chown=$MAMBA_USER:$MAMBA_USER app/docker-entrypoint.sh /app/docker-entrypoint.sh

USER root
RUN chmod +x /app/docker-entrypoint.sh
USER $MAMBA_USER

ENV MANIFEST=/app/artifacts/manifest.json

# BASE_URL_PATH lets a reverse proxy serve this behind a subpath (e.g.
# /portfolio/klebsiella-amr/app) so Streamlit generates correctly prefixed
# asset and websocket URLs. Leave unset to serve from the domain root (the
# default in the container).
ENV BASE_URL_PATH=""

# The password. Set it in the platform's environment, never in the image.
ENV AMR_APP_PASSWORD=""

EXPOSE 8501

# HEALTHCHECK execs directly (no entrypoint), so the micromamba env is not
# activated: use the absolute path to the env's python, not the bare name.
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD /opt/conda/bin/python -c "import os,sys,urllib.request; p=os.environ.get('PORT','8501'); b=os.environ.get('BASE_URL_PATH','').rstrip('/'); sys.exit(0 if urllib.request.urlopen(f'http://localhost:{p}{b}/_stcore/health', timeout=3).status==200 else 1)"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/docker-entrypoint.sh"]

# enableCORS=false + enableXsrfProtection=false are required behind Dokploy's
# reverse proxy; the proxy terminates TLS and forwards the websocket.
CMD ["sh", "-c", "streamlit run app.py \
    --server.address=0.0.0.0 \
    --server.port=${PORT:-8501} \
    --server.baseUrlPath=\"$BASE_URL_PATH\" \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false \
    --browser.gatherUsageStats=false"]
