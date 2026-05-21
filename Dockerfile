# Antibody Humanization Advisor — production container
#
# Build:
#   docker build -t humanization-advisor:latest .
#
# Run (with OASis DB bind-mounted from host — don't bake the 23 GB into the image):
#   docker run -d --name humanization \
#     --restart always \
#     -p 5000:5000 \
#     -v /path/on/host/OASis_9mers_v1.db:/data/OASis_9mers_v1.db:ro \
#     humanization-advisor:latest
#
# OR use docker-compose.yml in the repo root.

FROM python:3.12-slim AS base

# OS-level deps: hmmer is required by ANARCI; build-essential needed by some
# wheels. git is occasionally pulled in by sapiens/biophi installs.
RUN apt-get update && apt-get install -y --no-install-recommends \
      hmmer \
      build-essential \
      git \
      curl \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Non-root user for the running app
RUN useradd --create-home --shell /bin/bash app

WORKDIR /app

# Install Python deps first (better layer caching)
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir -r requirements.txt \
 && pip install --no-cache-dir \
      flask==2.3.3 \
      openpyxl==3.1.5 \
      gunicorn==21.2.0

# Pre-cache Sapiens model weights so the first user request is fast.
# (Sapiens pulls ~75 chunks from HuggingFace on first call. Doing it in the
# build means the container starts ready-to-serve and is air-gappable.)
RUN python3 -c "from sapiens import predict_scores; \
  predict_scores('QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR', chain_type='H')" \
  || echo "(warning: Sapiens pre-warm failed; first user request will pay the download cost)"

# App code last so source changes don't bust the deps layer
COPY pipeline   /app/pipeline
COPY evaluation /app/evaluation
COPY web        /app/web
COPY run_web.py /app/run_web.py

# Where the bind-mounted OASis DB will appear
ENV OASIS_DB_PATH=/data/OASis_9mers_v1.db

# Drop privileges
RUN chown -R app:app /app
USER app

EXPOSE 5000

# Use gunicorn (production WSGI) instead of Flask's dev server.
# One worker because the pipeline holds a single Sapiens/CamSol/ABodyBuilder2
# state and is CPU/GPU-bound — concurrency wouldn't help.
CMD ["gunicorn", "--workers=1", "--threads=1", "--timeout=600", \
     "--bind=0.0.0.0:5000", "--access-logfile=-", "web.app:app"]
