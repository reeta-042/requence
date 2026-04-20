FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# ── System dependencies ────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    wget curl git build-essential \
    && rm -rf /var/lib/apt/lists/*

# ── Miniconda ──────────────────────────────────────────────────
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"
RUN conda init bash

# Accept conda TOS
RUN conda config --set channel_priority flexible

# ── Bioinformatics conda environment ──────────────────────────
# abricate  – gene detection (CARD + ResFinder for Salmonella; CARD-only for E. coli)
# blast     – tBLASTn k-mer extraction
# snippy    – SNP calling against reference genome
RUN conda create -n requence python=3.10.14 -y && \
    conda install -n requence -c bioconda -c conda-forge \
        abricate=1.0.1 blast=2.16.0 snippy=4.6.0 -y

SHELL ["conda", "run", "-n", "requence", "/bin/bash", "-c"]

WORKDIR /app

# ── Environment variables ──────────────────────────────────────
ENV WORK_DIR=/app/work \
    MODELS_DIR=/app/models \
    SCRIPTS_DIR=/app/scripts \
    FEATURE_TEMPLATES_DIR=/app/feature_templates \
    # Salmonella reference + resources
    SALMONELLA_REFERENCE=/app/reference/salmonella_LT2.gbff \
    SALMONELLA_SCRIPTS=/app/scripts/salmonella \
    SALMONELLA_TEMPLATES=/app/feature_templates/salmonella \
    SALMONELLA_MODELS=/app/models/salmonella \
    # E. coli reference + resources
    ECOLI_REFERENCE=/app/reference/escherichia_coli.gbff \
    ECOLI_SCRIPTS=/app/scripts/e_coli \
    ECOLI_TEMPLATES=/app/feature_templates/e_coli \
    ECOLI_MODELS=/app/models/e_coli \
    # Shared CARD protein database
    CARD_PROTEIN_FILE=/app/card_db/card_all_proteins.fasta

# Reset SHELL for pip installs
SHELL ["/bin/bash", "-c"]

# ── Python dependencies ────────────────────────────────────────
RUN conda install -n requence -c conda-forge numpy=1.26.4 pandas=2.2.2 scikit-learn=1.4.2 -y && \
    /opt/conda/envs/requence/bin/pip install --no-cache-dir --default-timeout=3600 \
    fastapi==0.115.0 \
    uvicorn==0.32.0 \
    uvloop==0.22.1 \
    httptools==0.7.1 \
    python-multipart==0.0.12 \
    requests==2.32.3

COPY requirements.txt /app/
RUN /opt/conda/envs/requence/bin/pip install --no-cache-dir --default-timeout=3600 -r requirements.txt

# ── Layer-cache-friendly cleanup ──────────────────────────────
RUN conda clean --all -f -y && \
    rm -rf /opt/conda/pkgs/* /root/.cache/pip && \
    find /opt/conda -name "*.pyc" -delete && \
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

# Set SHELL back for final commands
SHELL ["conda", "run", "-n", "requence", "/bin/bash", "-c"]

# ── Copy project assets ────────────────────────────────────────
# scripts/salmonella/  and  scripts/e_coli/
COPY scripts/ /app/scripts/
# models/salmonella/   and  models/e_coli/
COPY models/ /app/models/
# feature_templates/salmonella/  and  feature_templates/e_coli/
COPY feature_templates/ /app/feature_templates/
# reference/salmonella_LT2.gbff  and  reference/escherichia_coli.gbff
COPY reference/ /app/reference/
COPY card_db/ /app/card_db/
COPY api/ /app/api/

# ── Permissions + runtime directories ─────────────────────────
RUN mkdir -p /app/logs /app/work /app/uploads /app/results && \
    find /app/scripts -name "*.sh" -exec chmod +x {} +

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=10s --retries=3 \
    CMD python -c \
        "import requests; r=requests.get('http://localhost:8000/health'); exit(0 if r.status_code==200 else 1)" \
    || exit 1

CMD ["conda", "run", "--no-capture-output", "-n", "requence", \
     "python", "-m", "uvicorn", "api.app:app", \
     "--host", "0.0.0.0", "--port", "8000"]