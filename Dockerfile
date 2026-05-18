# ─────────────────────────────────────────────────────────────────────────────
# WES Pipeline Docker Image
# Base: mambaforge (conda-forge + bioconda ready)
#
# Build:
#   docker build -t wes-pipeline:1.0 .
#
# Run (interactive):
#   docker run --rm -it \
#     -v $(pwd):/workspace \
#     -w /workspace \
#     wes-pipeline:1.0 bash
#
# Run (snakemake):
#   docker run --rm \
#     -v $(pwd):/workspace \
#     -w /workspace \
#     wes-pipeline:1.0 \
#     snakemake --cores 8 --use-conda
# ─────────────────────────────────────────────────────────────────────────────

FROM condaforge/mambaforge:23.11.0-0

LABEL maintainer="suraj-chauhan-21"
LABEL description="Whole-Exome Sequencing analysis pipeline"
LABEL version="1.0"

# ── System dependencies ──────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        wget \
        unzip \
        pigz \
        default-jdk \
    && rm -rf /var/lib/apt/lists/*

# ── Create working directory ─────────────────────────────────────────────────
WORKDIR /workspace

# ── Copy environment file and install all tools ──────────────────────────────
COPY envs/environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml -n wes-pipeline \
    && conda clean -afy

# ── Make conda environment the default shell environment ─────────────────────
SHELL ["conda", "run", "-n", "wes-pipeline", "/bin/bash", "-c"]

# ── Copy workflow files ───────────────────────────────────────────────────────
COPY workflow/  /workspace/workflow/
COPY config/    /workspace/config/
COPY scripts/   /workspace/scripts/
COPY test_data/ /workspace/test_data/

# ── Set up PATH to always use the wes-pipeline conda env ─────────────────────
ENV PATH /opt/conda/envs/wes-pipeline/bin:$PATH
ENV CONDA_DEFAULT_ENV wes-pipeline

# ── Verify key tools are available ───────────────────────────────────────────
RUN bwa 2>&1 | head -3 && \
    samtools --version | head -1 && \
    gatk --version && \
    fastqc --version && \
    snakemake --version && \
    echo "All tools verified OK"

# ── Default command: show pipeline help ──────────────────────────────────────
CMD ["snakemake", "--help"]
