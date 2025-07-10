FROM nvidia/cuda:12.3.0-runtime-ubuntu22.04

# System-Pakete
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        ca-certificates \
        && rm -rf /var/lib/apt/lists/*

# Miniconda installieren
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/conda.sh && \
    bash /tmp/conda.sh -b -p $CONDA_DIR && rm /tmp/conda.sh && \
    $CONDA_DIR/bin/conda clean -afy
ENV PATH=$CONDA_DIR/bin:$PATH

# Abhängigkeiten
COPY environment-linux.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy
ENV CONDA_DEFAULT_ENV=amp-rbc-md
ENV PATH=$CONDA_DIR/envs/amp-rbc-md/bin:$PATH

# Projektcode
WORKDIR /workspace
COPY . /workspace
RUN pip install -e ./ && \
    pip install martinize2==2.2.0 insane

# Standardbefehl
ENTRYPOINT ["amp-rbc-md"]
CMD ["--help"] 