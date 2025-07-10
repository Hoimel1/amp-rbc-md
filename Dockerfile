# amp-rbc-md Docker Image
# PyTorch 2.3.0+cu121 + nvidia-cudnn-cu12==8.9.2.26 + JAX 0.4.25
# Basierend auf NVIDIA CUDA 12.1

FROM nvidia/cuda:12.1-devel-ubuntu22.04

# Metadaten
LABEL maintainer="amp-rbc-md Team"
LABEL description="Batch-fähige Martini-3-Simulation roter Blutkörperchen-Peptide"
LABEL version="1.0.0"

# Umgebungsvariablen
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
ENV PYTHONUNBUFFERED=1

# System-Abhängigkeiten
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Miniconda installieren
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p $CONDA_DIR \
    && rm /tmp/miniconda.sh

# Arbeitsverzeichnis
WORKDIR /app

# Environment-Datei kopieren
COPY environment.yml /tmp/environment.yml

# Conda-Umgebung erstellen
RUN conda env create -f /tmp/environment.yml

# Umgebung aktivieren
SHELL ["conda", "run", "-n", "amp-rbc-md", "/bin/bash", "-c"]

# GPU-spezifische Pakete installieren
RUN pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html \
    && pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121 \
    && pip install nvidia-cudnn-cu12==8.9.2.26 \
    && pip install colabfold==1.5.5 --no-deps \
    && pip install absl-py appdirs biopython matplotlib numpy pandas py3Dmol requests tqdm dm-haiku

# Projektcode kopieren
COPY . /app

# Projekt installieren
RUN pip install -e .

# Installation testen
RUN python -c "import jax, torch, colabfold; print('✅ Installation erfolgreich')"

# Standard-Shell
SHELL ["/bin/bash", "-c"]

# Conda-Umgebung automatisch aktivieren
RUN echo "conda activate amp-rbc-md" >> ~/.bashrc

# Ports exponieren
EXPOSE 8888

# Standard-Befehl
CMD ["conda", "run", "-n", "amp-rbc-md", "bash"] 