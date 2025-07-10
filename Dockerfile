# Dockerfile für amp-rbc-md mit exakten Versionen
# PyTorch 2.3.0+cu121 + nvidia-cudnn-cu12==8.9.2.26 + JAX 0.4.25

FROM nvidia/cuda:12.1-devel-ubuntu22.04

# Setze Umgebungsvariablen
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Installiere System-Abhängigkeiten
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Installiere Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p $CONDA_DIR \
    && rm /tmp/miniconda.sh

# Kopiere Environment-Datei
COPY environment-linux.yml /tmp/environment.yml

# Erstelle conda-Umgebung
RUN conda env create -f /tmp/environment.yml

# Aktiviere Umgebung
SHELL ["conda", "run", "-n", "amp-rbc-md", "/bin/bash", "-c"]

# Installiere JAX 0.4.25 mit CUDA-Unterstützung (VOR ColabFold)
RUN pip install jax==0.4.25 jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Installiere PyTorch 2.3.0+cu121
RUN pip install torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Installiere spezifische cudnn-Version
RUN pip install nvidia-cudnn-cu12==8.9.2.26

# Installiere ColabFold OHNE Abhängigkeiten zu überschreiben
RUN pip install colabfold==1.5.5 --no-deps

# Installiere fehlende ColabFold-Abhängigkeiten manuell
RUN pip install absl-py appdirs biopython matplotlib numpy pandas py3Dmol requests tqdm dm-haiku

# Kopiere Projekt
COPY . /app
WORKDIR /app

# Installiere Projekt
RUN pip install -e .

# Teste Installation
RUN python -c "import jax; import jaxlib; import torch; import colabfold; print('=== DOCKER INSTALLATION ERFOLGREICH ==='); print(f'JAX Version: {jax.__version__}'); print(f'JAXlib Version: {jaxlib.__version__}'); print(f'JAX Devices: {jax.devices()}'); print(f'JAX linear_util verfügbar: {hasattr(jax, \"linear_util\")}'); print(f'PyTorch Version: {torch.__version__}'); print(f'PyTorch CUDA verfügbar: {torch.cuda.is_available()}'); print(f'ColabFold verfügbar: {hasattr(colabfold, \"batch\")}')"

# Setze Standard-Shell
SHELL ["/bin/bash", "-c"]

# Aktiviere conda-Umgebung automatisch
RUN echo "conda activate amp-rbc-md" >> ~/.bashrc

# Exponiere Ports falls benötigt
EXPOSE 8888

# Standard-Befehl
CMD ["conda", "run", "-n", "amp-rbc-md", "bash"] 