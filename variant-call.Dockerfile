FROM ubuntu:22.04

SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
 && apt-get install -y wget bzip2 ca-certificates curl git \
 && apt-get clean

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
     -O /tmp/miniconda.sh \
 && bash /tmp/miniconda.sh -b -p $CONDA_DIR \
 && rm /tmp/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH

# Install Mamba
RUN conda install -y -c conda-forge mamba \
 && conda clean -afy

# Install Sniffles and plotting dependencies via Conda/Mamba
RUN mamba install -y \
      -c bioconda \
      -c conda-forge \
        python=3.12 \
        sniffles=2.6.2 \
        numpy \
        pandas \
        matplotlib \
        seaborn \
        upsetplot \
    && mamba clean --all --yes

RUN pip install --no-deps sniffles2-plot

# Smoke-test both tools
RUN sniffles --version \
 && python3 -m sniffles2_plot --help

CMD ["bash"]
