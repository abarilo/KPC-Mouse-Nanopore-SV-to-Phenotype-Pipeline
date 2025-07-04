FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
      wget bzip2 ca-certificates curl git \
      samtools \
 && apt-get clean

#  Miniconda install
ENV CONDA_DIR=/opt/conda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O /tmp/miniconda.sh \
 && bash /tmp/miniconda.sh -b -p $CONDA_DIR \
 && rm /tmp/miniconda.sh \
 && $CONDA_DIR/bin/conda clean -afy

ENV PATH=$CONDA_DIR/bin:$PATH

#  Install Mamba
RUN conda install -y -c conda-forge mamba \
 && conda clean -afy

#  Install QC & mapping tools
RUN mamba install -y \
      -c bioconda \
      -c conda-forge \
      python=3.10 \
      minimap2=2.30 \
      nanoplot=1.44.1 \
    && mamba clean --all --yes

#  Verify installation
RUN NanoPlot --version \
 && minimap2 --version \
 && samtools --version

CMD ["bash"]
