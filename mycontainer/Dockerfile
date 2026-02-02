FROM continuumio/miniconda3:latest

# Copy conda env
COPY dvg_env_clean.yml /tmp/env.yml

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda env create -f /tmp/env.yml

RUN echo "source activate dvg_env" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Workdir
WORKDIR /pipeline

# Copy automation scripts
COPY scripts /pipeline/scripts

# Copy SARS-CoV-2 pipeline resources
COPY resources /pipeline/resources