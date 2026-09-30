FROM continuumio/miniconda3:26.7.1-1

# Ensure conda environment binaries are directly on PATH for non-interactive HTCondor execution
ENV PATH="/opt/conda/envs/dvg_env/bin:$PATH"

# Copy conda env
COPY dvg_env_clean.yml /tmp/env.yml

RUN conda env create -f /tmp/env.yml && \
    conda clean -afy

# Workdir
WORKDIR /pipeline

# Copy SARS-CoV-2 pipeline resources
COPY resources /pipeline/resources

# Copy automation scripts
COPY scripts /pipeline/scripts

# Re-generate STAR genome index with installed STAR binary to guarantee version compatibility
RUN rm -rf /pipeline/resources/Genome_Dir && \
    mkdir -p /pipeline/resources/Genome_Dir && \
    STAR --runMode genomeGenerate \
         --runThreadN 4 \
         --genomeDir /pipeline/resources/Genome_Dir \
         --genomeFastaFiles /pipeline/resources/NC_045512.2.fasta \
         --genomeSAindexNbases 6

