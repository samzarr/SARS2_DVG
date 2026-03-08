# SARS-CoV-2 Deletion Detection Pipeline (Docker)

This repository contains a Dockerized workflow for running a SARS-CoV-2 deletion-detection pipeline.

The container includes all required dependencies (Python tools, STAR, samtools, and R packages) and allows the pipeline to be executed reproducibly on different systems.

This setup is intended for local testing and for future execution on CHTC.

---

## Repository Structure

```text
mycontainer/
├── Dockerfile
├── dvg_env_clean.yml
├── resources/
│   ├── Genome_Dir/
│   ├── NC_045512.2.fasta
│   ├── ARTIC_primers_v3.bed
│   ├── ARTIC_primers_v4.bed
│   ├── ARTIC_primers_v4.1.bed
│   ├── standardize_alignments.py
│   ├── annotate_alignment_with_primers.py
│   ├── filter_aligned_reads.py
│   ├── extract_deletions.py
│   ├── deletion_summarization_and_further_filtration.R
│   └── other resource files
└── scripts/
    ├── master_automation.sh
    ├── run_pipeline.sh
    └── coguk_illumina_paired_artic.tsv
```

## Prerequisites

Before running the pipeline, install the following:

- Docker Desktop  
- Git  
- Terminal access  

Verify Docker installation:

```bash
docker --version
docker info
```
---

## Step 1: Enter the Project Directory

Open a terminal and navigate to the container directory:

```bash
cd path/to/mycontainer
```

Check the directory contents:

```bash
ls
```

You should see at least:

- `Dockerfile`
- `dvg_env_clean.yml`
- `resources/`
- `scripts/`

## Step 2: Build the Docker Image

Build the container image:

```bash
docker build -t dvg_pipeline .
```

Sometimes the conda environment installation may fail due to network issues.

If that happens, simply rerun the command:

```bash
docker build -t dvg_pipeline .
```

Verify that the image was created:

```bash
docker images | grep dvg_pipeline
```

## Step 3: Start an Interactive Container

Run the container and mount the local directory:

```bash
docker run --rm -it -v "$(pwd)":/pipeline dvg_pipeline bash
```

This mount is important.

Without the mount, output files will only exist inside the container and disappear after exit.

With the mount, files created inside `/pipeline` will appear in the local directory.

After starting the container, the prompt should look similar to:

```bash
(dvg_env) root@<container_id>:/pipeline#
```
## Step 4: Verify Software Inside the Container

Check that the required tools are available:

```bash
which prefetch fasterq-dump trim_galore STAR samtools python3 Rscript
```

Optional version checks:

```bash
STAR --version
samtools --version | head
python3 --version
Rscript --version
prefetch -V
fasterq-dump -V
```
## Step 5: Running the Pipeline

Two scripts are provided:

`master_automation.sh`  
Runs the pipeline for multiple samples listed in a TSV file.

`run_pipeline.sh`  
Helper script for single sample running.

### Run Using the Automation Script

Inside the container:

```bash
cd /pipeline/scripts
bash master_automation.sh
```

This reads the sample list from:

```bash
/pipeline/scripts/coguk_illumina_paired_artic.tsv
```

The script will automatically iterate through the samples listed in the TSV file and run the pipeline for each sample.
## Step 6: Expected Outputs

Outputs are generated inside the sample working directory.

Example:

```bash
/pipeline/scripts/ERR7512327/
```

Typical outputs include:

```text
ERR7512327Aligned.out.sam
ERR7512327_standardizationAligned.out.bam
ERR7512327_annotated_standardization_Aligned.out.bam
filtered_ERR7512327_annotated_standardizationAligned.out.bam
filtered_ERR7512327_annotated_standardizationAligned.deletion.sorted.txt
filtered_ERR7512327_annotated_standardizationAligned.out.sorted.bam
filtered_ERR7512327_annotated_standardizationAligned.out.sorted.coverage
```

The most important result file is:

```text
filtered_ERR7512327_0.01_5_20_2_34.txt
```

A log file should also be kept to record warnings or errors.

---

# Appendix

## A. Rebuilding the Container

If any of the following files change:

- `Dockerfile`
- environment YAML
- anything under `resources/`
- anything under `scripts/`

rebuild the container:

```bash
docker build -t dvg_pipeline .
```

Then rerun:

```bash
docker run --rm -it -v "$(pwd)":/pipeline dvg_pipeline bash
```


## B. Exiting the Container

Exit the container with:

```bash
exit
```

or press:

```
Ctrl + D
```


## C. Troubleshooting

### Docker command not found

Docker Desktop is not installed.

Install Docker Desktop and restart the terminal.


### Cannot connect to Docker daemon

Docker Desktop is installed but not running.

Open Docker Desktop and wait for it to start.

Check:

```bash
docker info
```


### Image not found

The Docker image was not built on the current system.

Run:

```bash
docker build -t dvg_pipeline .
```

### Conda environment build failure

Occasionally the conda environment installation fails due to network interruptions.

Retry the build command.


### Rscript command not found

The environment does not include `r-base`.

Add R to the environment file and rebuild the container.


## D. Current Status

At the time of writing:

- Docker image builds successfully
- Pipeline runs inside the container
- Test runs produce expected output files
- The next step is to run the workflow on CHTC

Suggested workflow for scaling:

```
1 sample → 10 samples → 100 samples → larger batches
```

This allows verification of runtime behavior and output consistency.
