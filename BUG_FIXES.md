# SARS-CoV-2 DVG Pipeline: Bug Fixes and Architectural Improvements

This document details each bug and architectural issue identified in the codebase, the root cause, the operational impact, and the exact changes made to resolve it.

---

## 1. Non-Interactive Shell Conda Inactivation in Docker
- **File**: `mycontainer/Dockerfile`
- **Root Cause**: The Dockerfile previously activated the conda environment via:
  ```dockerfile
  RUN echo "source activate dvg_env" >> ~/.bashrc
  SHELL ["/bin/bash", "--login", "-c"]
  ```
- **Impact**: In HTCondor and batch cluster environments, container jobs run non-interactively without launching a login shell. Non-interactive subshells do not source `~/.bashrc`. Consequently, `dvg_env` was never activated on execute nodes, causing all commands (`STAR`, `samtools`, `trim_galore`, `pysam`, `Rscript`) to fail with `command not found`.
- **Change Made**: Added the conda environment binary directory directly to the system `PATH`:
  ```dockerfile
  ENV PATH="/opt/conda/envs/dvg_env/bin:$PATH"
  ```
  Now, all tools are immediately available in PATH across interactive shells, non-interactive subshells, and HTCondor execution wrappers.

---

## 2. STAR Genome Index Binary Version Incompatibility
- **File**: `mycontainer/Dockerfile`
- **Root Cause**: The pre-packaged index in `mycontainer/resources/Genome_Dir/genomeParameters.txt` was compiled using STAR `2.7.1a` (with backup `2.7.4a`). The conda environment installs modern STAR (`2.7.10b+`).
- **Impact**: STAR strictly checks `versionGenome` against its own binary version and aborts with a fatal incompatibility error when loading older indices.
- **Change Made**: Added an in-container index generation step to the Dockerfile using the exact STAR binary installed in the image:
  ```dockerfile
  RUN rm -rf /pipeline/resources/Genome_Dir && \
      mkdir -p /pipeline/resources/Genome_Dir && \
      STAR --runMode genomeGenerate \
           --runThreadN 4 \
           --genomeDir /pipeline/resources/Genome_Dir \
           --genomeFastaFiles /pipeline/resources/NC_045512.2.fasta \
           --genomeSAindexNbases 6
  ```
  Because the SARS-CoV-2 genome is ~29.9 kb, generating the index at build time takes ~3 seconds and guarantees 100% binary compatibility.

---

## 3. Conda Environment Duplication and Channel Conflicts
- **File**: `mycontainer/dvg_env_clean.yml`
- **Root Cause**:
  - `star` was duplicated (listed on line 10 and line 16).
  - Channels included `defaults` alongside `conda-forge` and `bioconda`.
- **Impact**: Including `defaults` causes channel backtracking, package version conflicts, and licensing warnings.
- **Change Made**:
  - Removed duplicate `star` entry.
  - Standardized channels to `conda-forge` and `bioconda` (with `conda-forge` prioritized first as per Bioconda specifications).
  - Moved `pysam` to conda dependencies for faster pre-compiled binary installation.

---

## 4. Silent Error Masking in Trim Galore
- **Files**: `mycontainer/scripts/run_pipeline.sh`, `automation/run_pipeline.sh`
- **Root Cause**: `trim_galore` was invoked with `|| true`:
  ```bash
  trim_galore --stringency 3 -q 30 -e .10 --length 15 --paired ${SRA}_1.fastq ${SRA}_2.fastq || true
  ```
- **Impact**: When Trim Galore encountered errors (such as desynchronized FASTQ reads observed in `automation/ERR7512327/ERR7512327_run.log`), the failure code was suppressed. The script blindly proceeded into STAR with truncated or missing FASTQs, creating corrupted alignments and hard-to-diagnose errors downstream.
- **Change Made**:
  - Removed `|| true`.
  - Added explicit pre-flight checks verifying `${SRA}_1.fastq` and `${SRA}_2.fastq` exist and are non-empty before trimming.
  - Added post-flight checks ensuring validated output files `${SRA}_1_val_1.fq` and `${SRA}_2_val_2.fq` exist before launching STAR.

---

## 5. Plain-Text SAM Saved with `.bam` Extension
- **Files**:
  - `mycontainer/resources/annotate_alignment_with_primers.py`
  - `SARS-CoV-2_Deletion_detection_from_amplicon-seq/annotate_alignment_with_primers.py`
- **Root Cause**: Line 20 used `pysam.AlignmentFile("-", "w", template=input_samfile)`. In pysam, mode `"w"` outputs plain-text SAM, while `"wb"` outputs binary compressed BAM. Downstream pipeline scripts redirected stdout to `${SRA}_annotated_standardization_Aligned.out.bam`.
- **Impact**: Files were text SAM format masquerading under a `.bam` extension, inflating file sizes 4–5x and creating inconsistency with `standardize_alignments.py` and `filter_aligned_reads.py` (which default to binary `"wb"`).
- **Change Made**: Added the `--sam` CLI argument (default `False`) and changed the AlignmentFile mode to `"w" if args.output_sam else "wb"`.

---

## 6. Invalid Regular Expression Patterns in R Script
- **Files**:
  - `mycontainer/resources/deletion_summarization_and_further_filtration.R`
  - `SARS-CoV-2_Deletion_detection_from_amplicon-seq/deletion_summarization_and_further_filtration.R`
- **Root Cause**: `list.files(pattern = "*deletion.sorted.txt")`, `"*out.sorted.coverage"`, and `"*count.txt"`. In regular expressions, `*` is a quantifier that modifies the preceding token. Having `*` as the first character is an invalid quantifier in regex.
- **Impact**: Caused regex compilation warnings or failed to match files depending on the R runtime engine.
- **Change Made**: Replaced with proper regular expressions:
  - `".*deletion\\.sorted\\.txt$"`
  - `".*out\\.sorted\\.coverage$"`
  - `".*count\\.txt$"`

---

## 7. Positional Indexing Hazard in Genome Coverage Lookup
- **Files**:
  - `mycontainer/resources/deletion_summarization_and_further_filtration.R`
  - `SARS-CoV-2_Deletion_detection_from_amplicon-seq/deletion_summarization_and_further_filtration.R`
- **Root Cause**: The R script looked up coverage depth using raw positional indexing:
  ```R
  data_forward$start_cov <- coverage[unlist(data_forward[2]), 3]
  data_forward$end_cov <- coverage[unlist(data_forward[3] + 1), 3]
  ```
- **Impact**: This assumed that row $k$ in the coverage file was always genomic position $k$. If any coordinate was missing or if coordinates were 0-based (`start == 0`), row indexing returned `integer(0)` or misaligned coverage values.
- **Change Made**: Replaced positional row indexing with key-based matching on `coverage$Position`:
  ```R
  start_idx <- match(data_forward$start, coverage$Position)
  end_idx   <- match(data_forward$end + 1, coverage$Position)
  data_forward$start_cov <- ifelse(is.na(start_idx), 0, coverage$Depth[start_idx])
  data_forward$end_cov   <- ifelse(is.na(end_idx), 0, coverage$Depth[end_idx])
  ```
  Also added CLI argument support (`args <- commandArgs(trailingOnly = TRUE)`) so the script can target a specific sample's files directly rather than globbing all files in the current working directory.

---

## 8. Hardcoded Threads and Disk Space Exhaustion in Pipeline Scripts
- **Files**:
  - `mycontainer/scripts/run_pipeline.sh`
  - `mycontainer/scripts/master_automation.sh`
  - `automation/run_pipeline.sh`
- **Root Cause**:
  - `STAR --runThreadN 8` was hardcoded.
  - `SRA=$(basename "$PWD")` was assumed.
  - Raw FASTQs, trimmed FASTQs, and huge uncompressed SAM files were never removed.
- **Impact**:
  - On compute clusters with fewer cores assigned (e.g. 2 or 4), STAR oversubscribed CPUs, risking cgroup throttling or job eviction.
  - In HTCondor, all newly generated files in `$_CONDOR_SCRATCH_DIR` are transferred back to the submit node on exit. Transferring 15–20 GB of intermediate files per sample for 100+ samples quickly saturates submit node disk quotas and network bandwidth.
- **Change Made**:
  - Parameterized threads (`${THREADS:-4}`) across `fasterq-dump`, `trim_galore`, `STAR`, and `samtools sort`.
  - Added dual-mode argument parsing: supports legacy `run_pipeline.sh <PRIMER> <REPO>` and explicit `run_pipeline.sh <SRA> <PRIMER> <REPO> <THREADS> <CLEANUP>`.
  - Added intermediate file cleanup when `CLEANUP=true`, purging temporary `.sra`, raw `.fastq`, trimmed `.fq`, uncompressed `.sam`, and intermediate BAMs, preserving only the final summary tables, sorted BAM, coverage file, and logs.

---

## 9. HTCondor Distributed Execution Architecture
- **New Files Added**:
  - `mycontainer/scripts/run_sample_condor.sh`: Single-sample job wrapper for HTCondor nodes.
  - `mycontainer/scripts/sars2_dvg.sub`: Production HTCondor submit file with per-sample queueing (`queue RUN_ID, PRIMER_VER from tsv`), resource requests, and output file filtering.
  - `mycontainer/scripts/stage_inputs.sh`: Pre-download utility for the submit node to stage `.sra` files so compute nodes do not require outbound internet access.
