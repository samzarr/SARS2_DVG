#!/bin/bash
set -eo pipefail

# Argument parsing: supports both legacy (PRIMER REPO) and explicit (SRA PRIMER REPO THREADS CLEANUP)
if [[ "$1" == *.bed ]] || [ -f "$1" ] || [ -z "$1" ]; then
  # Legacy mode: current directory is sample directory
  SRA="$(basename "$PWD")"
  PRIMER_INPUT="$1"
  REPO="${2:-/pipeline/resources}"
  THREADS="${3:-4}"
  CLEANUP="${4:-false}"
else
  # Direct / HTCondor mode
  SRA="$1"
  PRIMER_INPUT="$2"
  REPO="${3:-/pipeline/resources}"
  THREADS="${4:-4}"
  CLEANUP="${5:-true}"
fi

# Ensure THREADS is a valid positive integer
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  THREADS=1
fi

name="${SRA}"
REF="${REPO}/NC_045512.2.fasta"
GENOME_DIR="${REPO}/Genome_Dir"

# Resolve primer BED file
case "$PRIMER_INPUT" in
  3|v3|V3)            PRIMER="${REPO}/ARTIC_primers_v3.bed" ;;
  4|v4|V4)            PRIMER="${REPO}/ARTIC_primers_v4.bed" ;;
  4.1|v4.1|4.1alt)    PRIMER="${REPO}/ARTIC_primers_v4.1.bed" ;;
  *)
    if [ -n "$PRIMER_INPUT" ] && [ -f "$PRIMER_INPUT" ]; then
      PRIMER="$PRIMER_INPUT"
    else
      PRIMER="${REPO}/ARTIC_primers_v4.bed"
    fi
    ;;
esac


# ==============================================================================
# Helper Functions: Step Timing and Worker Resource Tracking
# ==============================================================================
PIPELINE_START=$(date +%s)
STEP_START=0
CURRENT_STEP=""

get_memory_peak() {
  local mem_bytes=""
  if [ -r "/sys/fs/cgroup/memory.peak" ]; then
    mem_bytes=$(cat /sys/fs/cgroup/memory.peak 2>/dev/null || true)
  elif [ -r "/sys/fs/cgroup/memory/memory.max_usage_in_bytes" ]; then
    mem_bytes=$(cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes 2>/dev/null || true)
  fi
  if [[ "$mem_bytes" =~ ^[0-9]+$ ]]; then
    echo "$((mem_bytes / 1024 / 1024)) MB"
  else
    echo "N/A"
  fi
}

log_step() {
  local NEXT_STEP="$1"
  local NOW=$(date +%s)
  if [ "$STEP_START" -ne 0 ]; then
    local ELAPSED=$((NOW - STEP_START))
    local DISK_NOW=$(du -sh . 2>/dev/null | cut -f1)
    local RAM_NOW=$(get_memory_peak)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [TIMING] Completed '${CURRENT_STEP}' in ${ELAPSED}s ($((ELAPSED / 60))m $((ELAPSED % 60))s) | Scratch: ${DISK_NOW} | Peak RAM: ${RAM_NOW}"
  fi
  STEP_START=$NOW
  CURRENT_STEP="$NEXT_STEP"
  echo "=========================================================="
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: ${NEXT_STEP}"
  echo "=========================================================="
}

# Helper filter: drop high-volume per-read stderr prints while preserving genuine warnings and errors
filter_verbose_stderr() {
  sed -e '/^Shifted alignment read/d' \
      -e '/^No closest primer found for read/d' \
      -e '/has multiple alignments with different primer pools/d' >&2
}

echo "=========================================================="
echo "Running DVG pipeline for ${SRA}"
echo "Primer BED: ${PRIMER}"
echo "Threads:    ${THREADS}"
echo "Cleanup:    ${CLEANUP}"
echo "Initial Scratch Disk: $(du -sh . 2>/dev/null | cut -f1)"
echo "Initial Peak RAM:     $(get_memory_peak)"
echo "=========================================================="

# 1. Input acquisition and FASTQ extraction
log_step "1. Input Acquisition and FASTQ Extraction"

# Check if FASTQs already exist
if [ ! -s "${SRA}_1.fastq" ] || [ ! -s "${SRA}_2.fastq" ]; then
  # If .sra archive is not already present, download it on the compute node
  if [ ! -s "./${SRA}.sra" ]; then
    echo "Input ./${SRA}.sra not found in scratch. Downloading ${SRA} on compute node..."
    DOWNLOADED=false

    # Method A: prefetch from SRA Toolkit (installed in container)
    if command -v prefetch >/dev/null 2>&1; then
      echo "  Attempting download via prefetch..."
      if prefetch "${SRA}" -O ./ 2>/dev/null; then
        if [ -f "./${SRA}/${SRA}.sra" ]; then
          mv "./${SRA}/${SRA}.sra" "./${SRA}.sra"
          rmdir "./${SRA}" 2>/dev/null || true
        elif [ -f "./${SRA}.sra" ]; then
          :
        fi
        [ -s "./${SRA}.sra" ] && DOWNLOADED=true
      fi
    fi

    # Method B: Direct HTTPS download via curl from NCBI SRA AWS Open Data
    if [ "$DOWNLOADED" = false ]; then
      S3_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRA}/${SRA}"
      if command -v curl >/dev/null 2>&1; then
        echo "  Attempting direct HTTPS download from NCBI Open Data via curl..."
        if curl -f -L --retry 3 -sS -o "./${SRA}.sra" "$S3_URL"; then
          [ -s "./${SRA}.sra" ] && DOWNLOADED=true
        fi
      fi
    fi

    # Method C: Direct HTTPS download via wget
    if [ "$DOWNLOADED" = false ]; then
      S3_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRA}/${SRA}"
      if command -v wget >/dev/null 2>&1; then
        echo "  Attempting direct HTTPS download from NCBI Open Data via wget..."
        if wget -q -O "./${SRA}.sra" "$S3_URL"; then
          [ -s "./${SRA}.sra" ] && DOWNLOADED=true
        fi
      fi
    fi

    # Method D: Python standard library urllib (always available in container)
    if [ "$DOWNLOADED" = false ]; then
      echo "  Attempting direct HTTPS download via Python urllib..."
      if python3 -c "import urllib.request; urllib.request.urlretrieve('https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRA}/${SRA}', './${SRA}.sra')" 2>/dev/null; then
        [ -s "./${SRA}.sra" ] && DOWNLOADED=true
      fi
    fi

    if [ "$DOWNLOADED" = false ] || [ ! -s "./${SRA}.sra" ]; then
      echo "Error: Failed to download SRA archive for ${SRA} on compute node." >&2
      exit 1
    fi
    echo "Successfully downloaded ${SRA}.sra ($(du -h "./${SRA}.sra" 2>/dev/null | cut -f1)) on compute node."
  fi

  # Extract FASTQs from .sra archive
  echo "Extracting FASTQs from ./${SRA}.sra with fasterq-dump..."
  fasterq-dump --split-3 "./${SRA}.sra" -O ./ -e "${THREADS}"

  # Immediately remove .sra to free scratch disk space for alignment steps
  rm -f "./${SRA}.sra"
fi

if [ ! -s "${SRA}_1.fastq" ] || [ ! -s "${SRA}_2.fastq" ]; then
  echo "Error: Required input FASTQs (${SRA}_1.fastq, ${SRA}_2.fastq) missing or empty." >&2
  exit 1
fi

# 2. Adapter and quality trimming
log_step "2. Adapter and Quality Trimming (Trim Galore)"
echo "Trimming adapters and filtering quality with Trim Galore..."
trim_galore --cores "${THREADS}" --stringency 3 -q 30 -e .10 --length 15 --paired "${SRA}_1.fastq" "${SRA}_2.fastq"

if [ ! -s "${SRA}_1_val_1.fq" ] || [ ! -s "${SRA}_2_val_2.fq" ]; then
  echo "Error: Trim Galore failed to produce validated paired FASTQs." >&2
  exit 1
fi

# 3. STAR chimeric alignment
log_step "3. STAR Chimeric Alignment"
echo "Aligning reads with STAR..."
STAR --runThreadN "${THREADS}" \
  --readFilesIn "./${SRA}_1_val_1.fq" "./${SRA}_2_val_2.fq" \
  --outFileNamePrefix "${name}" \
  --genomeDir "${GENOME_DIR}" \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outSJfilterOverhangMin 12 12 12 12 \
  --outSJfilterCountUniqueMin 1 1 1 1 \
  --outSJfilterCountTotalMin 1 1 1 1 \
  --outSJfilterDistToOtherSJmin 0 0 0 0 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --scoreGapNoncan -4 --scoreGapATAC -4 \
  --chimScoreJunctionNonGTAG 0 \
  --chimOutType Junctions WithinBAM HardClip \
  --alignSJstitchMismatchNmax -1 -1 -1 -1 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000

# 4. Standardization, annotation, filtering, and deletion extraction
log_step "4. Standardization, Annotation, and Deletion Extraction"
echo "Standardizing alignments..."
{
  python3 "${REPO}/standardize_alignments.py" "${REF}" "${SRA}Aligned.out.sam" \
    > "${SRA}_standardizationAligned.out.bam"
} 2>&1 | filter_verbose_stderr

echo "Annotating alignments with primers..."
{
  python3 "${REPO}/annotate_alignment_with_primers.py" "${PRIMER}" \
    "${SRA}_standardizationAligned.out.bam" \
    > "${SRA}_annotated_standardization_Aligned.out.bam"
} 2>&1 | filter_verbose_stderr

echo "Filtering aligned reads..."
{
  python3 "${REPO}/filter_aligned_reads.py" "${PRIMER}" \
    "${SRA}_annotated_standardization_Aligned.out.bam" \
    --min-deletion-length 20 --max-overhang-primer-frac 1 --min-aligned-length 75 \
    --primer-pool-matching --max-primer-dist 1 \
    > "filtered_${SRA}_annotated_standardizationAligned.out.bam"
} 2>&1 | filter_verbose_stderr

echo "Extracting deletions..."
python3 "${REPO}/extract_deletions.py" \
  "filtered_${SRA}_annotated_standardizationAligned.out.bam" \
  --primer-bed "${PRIMER}" --min-deletion-length 20 --ignore-secondary \
  > "filtered_${SRA}_annotated_standardizationAligned.deletion.sorted.txt"

# 5. Sorting, depth, and R summarization
log_step "5. Sorting, Depth, and R Summarization"
echo "Sorting filtered BAM..."
samtools sort -@ "${THREADS}" "filtered_${name}_annotated_standardizationAligned.out.bam" \
  -o "filtered_${name}_annotated_standardizationAligned.out.sorted.bam"

echo "Calculating coverage..."
samtools depth -a -m 0 "filtered_${name}_annotated_standardizationAligned.out.sorted.bam" \
  > "filtered_${name}_annotated_standardizationAligned.out.sorted.coverage"

echo "Running R summarization and filtration..."
Rscript "${REPO}/deletion_summarization_and_further_filtration.R" \
  "filtered_${name}_annotated_standardizationAligned.deletion.sorted.txt" \
  "filtered_${name}_annotated_standardizationAligned.out.sorted.coverage"

# Wrap up Step 5 timing and report pre-cleanup peak resources
NOW=$(date +%s)
if [ "$STEP_START" -ne 0 ]; then
  ELAPSED=$((NOW - STEP_START))
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [TIMING] Completed '${CURRENT_STEP}' in ${ELAPSED}s ($((ELAPSED / 60))m $((ELAPSED % 60))s)"
fi

PEAK_DISK=$(du -sh . 2>/dev/null | cut -f1)
PEAK_RAM=$(get_memory_peak)
echo "=========================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WORKER RESOURCE USAGE: PRE-CLEANUP PEAK]"
echo "  Peak Scratch Storage:  ${PEAK_DISK}"
echo "  Peak Container Memory: ${PEAK_RAM}"
echo "  Largest scratch files before cleanup:"
du -ah . 2>/dev/null | sort -rh | head -n 8 | sed 's/^/    /'
echo "=========================================================="

# 6. Cleanup large intermediate files if requested
if [ "${CLEANUP}" = "true" ]; then
  STEP_START=$(date +%s)
  CURRENT_STEP="6. Intermediate File Cleanup"
  echo "=========================================================="
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: ${CURRENT_STEP}"
  echo "=========================================================="
  echo "Cleaning up large intermediate files for ${SRA}..."
  rm -f ./${SRA}.sra ./${SRA}*.fastq ./${SRA}*.fq ./${SRA}Aligned.out.sam \
        ./${SRA}_standardizationAligned.out.bam ./${SRA}_annotated_standardization_Aligned.out.bam \
        ./filtered_${SRA}_annotated_standardizationAligned.out.bam \
        ./filtered_${SRA}_annotated_standardizationAligned.count.txt \
        ./filtered_${name}_annotated_standardizationAligned.out.sorted.bam \
        ./filtered_${name}_annotated_standardizationAligned.deletion.sorted.txt \
        ./filtered_${name}_annotated_standardizationAligned.out.sorted.coverage \
        ./${name}Log.progress.out ./${name}SJ.out.tab
  NOW=$(date +%s)
  ELAPSED=$((NOW - STEP_START))
  FINAL_DISK=$(du -sh . 2>/dev/null | cut -f1)
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [TIMING] Completed '${CURRENT_STEP}' in ${ELAPSED}s | Scratch Storage after cleanup: ${FINAL_DISK}"
fi

PIPELINE_END=$(date +%s)
TOTAL_ELAPSED=$((PIPELINE_END - PIPELINE_START))
FINAL_PEAK_RAM=$(get_memory_peak)
FINAL_SCRATCH=$(du -sh . 2>/dev/null | cut -f1)

echo "=========================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline finished successfully for ${SRA}."
echo "Total Execution Duration:     ${TOTAL_ELAPSED}s ($((TOTAL_ELAPSED / 60))m $((TOTAL_ELAPSED % 60))s)"
echo "Final Scratch Disk Size:      ${FINAL_SCRATCH}"
echo "Maximum Memory Recorded:      ${FINAL_PEAK_RAM}"
echo "=========================================================="