#!/bin/bash
set -e

# Input table of sample list
TSV="coguk_illumina_paired_artic.tsv"

PIPELINE="$(pwd)/run_pipeline.sh"
REPO="$(realpath ../SARS-CoV-2_Deletion_detection_from_amplicon-seq)"

# Primer versions
PRIMER_V3="${REPO}/ARTIC_primers_v3.bed"
PRIMER_V4="${REPO}/ARTIC_primers_v4.bed"
PRIMER_V41="${REPO}/ARTIC_primers_v4.1.bed"

# Read in TSV (skip header)
tail -n +2 "$TSV" | while IFS=$'\t' read -r RUN_ID PRIMER_VER PROTOCOL_VER; do
  echo "==============================="
  echo ">>> Processing ${RUN_ID} (primer: ${PRIMER_VER})"
  echo "==============================="

  # 1. Download raw .SRA file
  prefetch "$RUN_ID"

  # 2. Enter sample directory
  cd "$RUN_ID"

  # 3. Dynamically select primer version
  case "$PRIMER_VER" in
    3|v3|V3)
      PRIMER_FILE="$PRIMER_V3"
      ;;
    4|v4|V4)
      PRIMER_FILE="$PRIMER_V4"
      ;;
    4.1|v4.1|4.1alt|v4.1alt)
      PRIMER_FILE="$PRIMER_V41"
      ;;
    *)
      echo "Warning: Unknown primer version '${PRIMER_VER}' for ${RUN_ID} — skipping."
      cd ..
      continue
      ;;
  esac

  # 4. Call pipeline to finish all steps
  bash "$PIPELINE" "$PRIMER_FILE" "$REPO"

  # 5. Return to upper working directory
  cd ..
done