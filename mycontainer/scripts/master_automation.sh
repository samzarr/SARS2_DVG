#!/bin/bash
set -e

# Base directory (inside Docker)
ROOT_DIR="/pipeline"

# Scripts
SCRIPTS="${ROOT_DIR}/scripts"

# SARS-CoV-2 pipeline repo (renamed as resources)
REPO="${ROOT_DIR}/resources"

# Sample sheet
TSV="${SCRIPTS}/coguk_illumina_paired_artic.tsv"

# Pipeline script
PIPELINE="${SCRIPTS}/run_pipeline.sh"

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
  # bash "$PIPELINE" "$PRIMER_FILE" "$REPO"
  LOG_FILE="${RUN_ID}_pipeline.log"
  echo "Logging output to ${LOG_FILE}"
  bash "$PIPELINE" "$PRIMER_FILE" "$REPO" 2>&1 | tee "$LOG_FILE"      

  # 5. Return to upper working directory
  cd ..
done