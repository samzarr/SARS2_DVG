#!/bin/bash
set -eo pipefail

THREADS="${1:-4}"
CLEANUP="${2:-true}"

# Base directory (inside Docker or local workspace)
ROOT_DIR="${PIPELINE_ROOT:-/pipeline}"
SCRIPTS="${ROOT_DIR}/scripts"
REPO="${ROOT_DIR}/resources"
TSV="${SCRIPTS}/coguk_illumina_paired_artic.tsv"
PIPELINE="${SCRIPTS}/run_pipeline.sh"

PRIMER_V3="${REPO}/ARTIC_primers_v3.bed"
PRIMER_V4="${REPO}/ARTIC_primers_v4.bed"
PRIMER_V41="${REPO}/ARTIC_primers_v4.1.bed"

if [ ! -f "$TSV" ]; then
  echo "Error: TSV file not found at ${TSV}" >&2
  exit 1
fi

# Process samples sequentially
grep -v '^[[:space:]]*#' "$TSV" | while IFS=$'\t' read -r RUN_ID PRIMER_VER PROTOCOL_VER; do
  [ -z "$RUN_ID" ] && continue
  [[ "$RUN_ID" =~ ^[rR]un_[iI]d ]] && continue

  echo "=========================================================="
  echo ">>> Processing ${RUN_ID} (Primer: ${PRIMER_VER}, Threads: ${THREADS})"
  echo "=========================================================="

  # 1. Download raw .SRA file if not already present
  if [ ! -d "$RUN_ID" ] || [ ! -f "${RUN_ID}/${RUN_ID}.sra" ]; then
    prefetch "$RUN_ID"
  fi

  # 2. Enter sample directory
  cd "$RUN_ID"

  # 3. Dynamically select primer version
  case "$PRIMER_VER" in
    3|v3|V3)                 PRIMER_FILE="$PRIMER_V3" ;;
    4|v4|V4)                 PRIMER_FILE="$PRIMER_V4" ;;
    4.1|v4.1|4.1alt|v4.1alt) PRIMER_FILE="$PRIMER_V41" ;;
    *)
      echo "Warning: Unknown primer version '${PRIMER_VER}' for ${RUN_ID} — skipping." >&2
      cd ..
      continue
      ;;
  esac

  # 4. Call pipeline with logging
  LOG_FILE="${RUN_ID}_pipeline.log"
  echo "Logging output to ${LOG_FILE}"
  bash "$PIPELINE" "$PRIMER_FILE" "$REPO" "$THREADS" "$CLEANUP" 2>&1 | tee "$LOG_FILE"

  # 5. Return to upper working directory
  cd ..
done