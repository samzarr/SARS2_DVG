#!/bin/bash
set -eo pipefail

TSV="${1:-coguk_illumina_paired_artic.tsv}"
DEST_DIR="${2:-staged_sra}"

if [ ! -f "$TSV" ]; then
  echo "Error: Sample TSV file '${TSV}' not found." >&2
  echo "Usage: $0 [sample_table.tsv] [destination_directory]" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"
echo "Staging SRA files listed in ${TSV} into ${DEST_DIR}..."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

grep -v '^[[:space:]]*#' "$TSV" | while IFS=$'\t' read -r RUN_ID PRIMER_VER PROTOCOL_VER; do
  [ -z "$RUN_ID" ] && continue
  [[ "$RUN_ID" =~ ^[rR]un_[iI]d ]] && continue

  TARGET_SRA="${DEST_DIR}/${RUN_ID}.sra"
  if [ -s "$TARGET_SRA" ]; then
    echo ">>> ${RUN_ID}.sra already exists in ${DEST_DIR}, skipping."
    continue
  fi

  # 1. Check local repository directories for pre-existing files
  LOCAL_FOUND=false
  for cand in \
    "${SCRIPT_DIR}/../../automation/${RUN_ID}/${RUN_ID}.sra" \
    "${SCRIPT_DIR}/../automation/${RUN_ID}/${RUN_ID}.sra" \
    "automation/${RUN_ID}/${RUN_ID}.sra" \
    "${RUN_ID}/${RUN_ID}.sra" \
    "${RUN_ID}.sra"; do
    if [ -s "$cand" ]; then
      echo ">>> Found existing local file: ${cand}"
      cp "$cand" "$TARGET_SRA"
      LOCAL_FOUND=true
      break
    fi
  done

  if [ "$LOCAL_FOUND" = true ]; then
    echo ">>> Successfully staged ${TARGET_SRA} ($(du -h "$TARGET_SRA" | cut -f1))"
    continue
  fi

  # 2. Download from NCBI SRA
  echo ">>> Downloading ${RUN_ID}..."
  DOWNLOADED=false

  # Method A: prefetch (if available in PATH)
  if command -v prefetch >/dev/null 2>&1; then
    echo "    Attempting download with prefetch..."
    if prefetch "$RUN_ID" --output-directory "$DEST_DIR" 2>/dev/null; then
      if [ -f "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" ]; then
        mv "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" "$TARGET_SRA"
        rmdir "${DEST_DIR}/${RUN_ID}" 2>/dev/null || true
      elif [ -f "${RUN_ID}/${RUN_ID}.sra" ]; then
        mv "${RUN_ID}/${RUN_ID}.sra" "$TARGET_SRA"
        rmdir "$RUN_ID" 2>/dev/null || true
      fi
      [ -s "$TARGET_SRA" ] && DOWNLOADED=true
    fi
  fi

  # Method B: Direct HTTPS download via curl/wget from NCBI AWS Open Data bucket
  if [ "$DOWNLOADED" = false ]; then
    S3_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${RUN_ID}/${RUN_ID}"
    if command -v curl >/dev/null 2>&1; then
      echo "    Attempting direct HTTPS download from NCBI SRA Open Data via curl..."
      if curl -f -L --retry 3 -sS -o "$TARGET_SRA" "$S3_URL"; then
        [ -s "$TARGET_SRA" ] && DOWNLOADED=true
      fi
    elif command -v wget >/dev/null 2>&1; then
      echo "    Attempting direct HTTPS download from NCBI SRA Open Data via wget..."
      if wget -q -O "$TARGET_SRA" "$S3_URL"; then
        [ -s "$TARGET_SRA" ] && DOWNLOADED=true
      fi
    fi
  fi

  # Method C: Apptainer/Singularity (standard on CHTC submit nodes)
  # if [ "$DOWNLOADED" = false ] && { command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1; }; then
  #   CONTAINER_RUNNER="apptainer"
  #   command -v apptainer >/dev/null 2>&1 || CONTAINER_RUNNER="singularity"
  #   CONTAINER_IMG="${CONTAINER_IMAGE:-docker://quay.io/biocontainers/sra-tools:3.1.1--h43eeafb_0}"
  #   echo "    Attempting download with ${CONTAINER_RUNNER} container..."
  #   if $CONTAINER_RUNNER exec "$CONTAINER_IMG" prefetch "$RUN_ID" -O "$DEST_DIR" 2>/dev/null; then
  #     if [ -f "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" ]; then
  #       mv "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" "$TARGET_SRA"
  #       rmdir "${DEST_DIR}/${RUN_ID}" 2>/dev/null || true
  #     fi
  #     [ -s "$TARGET_SRA" ] && DOWNLOADED=true
  #   fi
  # fi

  if [ -s "$TARGET_SRA" ]; then
    echo ">>> Successfully staged ${TARGET_SRA} ($(du -h "$TARGET_SRA" | cut -f1))"
  else
    echo "Warning: Failed to stage SRA file for ${RUN_ID}" >&2
    rm -f "$TARGET_SRA"
  fi
done

echo "Staging complete."
