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

tail -n +2 "$TSV" | while IFS=$'\t' read -r RUN_ID PRIMER_VER PROTOCOL_VER; do
  [ -z "$RUN_ID" ] && continue

  TARGET_SRA="${DEST_DIR}/${RUN_ID}.sra"
  if [ -s "$TARGET_SRA" ]; then
    echo ">>> ${RUN_ID}.sra already exists in ${DEST_DIR}, skipping."
    continue
  fi

  echo ">>> Downloading ${RUN_ID} via prefetch..."
  prefetch "$RUN_ID" --output-directory "$DEST_DIR" || {
    # Fallback if prefetch creates a directory <RUN_ID>/<RUN_ID>.sra
    if [ -f "${RUN_ID}/${RUN_ID}.sra" ]; then
      mv "${RUN_ID}/${RUN_ID}.sra" "$TARGET_SRA"
      rmdir "$RUN_ID" 2>/dev/null || true
    fi
  }

  # If prefetch placed it in a nested directory inside DEST_DIR
  if [ -f "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" ]; then
    mv "${DEST_DIR}/${RUN_ID}/${RUN_ID}.sra" "$TARGET_SRA"
    rmdir "${DEST_DIR}/${RUN_ID}" 2>/dev/null || true
  fi

  if [ -s "$TARGET_SRA" ]; then
    echo ">>> Successfully staged ${TARGET_SRA} ($(du -h "$TARGET_SRA" | cut -f1))"
  else
    echo "Warning: Failed to stage SRA file for ${RUN_ID}" >&2
  fi
done

echo "Staging complete."
