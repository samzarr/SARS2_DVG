#!/bin/bash
set -eo pipefail

RUN_ID="$1"
PRIMER_VER="$2"
CPUS="${3:-4}"
REPO="${PIPELINE_REPO:-/pipeline/resources}"
if [ -f "./run_pipeline.sh" ]; then
  PIPELINE="./run_pipeline.sh"
else
  PIPELINE="${PIPELINE_SCRIPT:-/pipeline/scripts/run_pipeline.sh}"
fi


# Ensure CPUS is a valid positive integer
if ! [[ "$CPUS" =~ ^[0-9]+$ ]]; then
  CPUS=1
fi

if [ -z "$RUN_ID" ]; then
  echo "Error: RUN_ID argument missing." >&2
  echo "Usage: $0 <RUN_ID> <PRIMER_VERSION> [CPUS]" >&2
  exit 1
fi

# Ensure execution is in HTCondor scratch directory where files were transferred
if [ -n "$_CONDOR_SCRATCH_DIR" ] && [ -d "$_CONDOR_SCRATCH_DIR" ]; then
  cd "$_CONDOR_SCRATCH_DIR"
fi

LOG_FILE="${RUN_ID}_pipeline.log"

{
  echo "=========================================================="
  echo "Starting HTCondor job on host: $(hostname)"
  echo "Current directory:             $(pwd)"
  echo "Sample ID:                     ${RUN_ID}"
  echo "Primer Version:                ${PRIMER_VER}"
  echo "CPUs allocated:                ${CPUS}"
  echo "=========================================================="

  # Run pipeline for single sample with intermediate cleanup enabled
  bash "$PIPELINE" "$RUN_ID" "$PRIMER_VER" "$REPO" "$CPUS" "true"

  echo "=========================================================="
  echo "HTCondor job finished successfully for ${RUN_ID}."
  echo "Remaining output files in scratch:"
  ls -lh
  echo "=========================================================="
} 2>&1 | tee "$LOG_FILE"
