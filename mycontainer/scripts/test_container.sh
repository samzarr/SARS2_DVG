#!/bin/bash
set -eo pipefail

# ==============================================================================
# SARS-CoV-2 DVG Pipeline: Container Verification & Test Script
# ==============================================================================
# This script tests your rebuilt Docker image to ensure it is healthy and ready
# for HTCondor cluster execution.
#
# Tests performed:
#   1. Container Environment Check (PATH, bioinformatic binaries, Python & R libs)
#   2. Pipeline Resource Check (STAR genome index, ARTIC primer beds, scripts)
#   3. Full End-to-End Test (Runs a real sample in an isolated scratch folder,
#      simulating the exact HTCondor execution model)
#
# Usage:
#   ./test_container.sh [IMAGE_NAME] [--quick | --full]
# Examples:
#   ./test_container.sh                      # Test default 'dvg_pipeline' image
#   ./test_container.sh myuser/sars2_dvg:v1  # Test custom image tag
#   ./test_container.sh dvg_pipeline --quick # Run environment check only
# ==============================================================================

IMAGE="dvg_pipeline"
MODE="full"

for arg in "$@"; do
  if [[ "$arg" == "--quick" ]]; then
    MODE="quick"
  elif [[ "$arg" == "--full" ]]; then
    MODE="full"
  else
    IMAGE="$arg"
  fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

pass() { echo -e "  [${GREEN}PASS${NC}] $1"; }
fail() { echo -e "  [${RED}FAIL${NC}] $1"; exit 1; }
info() { echo -e "${BLUE}=== $1 ===${NC}"; }
warn() { echo -e "${YELLOW}>>> WARNING: $1${NC}"; }

echo -e "${BOLD}===================================================================${NC}"
echo -e "${BOLD}SARS-CoV-2 DVG Pipeline: Docker Verification Suite${NC}"
echo -e "Target Docker Image: ${BLUE}${IMAGE}${NC}"
echo -e "Test Mode:           ${BLUE}${MODE}${NC}"
echo -e "${BOLD}===================================================================${NC}\n"

# ------------------------------------------------------------------------------
# 0. Pre-flight Checks (Docker host)
# ------------------------------------------------------------------------------
info "Checking Docker host availability"
if ! command -v docker >/dev/null 2>&1; then
  fail "docker command not found. Ensure Docker Desktop / Engine is installed and in PATH."
fi

if ! docker info >/dev/null 2>&1; then
  fail "Cannot connect to Docker daemon. Ensure Docker is running."
fi
pass "Docker daemon is active"

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  fail "Image '${IMAGE}' not found locally. Run: docker build -t ${IMAGE} -f mycontainer/Dockerfile mycontainer/"
fi
pass "Image '${IMAGE}' found locally"

# ------------------------------------------------------------------------------
# 1. Environment & Software Sanity Check (Inside Container)
# ------------------------------------------------------------------------------
echo ""
info "Phase 1: Verifying software tools, libraries & PATH inside container"

docker run --rm "$IMAGE" bash -c '
set -e
echo "Testing PATH: $PATH"

# 1. Check binaries
for cmd in prefetch fasterq-dump trim_galore cutadapt STAR samtools python3 Rscript; do
  if command -v "$cmd" >/dev/null 2>&1; then
    printf "  [PASS] Binary: %-15s (%s)\n" "$cmd" "$(command -v "$cmd")"
  else
    echo "  [FAIL] Missing binary on PATH: $cmd" >&2
    exit 1
  fi
done

# 2. Check Python packages
python3 -c "
import pysam, pandas, numpy, scipy
print(f\"  [PASS] Python libraries OK (pysam {pysam.__version__}, pandas {pandas.__version__})\")
"

# 3. Check R packages
Rscript -e "
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})
cat(\"  [PASS] R libraries loaded successfully (dplyr, readr, tidyr, ggplot2, tibble)\n\")
"

# 4. Check resources & scripts
[ -f "/pipeline/resources/NC_045512.2.fasta" ] || { echo "FAIL: NC_045512.2.fasta missing"; exit 1; }
[ -f "/pipeline/resources/ARTIC_primers_v3.bed" ] || { echo "FAIL: ARTIC_primers_v3.bed missing"; exit 1; }
[ -f "/pipeline/resources/ARTIC_primers_v4.bed" ] || { echo "FAIL: ARTIC_primers_v4.bed missing"; exit 1; }
[ -f "/pipeline/resources/ARTIC_primers_v4.1.bed" ] || { echo "FAIL: ARTIC_primers_v4.1.bed missing"; exit 1; }
[ -f "/pipeline/scripts/run_pipeline.sh" ] || { echo "FAIL: run_pipeline.sh missing"; exit 1; }
[ -f "/pipeline/scripts/run_sample_condor.sh" ] || { echo "FAIL: run_sample_condor.sh missing"; exit 1; }
echo "  [PASS] Pipeline scripts and primer BED files verified in /pipeline"

# 5. Check STAR Genome Index
if [ -d "/pipeline/resources/Genome_Dir" ] && [ -f "/pipeline/resources/Genome_Dir/genomeParameters.txt" ]; then
  INDEX_VER=$(grep "versionGenome" /pipeline/resources/Genome_Dir/genomeParameters.txt | cut -f2 || true)
  echo "  [PASS] STAR Genome_Dir found (Compiled for STAR version: ${INDEX_VER:-unknown})"
else
  echo "  [FAIL] STAR Genome_Dir missing or incomplete in /pipeline/resources/Genome_Dir" >&2
  exit 1
fi
'
pass "Container environment sanity check passed completely!"

if [[ "$MODE" == "quick" ]]; then
  echo ""
  echo -e "${GREEN}${BOLD}Quick test passed successfully!${NC}"
  echo "The container environment is properly configured."
  exit 0
fi

# ------------------------------------------------------------------------------
# 2. End-to-End Pipeline Execution Test (Simulating HTCondor)
# ------------------------------------------------------------------------------
echo ""
info "Phase 2: End-to-end execution test on real sample (ERR7512327)"

TEST_SRA=""
CANDIDATE_PATHS=(
  "${SCRIPT_DIR}/automation/ERR7512327/ERR7512327.sra"
  "${SCRIPT_DIR}/../automation/ERR7512327/ERR7512327.sra"
  "${SCRIPT_DIR}/../../automation/ERR7512327/ERR7512327.sra"
  "$(pwd)/automation/ERR7512327/ERR7512327.sra"
  "$(pwd)/../automation/ERR7512327/ERR7512327.sra"
  "${SCRIPT_DIR}/ERR7512327.sra"
  "$(pwd)/ERR7512327.sra"
  "${SCRIPT_DIR}/staged_sra/ERR7512327.sra"
  "$(pwd)/staged_sra/ERR7512327.sra"
)

for p in "${CANDIDATE_PATHS[@]}"; do
  if [ -s "$p" ]; then
    TEST_SRA="$p"
    break
  fi
done

TEST_SCRATCH="${SCRIPT_DIR}/test_run_tmp_$$"
mkdir -p "$TEST_SCRATCH"

cleanup_test_scratch() {
  if [ -d "$TEST_SCRATCH" ]; then
    rm -rf "$TEST_SCRATCH"
  fi
}
trap cleanup_test_scratch EXIT

if [ -n "$TEST_SRA" ]; then
  echo "Found local test SRA file: ${TEST_SRA}"
  cp "$TEST_SRA" "${TEST_SCRATCH}/ERR7512327.sra"
else
  echo "Local ERR7512327.sra not found in search paths."
  echo "Attempting to download ERR7512327 (~35 MB) using container prefetch..."
  if docker run --rm -v "${TEST_SCRATCH}:/scratch" -w /scratch "$IMAGE" prefetch ERR7512327 -O ./; then
    if [ -f "${TEST_SCRATCH}/ERR7512327/ERR7512327.sra" ]; then
      mv "${TEST_SCRATCH}/ERR7512327/ERR7512327.sra" "${TEST_SCRATCH}/ERR7512327.sra"
      rmdir "${TEST_SCRATCH}/ERR7512327" 2>/dev/null || true
    fi
  fi

  if [ ! -s "${TEST_SCRATCH}/ERR7512327.sra" ]; then
    warn "Could not find or download ERR7512327.sra. Skipping end-to-end test."
    echo "To run the full test, place ERR7512327.sra in automation/ERR7512327/ or provide internet access to prefetch."
    exit 0
  fi
  pass "Successfully downloaded ERR7512327.sra into test scratch directory"
fi

echo "Executing container simulating HTCondor isolated scratch execution..."
echo "Command: docker run --rm -v ${TEST_SCRATCH}:/scratch -w /scratch ${IMAGE} /pipeline/scripts/run_sample_condor.sh ERR7512327 4 2"
echo "-------------------------------------------------------------------"

docker run --rm \
  -v "${TEST_SCRATCH}:/scratch" \
  -w /scratch \
  "$IMAGE" \
  /pipeline/scripts/run_sample_condor.sh ERR7512327 4 2

echo "-------------------------------------------------------------------"
info "Validating generated outputs and intermediate cleanup"

# Check final result table
SUMMARY_FILE="${TEST_SCRATCH}/filtered_ERR7512327_0.01_5_20_2_34.txt"
if [ -s "$SUMMARY_FILE" ]; then
  N_RECORDS=$(wc -l < "$SUMMARY_FILE" | tr -d ' ')
  pass "Final summary table created: $(basename "$SUMMARY_FILE") (${N_RECORDS} lines)"
else
  fail "Final summary table missing or empty: ${SUMMARY_FILE}"
fi

# Check deletion candidate table
DELETIONS_FILE="${TEST_SCRATCH}/filtered_ERR7512327_annotated_standardizationAligned.deletion.sorted.txt"
if [ -s "$DELETIONS_FILE" ]; then
  pass "Deletion candidates table created: $(basename "$DELETIONS_FILE")"
else
  fail "Deletion candidates table missing or empty"
fi

# Check sorted BAM
BAM_FILE="${TEST_SCRATCH}/filtered_ERR7512327_annotated_standardizationAligned.out.sorted.bam"
if [ -s "$BAM_FILE" ]; then
  pass "Filtered sorted BAM created: $(basename "$BAM_FILE") ($(du -h "$BAM_FILE" | cut -f1))"
else
  fail "Filtered sorted BAM missing or empty"
fi

# Check coverage file
COV_FILE="${TEST_SCRATCH}/filtered_ERR7512327_annotated_standardizationAligned.out.sorted.coverage"
if [ -s "$COV_FILE" ]; then
  pass "Coverage profile created: $(basename "$COV_FILE") ($(wc -l < "$COV_FILE" | tr -d ' ') positions)"
else
  fail "Coverage profile missing or empty"
fi

# Check pipeline log
LOG_FILE="${TEST_SCRATCH}/ERR7512327_pipeline.log"
if [ -s "$LOG_FILE" ]; then
  pass "Pipeline log created: $(basename "$LOG_FILE")"
else
  fail "Pipeline log missing or empty"
fi

# Verify intermediate cleanup
INTERMEDIATES_LEAKED=0
for f in "${TEST_SCRATCH}/ERR7512327.sra" "${TEST_SCRATCH}"/*.fastq "${TEST_SCRATCH}"/*.fq "${TEST_SCRATCH}"/*Aligned.out.sam; do
  if [ -e "$f" ]; then
    warn "Intermediate file was not cleaned up: $(basename "$f")"
    INTERMEDIATES_LEAKED=$((INTERMEDIATES_LEAKED + 1))
  fi
done

if [ "$INTERMEDIATES_LEAKED" -eq 0 ]; then
  pass "Intermediate cleanup verified: raw .sra, fastqs, and SAM were removed as expected"
else
  warn "Some intermediate files remained in scratch directory."
fi

echo ""
echo -e "${GREEN}${BOLD}===================================================================${NC}"
echo -e "${GREEN}${BOLD}ALL TESTS PASSED SUCCESSFULLY!${NC}"
echo -e "Your Docker image '${IMAGE}' is fully verified and ready to run on HTCondor."
echo -e "${GREEN}${BOLD}===================================================================${NC}"
