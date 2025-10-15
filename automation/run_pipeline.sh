#!/bin/bash
set -e

SRA=$(basename "$PWD")
PRIMER="$1"
REPO="$2"

PRIMER=${1:-${REPO}/ARTIC_primers_v4.bed}

REF=${REPO}/NC_045512.2.fasta

GENOME_DIR=${REPO}/Genome_Dir

echo "Running pipeline for ${SRA} with primer: ${PRIMER}"

name=${SRA}

fasterq-dump --split-files ./${SRA}.sra -O ./

trim_galore --stringency 3 -q 30 -e .10 --length 15 --paired ${SRA}_1.fastq ${SRA}_2.fastq || true

STAR --runThreadN 8 --readFilesIn ./${SRA}_1_val_1.fq ./${SRA}_2_val_2.fq \
  --outFileNamePrefix ${name} \
  --genomeDir ${GENOME_DIR} \
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

python3 ${REPO}/standardize_alignments.py ${REF} ${SRA}Aligned.out.sam \
  > ${SRA}_standardizationAligned.out.bam

python3 ${REPO}/annotate_alignment_with_primers.py ${PRIMER} \
  ${SRA}_standardizationAligned.out.bam \
  > ${SRA}_annotated_standardization_Aligned.out.bam

python3 ${REPO}/filter_aligned_reads.py ${PRIMER} "${SRA}_annotated_standardization_Aligned.out.bam" --min-deletion-length 20 --max-overhang-primer-frac 1 --min-aligned-length 75 --primer-pool-matching --max-primer-dist 1 > "filtered_${SRA}_annotated_standardizationAligned.out.bam"

python3 ${REPO}/extract_deletions.py \
  filtered_${SRA}_annotated_standardizationAligned.out.bam \
  --primer-bed ${PRIMER} --min-deletion-length 20 --ignore-secondary \
  > filtered_${SRA}_annotated_standardizationAligned.deletion.sorted.txt

samtools sort filtered_${name}_annotated_standardizationAligned.out.bam -o filtered_${name}_annotated_standardizationAligned.out.sorted.bam

samtools depth -a -m 0 filtered_${name}_annotated_standardizationAligned.out.sorted.bam > filtered_${name}_annotated_standardizationAligned.out.sorted.coverage

Rscript ${REPO}/deletion_summarization_and_further_filtration.R