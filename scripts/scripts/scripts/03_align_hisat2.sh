#!/bin/bash
# 03_align_hisat2.sh - Align reads to reference genome
mkdir -p results/bam

# Build index (do this only once)
hisat2-build ref/genome.fa ref/genome_index

for r1 in data/trimmed/*_1_paired.fastq.gz; do
  base=$(basename "${r1}" "_1_paired.fastq.gz")
  r2="data/trimmed/${base}_2_paired.fastq.gz"
  sam="results/bam/${base}.sam"
  bam="results/bam/${base}.sorted.bam"

  hisat2 -p 4 -x ref/genome_index -1 "${r1}" -2 "${r2}" -S "${sam}"
  samtools view -@ 4 -bS "${sam}" | samtools sort -@ 4 -o "${bam}"
  samtools index "${bam}"
  rm "${sam}"
done
