#!/bin/bash
# 02_trim.sh - Paired-end trimming with Trimmomatic
mkdir -p data/trimmed

for r1 in data/raw/*_1.fastq.gz; do
  base=$(basename "${r1}" "_1.fastq.gz")
  r2="data/raw/${base}_2.fastq.gz"
  out1="data/trimmed/${base}_1_paired.fastq.gz"
  out2="data/trimmed/${base}_2_paired.fastq.gz"
  out1u="data/trimmed/${base}_1_unpaired.fastq.gz"
  out2u="data/trimmed/${base}_2_unpaired.fastq.gz"

  trimmomatic PE -threads 4 \
    "${r1}" "${r2}" \
    "${out1}" "${out1u}" "${out2}" "${out2u}" \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
