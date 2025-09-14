#!/bin/bash
# 00_download_sra.sh - Download FASTQ files from SRA
mkdir -p data/raw

# Loop through SRA accessions listed in sample_info.csv
while IFS=, read -r sample r1 r2 condition batch sra; do
  # Skip header
  [[ "$sample" == "sample" ]] && continue

  echo "Downloading $sra..."
  prefetch "$sra"
  fasterq-dump "$sra" --split-files -O data/raw
  gzip data/raw/"$sra"_1.fastq
  gzip data/raw/"$sra"_2.fastq
done < data/sample_info.csv
