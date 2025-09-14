#!/bin/bash
# 01_qc.sh - Run FastQC and MultiQC on raw FASTQ files
mkdir -p results/qc

fastqc -t 4 -o results/qc data/raw/*.fastq.gz
multiqc results/qc -o results/qc
