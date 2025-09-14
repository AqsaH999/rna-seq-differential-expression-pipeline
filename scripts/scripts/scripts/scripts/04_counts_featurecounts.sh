#!/bin/bash
# 04_counts_featurecounts.sh - Generate gene counts using featureCounts
mkdir -p results/counts

featureCounts -T 4 -p -t exon -g gene_id -a ref/annotation.gtf -o results/counts/featurecounts.txt results/bam/*.sorted.bam
