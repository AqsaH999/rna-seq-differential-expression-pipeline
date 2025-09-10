
# RNA-seq Differential Expression Pipeline

This repository implements a reproducible RNA-seq differential expression workflow using **Galaxy** and **R (DESeq2)**.  
The pipeline processes raw FASTQ reads through quality control, alignment, quantification, and statistical analysis to identify differentially expressed genes.

---

## 📌 Project Overview
- **Input**: Raw FASTQ RNA-seq data (small test dataset from [SRA](https://www.ncbi.nlm.nih.gov/sra) or Galaxy Training Network).  
- **Steps**:
  1. Quality Control (FastQC, Trimmomatic)
  2. Alignment (HISAT2)
  3. Quantification (FeatureCounts)
  4. Differential Expression Analysis (DESeq2 in R)
  5. Visualization (Volcano plot, PCA, Heatmap)

- **Output**: List of differentially expressed genes + plots for interpretation.

---

## 🛠️ Tools & Technologies
- **Linux/Unix** (file handling, pipeline scripting)
- **Galaxy** (workflow execution)
- **R / Bioconductor** (DESeq2, ggplot2, pheatmap)
- **Python** (optional: matplotlib/seaborn for visualization)
- **GitHub** (version control & reproducibility)

---

## 📂 Repository Structure
rna-seq-differential-expression-pipeline/
│── data/ # test FASTQ files or SRA accession IDs
│── scripts/ # shell scripts for preprocessing/QC
│── notebooks/ # R notebooks for DESeq2 analysis
│── workflow/ # Galaxy workflow export (.ga file)
│── results/ # plots (volcano, PCA, heatmap)
│── environment.yml # Conda environment file (dependencies)
│── README.md # project documentation
