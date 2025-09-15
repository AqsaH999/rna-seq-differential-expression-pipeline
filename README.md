# RNA-seq Differential Expression Pipeline

A reproducible RNA-seq differential expression workflow using **Galaxy** and **R (DESeq2)**.  
The pipeline processes raw FASTQ reads through quality control, alignment, quantification, and statistical analysis to identify differentially expressed genes.

---

##  Project Overview
- **Input**: Raw FASTQ RNA-seq data (test dataset from [SRA](https://www.ncbi.nlm.nih.gov/sra) or Galaxy Training Network)  
- **Pipeline Steps**:
  1. Quality Control — *FastQC, Trimmomatic*  
  2. Alignment — *HISAT2*  
  3. Quantification — *FeatureCounts*  
  4. Differential Expression Analysis — *DESeq2 (R)*  
  5. Visualization — *Volcano plot, PCA, Heatmap*  

- **Output**: Differentially expressed gene list + publication-ready plots  

---

## Tools & Technologies
- **Galaxy** — workflow execution  
- **Linux/Unix** — file handling, scripting  
- **R / Bioconductor** — DESeq2, ggplot2, pheatmap  
- **Python** *(optional)* — matplotlib, seaborn for custom visualization  
- **GitHub** — version control & reproducibility  

---


##  Requirements
- Python ≥ 3.8  
- R ≥ 4.0 with **DESeq2** installed  
- Galaxy or local tools: FastQC, HISAT2, FeatureCounts  

---
##  Repository Structure

The repository is organized as follows:

rna-seq-differential-expression-pipeline/

│── data/ # Raw FASTQ files, reference genome, metadata (sample_info.csv)

│── scripts/ # Shell scripts for QC, trimming, alignment, counting

│── notebooks/ # R/Python notebooks for analysis & visualization

│── workflow/ # Galaxy workflow exports (.ga files) [optional]

│── results/ # Output files: count matrix, plots (PCA, volcano, heatmap)

│── CITATION.cff # Citation metadata for this repository

│── LICENSE # Open-source license (MIT)

│── README.md # Project documentation

│── .gitignore # Ignored files for Python & R

│── environment.yml # Conda environment specification


##  Usage

Clone this repository and navigate to the project directory:

```bash
git clone https://github.com/your-username/rna-seq-differential-expression-pipeline.git
cd rna-seq-differential-expression-pipeline

1. Set up environment
Install required tools and R packages via Conda:
conda env create -f environment.yml
conda activate rnaseq-pipeline

2. Prepare input data
Place your raw FASTQ files in the data/ folder.
Update data/sample_info.csv with sample IDs, conditions, and file paths.

3. Run pipeline scripts
Execute the workflow step by step:
bash scripts/01_qc.sh         # Run FastQC on raw reads
bash scripts/02_trim.sh       # Trim adapters & low-quality bases
bash scripts/03_align.sh      # Align reads to reference genome (HISAT2/STAR)
bash scripts/04_counts.sh     # Generate gene counts (featureCounts/HTSeq)

4. Differential expression analysis
Use the R notebook for DESeq2 analysis:
Rscript notebooks/deseq2_analysis.R
This will produce:
Normalized count matrix
PCA plot of samples
Volcano plot of DEGs
Heatmap of top genes

5.View results
All output files will be saved in the results/ directory.
