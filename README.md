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

##  Repository Structure
rna-seq-differential-expression-pipeline/
│── data/ # test FASTQ files or SRA accession IDs
│── scripts/ # preprocessing/QC shell scripts
│── notebooks/ # R notebooks (DESeq2 analysis)
│── workflow/ # Galaxy workflow export (.ga file)
│── results/ # output plots (volcano, PCA, heatmap)
│── environment.yml # Conda environment file (dependencies)
│── README.md # project documentation


---

## ⚙ Requirements
- Python ≥ 3.8  
- R ≥ 4.0 with **DESeq2** installed  
- Galaxy or local tools: FastQC, HISAT2, FeatureCounts  

---

## � Usage
Clone this repository:
```bash
git clone https://github.com/AqsaH999/rna-seq-differential-expression-pipeline.git
cd rna-seq-differential-expression-pipeline

 License
This project is licensed under the MIT License.
See the LICENSE file for details.
 Author: Aqsa Hamdani
 Master's in Bioinformatics | RNA-seq & Differential Expression Analysis

