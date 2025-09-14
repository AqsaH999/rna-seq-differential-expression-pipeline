
# RNA-seq Differential Expression Pipeline

This repository implements a reproducible RNA-seq differential expression workflow using **Galaxy** and **R (DESeq2)**.  
The pipeline processes raw FASTQ reads through quality control, alignment, quantification, and statistical analysis to identify differentially expressed genes.



## 📌 Project Overview
- **Input**: Raw FASTQ RNA-seq data (small test dataset from [SRA](https://www.ncbi.nlm.nih.gov/sra) or Galaxy Training Network).  
- **Steps**:
  1. Quality Control (FastQC, Trimmomatic)
  2. Alignment (HISAT2)
  3. Quantification (FeatureCounts)
  4. Differential Expression Analysis (DESeq2 in R)
  5. Visualization (Volcano plot, PCA, Heatmap)

- **Output**: List of differentially expressed genes + plots for interpretation.



## 🛠️ Tools & Technologies
- **Linux/Unix** (file handling, pipeline scripting)
- **Galaxy** (workflow execution)
- **R / Bioconductor** (DESeq2, ggplot2, pheatmap)
- **Python** (optional: matplotlib/seaborn for visualization)
- **GitHub** (version control & reproducibility)



## 📂 Repository Structure
rna-seq-differential-expression-pipeline/
│── data/ # test FASTQ files or SRA accession IDs
│── scripts/ # shell scripts for preprocessing/QC
│── notebooks/ # R notebooks for DESeq2 analysis
│── workflow/ # Galaxy workflow export (.ga file)
│── results/ # plots (volcano, PCA, heatmap)
│── environment.yml # Conda environment file (dependencies)
│── README.md # project documentation


---

## Requirements
- Python ≥ 3.8
- R ≥ 4.0 with DESeq2
- Galaxy or local tools (FastQC, HISAT2, featureCounts)

## Usage
```bash
git clone https://github.com/your-username/your-repo.git
cd your-repo

### 2. Add `.gitignore` (Python + R)

Create a `.gitignore` file in the root:

```gitignore
# Python
__pycache__/
*.pyc
*.pyo
*.pyd
*.pkl
*.egg-info/
.env
.venv/
build/
dist/
.ipynb_checkpoints/

# R
.Rhistory
.RData
.Rproj.user/
.Ruserdata
*.Rproj
*.Rhistory
*.RData
*.Rout
.Rapp.history

# OS / Misc
.DS_Store
Thumbs.db

MIT License

Copyright (c) 2025 Aqsa Hamdani

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.

