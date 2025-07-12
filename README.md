# Masters-Thesis
Accessibility for all the scripts used for generation of plots for both Proteomics and RNA seq
## Organisation
- RNA seq file contains full pipeline of quantification of data as well as the scripts for the visualisation
- Proteomics file contains 
# RNA Sequencing & Proteomics Analysis Repository

Welcome to the **RNA Sequencing & Proteomics Analysis** repository!  
This repository contains scripts, notebooks, and documentation for performing separate analyses on RNA-seq and Proteomics datasets.

## Repository Structure

```
.
├── RNA_Sequencing/
│   ├── data/
│   ├── scripts/
│   └── results/
├── Proteomics/
│   ├── data/
│   ├── scripts/
│   └── results/
|── README.md
```

- **RNA_Sequencing/**: All files related to RNA sequencing analysis.
  - `data/`: Raw and processed RNA-seq data.
  - `scripts/`: Code for preprocessing, quality control, and downstream RNA-seq analysis.
  
- **Proteomics/**: All files related to proteomics analysis.
  - `data/`: processed proteomics data.
  - `scripts/`: Code for  quality control, and downstream proteomics analysis.

## Getting Started

### Prerequisites

- [R](https://www.r-project.org/) (for some analysis scripts)
- Package requirements for each analysis are listed in 'Software'.

### Software required
RNA Seq:
- HISAT2 v2.1.0
- Salmon v1.3.0
- Fast QC v0.11.7
- ClusterProfiler v4.6.2
- GSEA v4.3.3
- tidyverse v2.0.0
- Glimma
- **Bioconductor** 3.1.8:
- DESeq2 v1.48.1
- vsn v3.77.0
- rtracklayer
- GenomicFeatures
- reactome.db
- tximport v1.26.1
- variancePartition
- fgsea
- EnhancedVolcano

Proteomics
- DIA-nn v1.8.1/2.0.1
- limma v3.1.0
- ClusterProfiler v4.6.2
- tidyverse v2.0.0
- Bioconductor 3.1.8

### Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```

2. Install dependencies for each analysis as needed.


- Refer to the `RNA_Sequencing/README.md` and `Proteomics/README.md` for detailed instructions on running analyses in each section.
- Sample commands and workflows are provided in each subfolder.
