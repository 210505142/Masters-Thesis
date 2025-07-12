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
├── README.md
└── LICENSE
```

- **RNA_Sequencing/**: All files related to RNA sequencing analysis.
  - `data/`: Raw and processed RNA-seq data.
  - `scripts/`: Code for preprocessing, quality control, and downstream RNA-seq analysis.
  - `results/`: Output files and figures from RNA-seq analysis.
- **Proteomics/**: All files related to proteomics analysis.
  - `data/`: Raw and processed proteomics data.
  - `scripts/`: Code for preprocessing, quality control, and downstream proteomics analysis.
  - `results/`: Output files and figures from proteomics analysis.

## Getting Started

### Prerequisites

- [Python 3.x](https://www.python.org/)
- [R](https://www.r-project.org/) (for some analysis scripts)
- Package requirements for each analysis are listed in the respective `requirements.txt` or `environment.yml` files inside each folder.

### Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```

2. Install dependencies for each analysis as needed.

### Usage

- Refer to the `RNA_Sequencing/README.md` and `Proteomics/README.md` for detailed instructions on running analyses in each section.
- Sample commands and workflows are provided in each subfolder.
