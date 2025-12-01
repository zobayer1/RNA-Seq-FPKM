# CSE763 Bioinformatics Assignment – RNA-Seq FPKM Analysis

This repository contains work for a CSE763 bioinformatics assignment focusing on RNA-Seq gene expression analysis using FPKM values from a paired tumor vs normal dataset.

The overall goal is to:
- Load and understand the RNA-Seq FPKM dataset.
- Perform differential expression analysis between tumor and normal samples using R.
- Visualize results (e.g., volcano and MA plots).
- Conduct downstream functional/pathway analysis.

We will refine and expand this README as the project progresses.

## Prerequisites

To work with this repository you should have:

- **R** installed (version 4.x recommended).
- A basic R environment (R console or RStudio).
- Recommended R packages (to be installed as needed during the analysis):
  - `tidyverse` (or at least `readr`, `dplyr`, `ggplot2`)
  - `DESeq2`
  - `AnnotationDbi` and organism annotation packages (e.g., `org.Hs.eg.db`)
  - Optional: `clusterProfiler`, `enrichplot`, `fgsea`, `pheatmap`

## Setup with `renv` (recommended)

This project uses **renv** for per-project dependency management. After cloning the repository, you can restore the package environment recorded in `renv.lock`.

1. Start R from the project root (this folder), then run:
    ```r
    install.packages("renv")  # only if renv is not already installed
    renv::restore()           # installs all packages listed in renv.lock
    ```
2. Once `renv::restore()` completes, the project will use its own library. In future sessions, just start R (or open the project in RStudio) from the project root and `renv` will automatically activate.

## Data Files

### The Raw Data for DESeq2
File: [GSE183947_raw_counts_GRCh38.p13_NCBI.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE183947&format=file&file=GSE183947_raw_counts_GRCh38.p13_NCBI.tsv.gz)

- Required for differential expression analysis.
- Contains raw read counts (integer values) for each gene in each sample.
- DESeq2 requires raw counts, not normalized data like FPKM/TPM.

### The FPKM Data
File: [GSE183947_fpkm.csv](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947)

- Contains FPKM values (floating point numbers) - normalized for gene length and sequencing depth.
- FPKM is for visualization/exploration, not statistical testing.

### Sample Metadata
File: [GSE183947_series_matrix.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183947/matrix/GSE183947_series_matrix.txt.gz)

- Contains the sample information (metadata): which samples are cancer vs. normal.
- This file has both the expression matrix (FPKM values) AND the sample descriptions.
- Needed for extracting the sample conditions (cancer/normal).

### Gene Annotation
File: [Human.GRCh38.p13.annot.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz)

- Contains gene names, symbols, types, etc.
- Useful for converting Ensembl IDs to gene symbols for interpretation.

## Project Files
- **Assignment description:** [`task.md`](task/task.md)  
  Contains the detailed bioinformatics assignment tasks and questions to be addressed.

- **Analysis plan:** [`plan.md`](task/plan.md)  
  Step-by-step workflow plan for answering all assignment tasks using R.

- **Main analysis/report:** [`analysis.Rmd`](analysis.Rmd)  
  R Markdown document implementing the analysis pipeline and generating the report.

## License
This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
