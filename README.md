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

**Note (project decision):** During preliminary inspection, the available raw read-count information for this dataset was found to be incomplete and did not provide a single, consistent counts matrix covering all tumor and normal samples. As a result, standard count-based workflows (e.g., DESeq2 on a full raw counts matrix) are not applied here. Instead, this analysis uses the FPKM expression table together with a limma-based approach on log2(FPKM+1) values for the differential expression step, and this limitation is explicitly discussed in the analysis.

### FPKM Expression Data (used for visualization and limma-based DE)
File: `GSE183947_fpkm.csv`

- Contains FPKM values (floating point numbers) – normalized for gene length and sequencing depth.
- Used in this project for visualization, exploration, and differential expression using a limma-based approach on log2(FPKM+1).
- Note: Using normalized expression (FPKM) for differential testing is not ideal. This approach is a pragmatic alternative when a suitable complete raw counts matrix is not available; interpretation should mention this limitation.

### Sample Metadata
File: `GSE183947_series_matrix.txt`

- Contains the sample information (metadata): which samples are cancer vs. normal.
- This file has the sample descriptions and is used to extract sample condition (cancer/normal) and sample IDs for aligning with the expression matrix.

### Gene Annotation
File: `Human.GRCh38.p13.annot.tsv`

- Contains gene names, symbols, types, etc.
- Useful for converting gene IDs/symbols in the expression table to standardized symbols for interpretation and downstream enrichment analysis.

## Project Files
- **Assignment description:** [`task.md`](task/task.md)  
  Contains the detailed bioinformatics assignment tasks and questions to be addressed.

- **Analysis plan:** [`plan.md`](task/plan.md)  
  Step-by-step workflow plan for answering all assignment tasks using R.

- **Main analysis/report:** [`analysis.Rmd`](analysis.Rmd)  
  R Markdown document implementing the analysis pipeline and generating the report.

## License
This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
