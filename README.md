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

## Project Files

- **Data file:** [`data/GSE183947_fpkm.csv`](data/GSE183947_fpkm.csv)  
  FPKM expression matrix used as the primary input for the analyses.

- **Assignment description:** [`task.md`](task/task.md)  
  Contains the detailed bioinformatics assignment tasks and questions to be addressed.

- **Analysis plan:** [`plan.md`](task/plan.md)  
  Step-by-step workflow plan for answering all assignment tasks using R.

- **Main analysis/report:** [`analysis.Rmd`](analysis.Rmd)  
  R Markdown document implementing the analysis pipeline and generating the report.

## Next Steps

- Follow the workflow outlined in `plan.md` to implement the R analysis in `analysis.Rmd`.
- Knit `analysis.Rmd` after restoring the `renv` environment to produce the final report.
- Update this `README.md` with more detailed usage instructions, commands, and examples as the project evolves.

## License
This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
