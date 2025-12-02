# RNA-Seq FPKM Differential Expression & Downstream Analysis Plan

This plan outlines a detailed, step-by-step workflow for solving all 4 problems in the assignment using R and the datasets in the `data/` directory:
- FPKM-normalized expression values stored in `GSE183947_fpkm.csv` for exploration/visualization **and FPKM-based DE analysis (using limma on log2(FPKM+1)) as the primary approach due to the lack of a complete raw counts matrix for all samples**.
- GEO series matrix `GSE183947_series_matrix.txt` with sample metadata (e.g., tumor vs normal, pairing).
- Gene annotation table `Human.GRCh38.p13.annot.tsv` for mapping gene IDs/symbols to standardized symbols and other attributes.

It focuses on what needs to be studied and done, not on the final solutions.

---

## 1. Global Preparation and Background

### 1.1 Understand the Assignment and Data
- Read `task.md` carefully to identify the 4 specific questions/problems and any constraints (e.g., tools, packages, or plots required).
- Confirm how each file in `data/` is used:
  - FPKM expression file (`GSE183947_fpkm.csv`) → primary expression matrix used for exploratory plots / QC **and as the main quantitative input for limma-based differential expression**.
  - GEO series matrix (`GSE183947_series_matrix.txt`) → sample-level metadata, including which samples are tumor vs normal and how they are paired.
  - Gene annotation (`Human.GRCh38.p13.annot.tsv`) → mapping from gene IDs/symbols to human-readable gene symbols and possibly Entrez IDs.
- Note that no single, complete raw counts matrix covering all tumor and normal samples is available for this dataset; therefore, standard count-based DE workflows (e.g., DESeq2 on raw counts) cannot be applied to the full cohort.
- Identify which tasks correspond to:
  - Data loading and understanding.
  - Differential expression analysis.
  - Visualization (MA plot, volcano plot, etc.).
  - Downstream analysis (e.g., pathway/enrichment analysis).

### 1.2 Study Required Biological and Statistical Concepts
- Review RNA-Seq basics:
  - Sequencing reads → alignment → read counts → normalization → differential expression.
- Study what FPKM means:
  - FPKM = Fragments Per Kilobase of transcript per Million mapped reads.
  - Understand that FPKM is normalized within each sample and has limitations for between-sample DE testing compared with raw counts.
- Understand core differential expression (DE) concepts:
  - Null hypothesis: no difference in expression between tumor and normal.
  - Log2 fold change (log2FC) and what positive/negative values mean.
  - p-values and multiple testing corrections (FDR, Benjamini–Hochberg adjustment).
- Review experimental design concepts:
  - Unpaired design: comparing tumor vs normal as independent groups.
  - Paired design: accounting for patient-specific effects using patient/pair ID.

### 1.3 Set Up R Environment
- Ensure R and an IDE (e.g., RStudio) or terminal R session work correctly on Linux.
- Install core CRAN packages:
  - `tidyverse` (or at least `readr`, `dplyr`, `tidyr`, `ggplot2`) for data handling and plotting.
  - `data.table` (optional) for fast data loading.
  - `pheatmap` or `ComplexHeatmap` (optional) for QC/heatmaps.
- Install Bioconductor packages:
  - `DESeq2` for differential expression analysis.
  - `AnnotationDbi` and appropriate organism database (likely `org.Hs.eg.db` for human data).
  - Optionally `clusterProfiler`, `enrichplot`, `fgsea` for enrichment analysis in R.
- Verify installation by loading all required libraries without errors.

### 1.4 Organize Project Structure
- Use a clear directory layout:
  - `data/` (already exists) for input data.
  - `scripts/` for R scripts (e.g., `01_load_data.R`, `02_deseq_analysis.R`).
  - `results/` for intermediate outputs (e.g., DE result tables).
  - `figures/` for plots (MA plot, volcano plot, enrichment plots).
- Decide on main analysis entry points:
  - A single `analysis.R` script or an R Markdown file combining all steps for reproducibility.
- Decide on file naming conventions for outputs (e.g., include date, thresholds in filenames).

---

## 2. Task 1 – Load and Understand the RNA-Seq Data

### 2.1 Study Needed R and Data Concepts
- R data import and inspection:
  - How to use `readr::read_csv()`, `readr::read_tsv()`, or `data.table::fread()` to read large CSV/TSV files.
  - Functions to inspect data: `head()`, `str()`, `summary()`, `dim()`, `colnames()`, `rownames()`.
- Data structures for expression data:
  - Wide matrix format with genes in rows and samples in columns.
  - Distinction between expression matrix and sample metadata (phenotype/colData).
- Sample metadata (colData):
  - Understand that DE tools (e.g., DESeq2) expect a `colData` data frame with one row per sample.
  - Columns to include: sample ID, condition (tumor/normal), patient ID (for pairing), and any other covariates if required.
- Gene annotation:
  - Understand that a separate table (here `Human.GRCh38.p13.annot.tsv`) provides mapping from gene IDs to symbols/biotypes.
  - Know when to apply annotation (often after DE, but basic checks can be done during data understanding).

### 2.2 Inspect and Parse the FPKM CSV File
- Load `data/GSE183947_fpkm.csv` into R.
- Examine structure:
  - Determine which column(s) contain gene identifiers (e.g., gene symbol, Ensembl ID).
  - Identify how sample columns are named and how many samples there are.
  - Check if samples appear paired based on naming (e.g., normal/tumor suffixes).
- Check data quality at a basic level:
  - Look for missing values (NA) and consider how to handle them.
  - Check for duplicated gene IDs and decide how to resolve duplicates (e.g., keep max, mean, or first occurrence).
  - Look at ranges of FPKM values and presence of zeros/near-zeros.

### 2.3 Briefly Consider Raw Counts (conceptual, no full matrix)
- If raw read-count files or subsets are provided in the original resource, inspect their structure conceptually (rows as genes, columns as samples) and assess whether they form a complete matrix covering all tumor and normal samples.
- Document intended/limited use:
  - Note that for this assignment, the available raw read-count information is incomplete and does **not** provide a unified counts matrix for all 60 samples and both conditions.
  - Because of this, it is **not** used as the primary input for tumor vs normal DE analysis in this plan.
  - Instead, differential expression will be based on log2(FPKM+1) with limma (Task 2), and this limitation will be clearly documented in the analysis and report.

### 2.4 Build Expression Matrix, Sample Metadata, and Annotation Checks
- Construct the expression matrix from FPKM:
  - Set row names to gene identifiers (gene symbols from the first column of `GSE183947_fpkm.csv`).
  - Ensure all expression columns are numeric.
  - Confirm that rows = genes and columns = samples.
- Load sample metadata from GEO series matrix:
  - Read `data/GSE183947_series_matrix.txt`.
  - Extract sample IDs and relevant characteristics (e.g., tumor vs normal status, any available pairing/patient information).
- Create sample metadata (`colData`):
  - Derive a `sample_id` column matching the column names of the FPKM matrix.
  - Create a `condition` factor with levels `normal` and `tumor` based on the series matrix.
  - If possible, create a `patient_id` to represent paired samples (same patient in normal and tumor).
- Load and check gene annotation (`Human.GRCh38.p13.annot.tsv`):
  - Identify key columns, in particular the `Symbol` column used to map FPKM gene symbols.
  - Check compatibility by verifying what proportion of FPKM gene symbols are found in the annotation table, using both `Symbol` and `Synonyms` (including synonyms separated by `,`, `;`, or `|`).
  - Note any systematic differences or missing symbols and how they will be handled later (e.g., treating them as unannotated in DE results and enrichment).
- Verify consistency:
  - Confirm that the order of samples in `colData` matches the columns in the FPKM expression matrix.
  - Summarize the number of genes and samples and the condition breakdown (normal vs tumor).

### 2.5 Basic Exploratory Data Analysis (EDA)
- Compute simple summaries using the FPKM expression matrix:
  - Number of genes and samples.
  - Distribution of FPKM values for a few random samples.
- Optional quick visualizations for understanding:
  - Boxplots of log-transformed FPKM per sample to assess overall distribution.
  - Sample-sample correlation heatmap or PCA plot for global patterns (e.g., do tumors cluster separately from normals?).
- Cross-check metadata:
  - Summarize the number of `normal` vs `tumor` samples.
  - If `patient_id` is available, verify that each patient has both a normal and a tumor sample.

### 2.6 Document Assumptions and Limitations
- Clearly note:
  - The nature of the data (FPKM vs raw counts) and that a suitable, complete raw counts matrix for both tumor and normal samples is not available, so **DESeq2-based DE on counts is not feasible for the full dataset**.
  - As a result, differential expression will be performed using log2(FPKM+1) with limma as the main method, and this choice should be explicitly discussed as a limitation.
  - Any preprocessing steps applied at this stage (filtering low-expression genes, handling missing values, resolving duplicate IDs).
  - Any assumptions about pairing inferred from sample names or metadata in the series matrix.
  - Any ID-format or naming issues discovered between FPKM and annotation and how you plan to handle them.

---

## 3. Task 2 – Differential Gene Expression Analysis (FPKM + limma)

### 3.1 Study Statistical and Methodological Foundations
- limma fundamentals (FPKM-based DE as the main approach here):
  - Understand that limma was originally developed for microarray-like log-intensity data but is widely used with voom-transformed counts.
  - Recognize that, in this assignment, no complete raw counts matrix covering all tumor and normal samples is available. Therefore, DESeq2-based tumor vs normal DE on counts is not possible for all 60 samples.
  - Plan to use **log2(FPKM + 1)** with limma as the DE method, and clearly document this limitation in the report.
- Design formulas:
  - Simple design: `~ condition` (unpaired comparison of tumor vs normal).
  - Paired design: `~ patient_id + condition` to account for matched pairs, using the pairing information from the series matrix.
  - Understand how limma interprets factors and how the reference level of `condition` determines the sign of the estimated log2 fold changes.
- Multiple testing and thresholds:
  - Review Benjamini–Hochberg FDR control.
  - Decide on significance thresholds (e.g., `adj.P.Val < 0.05` or more stringent).
  - Decide on a biologically meaningful log2FC cutoff (e.g., |log2FC| > 1).

### 3.2 Prepare Data for DE Analysis (FPKM + limma)
- Ensure the FPKM-based expression matrix is in a suitable format for limma:
  - Use the `expr_fpkm` matrix (genes × samples) constructed in Task 1.
  - Compute log2-transformed expression (e.g., `log2(FPKM + 1)`) to obtain approximately homoscedastic data for linear modeling.
- Finalize `colData`:
  - Confirm `condition` is a factor and choose the reference level (e.g., `normal` as reference so positive log2FC means higher in tumor).
  - Include `patient_id` as a blocking factor if a paired analysis is desired.
- Filter genes (optional but common):
  - Remove genes with extremely low FPKM across all samples (e.g., near-zero expression) according to a chosen rule (e.g., require a minimum FPKM or log2(FPKM+1) threshold in a minimum number of samples).

### 3.3 Run Differential Expression with limma
- Construct the design matrix using `model.matrix()`:
  - Unpaired design: `design <- model.matrix(~ condition, data = colData)`.
  - Paired design: `design <- model.matrix(~ patient_id + condition, data = colData)`.
- Fit the linear model using limma:
  - Use `lmFit(log_expr_fpkm, design)` followed by `eBayes()`.
- Extract results for the tumor vs normal contrast:
  - Use `topTable()` (or `coef`/`contrast.fit` if needed) to obtain log2FC, p-values, and adjusted p-values for the condition coefficient.
- Inspect overall results:
  - Summarize numbers of significantly up- and down-regulated genes at chosen thresholds.
  - Examine distributions of log2FC and adjusted p-values.

### 3.4 Post-Processing and Annotation of DE Results
- Apply significance and fold-change thresholds:
  - Create logical flags for significant genes (e.g., `adj.P.Val < threshold`).
  - Create a `regulation` label: `up`, `down`, or `not_significant` based on sign and magnitude of log2FC.
- Annotate genes:
  - If gene identifiers are not human-readable, map them to gene symbols using `Human.GRCh38.p13.annot.tsv` (and/or organism-specific annotation packages).
  - Optionally add gene descriptions if available.
- Save DE results:
  - Export a complete results table (including all genes) to `results/` as a CSV.
  - Optionally, export a filtered table of only significant genes.

### 3.5 Quality Assessment of DE Analysis
- Check diagnostic outputs:
  - Plot p-value histogram and/or volcano plot to assess the distribution of test statistics.
  - Review the number of DE genes; too many or too few may indicate issues.
- Document methodological choices and limitations:
  - Clearly state that DE was performed using log2(FPKM+1) with limma because a full raw counts matrix for both tumor and normal samples was not available.
  - Design formula used and rationale (paired vs unpaired).
  - Threshold values for adjusted p-value and log2FC.
  - Any filters applied before DE analysis.

---

## 4. Task 3 – Visualization: Volcano Plot and MA Plot

### 4.1 Study Visualization Principles
- [x] Understand volcano plot axes: x = log2 fold change (tumor vs normal), y = -log10(adjusted p-value).
- [x] Understand MA plot axes: x = mean expression (e.g., average log2(FPKM+1)), y = log2 fold change.
- [x] Decide how to represent significance (color by `regulation`, reference lines for `alpha` and `logfc_cutoff`).

### 4.2 Prepare Data for Plots (from limma FPKM DE results)
- [x] Use the `res_limma_df` table from Task 2 as the single source of DE results.
- [x] Confirm that `logFC`, `adj.P.Val`, `AveExpr`, and `regulation` in `res_limma_df` are defined based on the limma model on log2(FPKM+1) and thresholds `alpha` and `logfc_cutoff`.
- [x] Derive `-log10(adj.P.Val)` for the volcano plot as `neg_log10_adjP`.

### 4.3 Create Volcano Plot in R (FPKM + limma)
- [x] Use `ggplot2` to create a volcano plot with x = `logFC`, y = `-log10(adj.P.Val)`, colored by `regulation`.
- [x] Add vertical lines at `x = ±logfc_cutoff` and a horizontal line at `y = -log10(alpha)` to show significance thresholds consistent with Task 2.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, add labels for a small number of highly significant genes (e.g., using `ggrepel`) if required by the assignment.
- [x] Save the volcano plot to `figures/` with a clear filename (e.g., `volcano_limma_fpkm_tumor_vs_normal.png`).

### 4.4 Create MA Plot in R (FPKM + limma)
- [x] Use `ggplot2` to create an MA plot with x = `AveExpr` (average log2 expression) and y = `logFC`, colored by `regulation`.
- [x] Add a horizontal reference line at `y = 0` to indicate no change in expression.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, highlight or label significantly differentially expressed genes if desired.
- [x] Save the MA plot to `figures/` with a descriptive filename (e.g., `ma_limma_fpkm_tumor_vs_normal.png`).

### 4.5 Interpret and Document Visualizations
- [ ] Write short, clear narrative text in the report explaining what the volcano plot reveals about the number, direction, and significance of differentially expressed genes.
- [ ] Write short, clear narrative text explaining what the MA plot reveals about how fold changes depend on average expression and whether any systematic patterns or biases are apparent.
- [ ] Explicitly mention that these visualizations are based on limma DE results using log2(FPKM+1) as input rather than DESeq2 on raw counts, due to the lack of a complete counts' matrix.

---

## 5. Task 4 – Downstream Functional and Pathway Analysis

### 5.1 Study Biological / Bioinformatics Background
- Gene sets and pathways:
  - Understand what GO terms (BP/MF/CC), KEGG pathways, Reactome pathways, and MSigDB collections represent.
  - Distinguish between over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
- Enrichment concepts:
  - ORA: uses a list of significant genes vs a background set.
  - GSEA: uses a ranked list of all genes with a continuous statistic.
  - Key outputs: enrichment score, normalized enrichment score (NES), p-value, FDR/q-value.
- Gene identifier types:
  - Study mapping between Ensembl IDs, Entrez IDs, and gene symbols.
  - Learn why consistent ID types are required for each tool.

### 5.2 Prepare Gene Lists for Enrichment
- From DE results, define:
  - List of significantly up-regulated genes (e.g., `padj < threshold` and log2FC > cutoff).
  - List of significantly down-regulated genes (e.g., `padj < threshold` and log2FC < -cutoff).
- Optionally create a ranked gene list:
  - Rank genes by log2FC or a signed statistic (e.g., sign(log2FC) * -log10(pvalue)).

### 5.3 Map Gene Identifiers
- Determine which ID type the enrichment tools expect:
  - DAVID, KEGG, GO tools, or R packages like `clusterProfiler` may require Entrez IDs or gene symbols.
- Use annotation resources:
  - Use `AnnotationDbi` and `org.Hs.eg.db` (if human) to map from gene IDs in the dataset to required IDs.
  - Track how many genes successfully map and how many are lost.

### 5.4 Choose Enrichment Tools and Strategy
- Option A: Web-based tools (e.g., DAVID, Enrichr, web-based GO/KEGG portals):
  - Prepare text files or copy-paste gene lists of up- and down-regulated genes.
  - Select species (likely Homo sapiens).
  - Choose annotation categories (GO BP/MF/CC, KEGG, Reactome, etc.).
  - Download or record top enriched terms (p-value, FDR, gene counts).
- Option B: R-based tools (preferred for reproducibility if allowed/expected):
  - Use `clusterProfiler` for functions like `enrichGO`, `enrichKEGG`, `enricher`, or `GSEA`.
  - Provide gene lists and universe/background as all expressed genes.
  - Optionally use MSigDB gene sets with `fgsea` or similar packages.

### 5.5 Run Enrichment / Pathway Analysis
- Run ORA for up-regulated genes:
  - Identify significantly enriched GO terms and pathways.
- Run ORA for down-regulated genes:
  - Identify enriched terms/pathways potentially associated with suppressed functions in tumors.
- If performing GSEA:
  - Use ranked gene list and chosen gene sets (e.g., Hallmark, KEGG).
  - Extract top positively and negatively enriched pathways.

### 5.6 Summarize and Visualize Enrichment Results
- Visualization:
  - Create bar plots or dot plots of top enriched GO terms/pathways.
  - If using R tools, consider enrichment maps or network plots (optional, depending on assignment scope).
- Interpretation:
  - Identify major biological themes (e.g., cell cycle, immune response, apoptosis, metabolism).
  - Distinguish processes/pathways enriched among up-regulated vs down-regulated genes.
- Documentation:
  - Record thresholds, databases, and tools used for enrichment.
  - Note limitations (e.g., incomplete mapping, redundancy of GO terms, FPKM-based DE).

---

## 6. Integration, Reporting, and Final Checks

### 6.1 Plan the Structure of the Final Report / Answers
- Organize content around the 4 assignment tasks:
  - Task 1: Data loading and understanding.
  - Task 2: Differential expression analysis.
  - Task 3: Visualization (volcano & MA plots).
  - Task 4: Downstream functional and pathway analysis.
- For each task, plan to include:
  - Brief description of methods used.
  - Key results (tables, counts, or summaries).
  - Representative figures if applicable.
  - Short interpretation statements answering the specific question.

### 6.2 Ensure Reproducibility
- Keep all R commands in one or more scripts or an R Markdown file.
- Save intermediate objects (e.g., DE results) if re-running analyses is expensive.
- Set a random seed where needed (e.g., for GSEA implementations that involve randomness).

### 6.3 Final Sanity Checks
- Verify that all 4 assignment questions are explicitly answered using the planned workflow output.
- Cross-check that:
  - Data loaded and preprocessed correctly.
  - DE analysis steps and thresholds match what is described in the text.
  - Plots are consistent with DE results.
  - Enrichment results logically follow from the DE gene lists.
- Note any caveats and limitations clearly (especially the use of FPKM, not raw counts).

---

## 7. Report Writing Guide (5 Main Sections)

This section describes how to structure the written report around the 4 analysis tasks and how to integrate results, figures, and references.

### 7.1 Suggested Report Structure (5 Sections)

1. **Introduction**
   - Brief biological and experimental context:
     - What dataset GSE183947 represents (tumor vs normal, tissue type, number of samples/patients).
     - High-level goal of the analysis (identify differentially expressed genes and pathways between tumor and normal).
   - Motivation:
     - Why differential expression and pathway analysis are useful in this context.
   - Short overview of the workflow:
     - 1–2 sentences summarizing data loading, DE analysis, visualization, and downstream enrichment.

2. **Materials and Methods**
   - **Data description**:
     - Source: `GSE183947_fpkm.csv` (mention GEO if relevant).
     - Variables: what rows (genes) and columns (samples and any metadata) represent.
   - **Preprocessing**:
     - How data were read into R (tools like `readr`, `data.table`).
     - Any filtering of low-expression genes.
     - Handling of missing values or duplicated IDs.
   - **Differential expression analysis**:
     - Statistical model and design formula (e.g., `~ condition` or `~ patient_id + condition`).
     - Software and versions (e.g., R version, DESeq2 version) if used.
     - Thresholds for calling genes differentially expressed (FDR and log2FC cutoffs).
   - **Visualization**:
     - Methods used for volcano and MA plots (e.g., `ggplot2`, `plotMA` from DESeq2).
   - **Downstream/Enrichment analysis**:
     - Tools (e.g., DAVID, clusterProfiler) and gene sets/databases used (GO, KEGG, MSigDB).
     - How gene lists were defined (up/down lists, ranked lists).

3. **Results**
   - **Task 1 – Data understanding**:
     - Number of genes and samples.
     - Summary statistics of expression (e.g., distribution of FPKM, any outliers).
     - Any major QC findings (e.g., clustering of tumor vs normal).  
   - **Task 2 – Differential expression**:
     - Total number of significantly up- and down-regulated genes at chosen thresholds.
     - Maybe a short table of top N DE genes with log2FC and FDR.
   - **Task 3 – Visualization**:
     - Include volcano plot(s) with a brief description of what they show.
     - Include MA plot with commentary on how fold changes vary with expression level.
   - **Task 4 – Downstream analysis**:
     - Tables or bullet lists of top enriched GO terms/pathways for up-regulated and down-regulated genes.
     - 1–2 plots summarizing enrichment (e.g., bar plot of top GO terms).

4. **Discussion**
   - Interpret the main findings:
     - What biological processes/pathways appear activated in tumor vs normal.
     - How these results relate to general cancer biology (e.g., proliferation, apoptosis, immune response).
   - Reflect on consistency between DE results, plots, and enrichment:
     - Do the enriched pathways match the direction and magnitude seen in the DE results and plots?
   - Discuss limitations:
     - Use of FPKM (not raw counts) for DE.
     - Any potential batch effects or unmodeled covariates.
     - Limited sample size, multiple testing issues.

5. **Conclusion and Future Work**
   - Concise summary of the most important results:
     - 2–3 sentences on key DE patterns and enriched pathways.
   - Brief statement on biological implications.
   - Optional: suggestions for follow-up analyses (e.g., validation experiments, additional datasets, deeper pathway exploration).

### 7.2 Integrating Outputs from the 4 Tasks into the Report

- **From Task 1 (Data loading & understanding):**
  - Use key descriptive statistics and any QC plots in the *Materials and Methods* (data description) and *Results* (data overview).
  - Report exact numbers: how many genes, how many samples in each group (tumor vs normal), how pairs were defined.
- **From Task 2 (Differential expression):**
  - Present summary counts of DE genes and top DE genes in the *Results* section.
  - In *Methods*, clearly document the design formula and thresholds used.
  - Use DE tables to support interpretations in the *Discussion* (e.g., highlight important DE genes with known roles in cancer).
- **From Task 3 (Volcano & MA plots):**
  - Include selected figures in the *Results* section with informative captions (what the axes represent, thresholds, and main patterns).
  - Refer to these plots in the text when describing the distribution of DE genes (e.g., “As shown in Figure 1…”).
- **From Task 4 (Enrichment / pathway analysis):**
  - Place tables/plots of enriched pathways in the *Results* section.
  - Use these results heavily in the *Discussion* to build a biological narrative (e.g., “Pathways related to cell cycle and DNA repair were enriched among up-regulated genes, indicating…”).
- Ensure cross-referencing:
  - Number tables and figures (e.g., Table 1, Figure 1) and refer to them by number in the text.
  - Make sure all key outputs from the 4 tasks have a clear place and purpose in the narrative.

### 7.3 Adding Study References to the Report

- **Reference types to include:**
  - Primary publications for the dataset (GSE183947 original study, if available).
  - Methodological references:
    - RNA-Seq and DE methodology papers.
    - DESeq2 or other tools used.
  - Biological background sources:
    - Review articles or textbook chapters describing the relevant cancer type, pathways, or biological processes.
- **Referencing style:**
  - Follow the citation style specified in your course/assignment (e.g., numeric [1], author-year (Smith et al., 2020)).
  - Be consistent throughout the report.
- **Where to cite:**
  - Introduction:
    - Cite sources for biological background and dataset description.
  - Methods:
    - Cite software/tools and methods (e.g., DESeq2, enrichment tools, gene set databases).
  - Discussion:
    - Cite studies that support or contrast your findings.
- **Practical tips:**
  - Keep a small `References` or `Bibliography` section at the end with complete citations.
  - As you work on Tasks 2–4, note down any relevant papers or manuals you consult so they can be cited properly later.

---

## 8. Task Checklist for Progress Tracking

Use this checklist to track your progress through the plan. You can duplicate it into your notes and mark items as you complete them.

### 8.1 Global Preparation
- [x] Install/verify R and required packages (`tidyverse`, `DESeq2`, annotation and enrichment packages).
- [x] Set up project structure (`scripts/`, `results/`, `figures/`).
- [x] Review RNA-Seq, FPKM, and basic DE concepts.
- [x] Understand the assignment requirements and download required data files.

### 8.2 Task 1 – Data Loading and Understanding
- [x] Load `data/GSE183947_fpkm.csv` into R and inspect its structure.
- [x] Load `data/GSE183947_series_matrix.txt` and extract sample metadata (tumor/normal status, pairing if available).
- [x] Load `data/Human.GRCh38.p13.annot.tsv` and verify gene ID compatibility with FPKM symbols using the `Symbol` and `Synonyms` columns.
- [x] Construct the FPKM-based expression matrix (genes in rows, samples in columns) and corresponding `colData` (`sample_id`, `condition`, `patient_id` if possible).
- [x] Perform basic EDA (summaries, distributions, PCA/boxplots) using FPKM expression and sample metadata.
- [x] Document assumptions (e.g., how pairing was inferred, handling of missing values, low-expression filtering, ID mapping decisions, and the decision to not use the partial raw counts file for DE).

### 8.3 Task 2 – Differential Gene Expression Analysis (FPKM + limma)
- [x] Decide on the statistical design for DE (unpaired `~ condition` vs paired `~ patient_id + condition`) based on metadata, using FPKM-based expression and the pairing information from the series matrix.
- [x] Ensure the FPKM-based expression matrix is suitably transformed for limma (e.g., compute `log2(FPKM + 1)`) and that low-expression genes are filtered according to a chosen rule (e.g., require minimum log2(FPKM+1) in a minimum number of samples).
- [x] Construct a design matrix using `model.matrix()` with the selected design formula and the `colData` (including `condition` and optionally `patient_id`).
- [x] Fit the limma linear model to the log-transformed FPKM data (`lmFit` + `eBayes`) and verify that the coefficient corresponding to the tumor vs normal contrast (e.g., `conditiontumor`) is correctly identified.
- [x] Extract DE results for the tumor vs normal contrast using `topTable()` (log2FC, p-values, adjusted p-values) and summarize the number of significantly up- and down-regulated genes at selected thresholds.
- [x] Apply significance and log2 fold-change cutoffs (e.g., FDR < 0.05 and |log2FC| > 1) and add a `regulation` label (up/down/not_significant) to the results.
- [x] Save the complete DE results table (all genes) to `results/` (e.g., `de_results_limma_fpkm.csv`) for later visualization and downstream analysis.
- [ ] Document key analysis choices (design formula, thresholds, filtering criteria) and explicitly state that FPKM-based limma was used as the DE method because a complete raw counts matrix was not available.

### 8.4 Task 3 – Visualization: Volcano Plot and MA Plot
- [x] Understand volcano plot axes: x = log2 fold change (tumor vs normal), y = -log10(adjusted p-value).
- [x] Understand MA plot axes: x = mean expression (e.g., average log2(FPKM+1)), y = log2 fold change.
- [x] Decide how to represent significance (color by `regulation`, reference lines for `alpha` and `logfc_cutoff`).
- [x] Use the `res_limma_df` table from Task 2 as the single source of DE results.
- [x] Confirm that `logFC`, `adj.P.Val`, `AveExpr`, and `regulation` in `res_limma_df` are defined based on the limma model on log2(FPKM+1) and thresholds `alpha` and `logfc_cutoff`.
- [x] Derive `-log10(adj.P.Val)` for the volcano plot as `neg_log10_adjP`.
- [x] Use `ggplot2` to create a volcano plot with x = `logFC`, y = `-log10(adj.P.Val)`, colored by `regulation`.
- [x] Add vertical lines at `x = ±logfc_cutoff` and a horizontal line at `y = -log10(alpha)` to show significance thresholds consistent with Task 2.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, add labels for a small number of highly significant genes (e.g., using `ggrepel`) if required by the assignment.
- [x] Save the volcano plot to `figures/` with a clear filename (e.g., `volcano_limma_fpkm_tumor_vs_normal.png`).
- [x] Use `ggplot2` to create an MA plot with x = `AveExpr` (average log2 expression) and y = `logFC`, colored by `regulation`.
- [x] Add a horizontal reference line at `y = 0` to indicate no change in expression.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, highlight or label significantly differentially expressed genes if desired.
- [x] Save the MA plot to `figures/` with a descriptive filename (e.g., `ma_limma_fpkm_tumor_vs_normal.png`).
- [ ] Write short, clear narrative text in the report explaining what the volcano plot reveals about the number, direction, and significance of differentially expressed genes.
- [ ] Write short, clear narrative text explaining what the MA plot reveals about how fold changes depend on average expression and whether any systematic patterns or biases are apparent.
- [ ] Explicitly mention that these visualizations are based on limma DE results using log2(FPKM+1) as input rather than DESeq2 on raw counts, due to the lack of a complete counts' matrix.

---

## 5. Task 4 – Downstream Functional and Pathway Analysis

### 5.1 Study Biological / Bioinformatics Background
- Gene sets and pathways:
  - Understand what GO terms (BP/MF/CC), KEGG pathways, Reactome pathways, and MSigDB collections represent.
  - Distinguish between over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
- Enrichment concepts:
  - ORA: uses a list of significant genes vs a background set.
  - GSEA: uses a ranked list of all genes with a continuous statistic.
  - Key outputs: enrichment score, normalized enrichment score (NES), p-value, FDR/q-value.
- Gene identifier types:
  - Study mapping between Ensembl IDs, Entrez IDs, and gene symbols.
  - Learn why consistent ID types are required for each tool.

### 5.2 Prepare Gene Lists for Enrichment
- From DE results, define:
  - List of significantly up-regulated genes (e.g., `padj < threshold` and log2FC > cutoff).
  - List of significantly down-regulated genes (e.g., `padj < threshold` and log2FC < -cutoff).
- Optionally create a ranked gene list:
  - Rank genes by log2FC or a signed statistic (e.g., sign(log2FC) * -log10(pvalue)).

### 5.3 Map Gene Identifiers
- Determine which ID type the enrichment tools expect:
  - DAVID, KEGG, GO tools, or R packages like `clusterProfiler` may require Entrez IDs or gene symbols.
- Use annotation resources:
  - Use `AnnotationDbi` and `org.Hs.eg.db` (if human) to map from gene IDs in the dataset to required IDs.
  - Track how many genes successfully map and how many are lost.

### 5.4 Choose Enrichment Tools and Strategy
- Option A: Web-based tools (e.g., DAVID, Enrichr, web-based GO/KEGG portals):
  - Prepare text files or copy-paste gene lists of up- and down-regulated genes.
  - Select species (likely Homo sapiens).
  - Choose annotation categories (GO BP/MF/CC, KEGG, Reactome, etc.).
  - Download or record top enriched terms (p-value, FDR, gene counts).
- Option B: R-based tools (preferred for reproducibility if allowed/expected):
  - Use `clusterProfiler` for functions like `enrichGO`, `enrichKEGG`, `enricher`, or `GSEA`.
  - Provide gene lists and universe/background as all expressed genes.
  - Optionally use MSigDB gene sets with `fgsea` or similar packages.

### 5.5 Run Enrichment / Pathway Analysis
- Run ORA for up-regulated genes:
  - Identify significantly enriched GO terms and pathways.
- Run ORA for down-regulated genes:
  - Identify enriched terms/pathways potentially associated with suppressed functions in tumors.
- If performing GSEA:
  - Use ranked gene list and chosen gene sets (e.g., Hallmark, KEGG).
  - Extract top positively and negatively enriched pathways.

### 5.6 Summarize and Visualize Enrichment Results
- Visualization:
  - Create bar plots or dot plots of top enriched GO terms/pathways.
  - If using R tools, consider enrichment maps or network plots (optional, depending on assignment scope).
- Interpretation:
  - Identify major biological themes (e.g., cell cycle, immune response, apoptosis, metabolism).
  - Distinguish processes/pathways enriched among up-regulated vs down-regulated genes.
- Documentation:
  - Record thresholds, databases, and tools used for enrichment.
  - Note limitations (e.g., incomplete mapping, redundancy of GO terms, FPKM-based DE).

---

## 6. Integration, Reporting, and Final Checks

### 6.1 Plan the Structure of the Final Report / Answers
- Organize content around the 4 assignment tasks:
  - Task 1: Data loading and understanding.
  - Task 2: Differential expression analysis.
  - Task 3: Visualization (volcano & MA plots).
  - Task 4: Downstream functional and pathway analysis.
- For each task, plan to include:
  - Brief description of methods used.
  - Key results (tables, counts, or summaries).
  - Representative figures if applicable.
  - Short interpretation statements answering the specific question.

### 6.2 Ensure Reproducibility
- Keep all R commands in one or more scripts or an R Markdown file.
- Save intermediate objects (e.g., DE results) if re-running analyses is expensive.
- Set a random seed where needed (e.g., for GSEA implementations that involve randomness).

### 6.3 Final Sanity Checks
- Verify that all 4 assignment questions are explicitly answered using the planned workflow output.
- Cross-check that:
  - Data loaded and preprocessed correctly.
  - DE analysis steps and thresholds match what is described in the text.
  - Plots are consistent with DE results.
  - Enrichment results logically follow from the DE gene lists.
- Note any caveats and limitations clearly (especially the use of FPKM, not raw counts).

---

## 7. Report Writing Guide (5 Main Sections)

This section describes how to structure the written report around the 4 analysis tasks and how to integrate results, figures, and references.

### 7.1 Suggested Report Structure (5 Sections)

1. **Introduction**
   - Brief biological and experimental context:
     - What dataset GSE183947 represents (tumor vs normal, tissue type, number of samples/patients).
     - High-level goal of the analysis (identify differentially expressed genes and pathways between tumor and normal).
   - Motivation:
     - Why differential expression and pathway analysis are useful in this context.
   - Short overview of the workflow:
     - 1–2 sentences summarizing data loading, DE analysis, visualization, and downstream enrichment.

2. **Materials and Methods**
   - **Data description**:
     - Source: `GSE183947_fpkm.csv` (mention GEO if relevant).
     - Variables: what rows (genes) and columns (samples and any metadata) represent.
   - **Preprocessing**:
     - How data were read into R (tools like `readr`, `data.table`).
     - Any filtering of low-expression genes.
     - Handling of missing values or duplicated IDs.
   - **Differential expression analysis**:
     - Statistical model and design formula (e.g., `~ condition` or `~ patient_id + condition`).
     - Software and versions (e.g., R version, DESeq2 version) if used.
     - Thresholds for calling genes differentially expressed (FDR and log2FC cutoffs).
   - **Visualization**:
     - Methods used for volcano and MA plots (e.g., `ggplot2`, `plotMA` from DESeq2).
   - **Downstream/Enrichment analysis**:
     - Tools (e.g., DAVID, clusterProfiler) and gene sets/databases used (GO, KEGG, MSigDB).
     - How gene lists were defined (up/down lists, ranked lists).

3. **Results**
   - **Task 1 – Data understanding**:
     - Number of genes and samples.
     - Summary statistics of expression (e.g., distribution of FPKM, any outliers).
     - Any major QC findings (e.g., clustering of tumor vs normal).  
   - **Task 2 – Differential expression**:
     - Total number of significantly up- and down-regulated genes at chosen thresholds.
     - Maybe a short table of top N DE genes with log2FC and FDR.
   - **Task 3 – Visualization**:
     - Include volcano plot(s) with a brief description of what they show.
     - Include MA plot with commentary on how fold changes vary with expression level.
   - **Task 4 – Downstream analysis**:
     - Tables or bullet lists of top enriched GO terms/pathways for up-regulated and down-regulated genes.
     - 1–2 plots summarizing enrichment (e.g., bar plot of top GO terms).

4. **Discussion**
   - Interpret the main findings:
     - What biological processes/pathways appear activated in tumor vs normal.
     - How these results relate to general cancer biology (e.g., proliferation, apoptosis, immune response).
   - Reflect on consistency between DE results, plots, and enrichment:
     - Do the enriched pathways match the direction and magnitude seen in the DE results and plots?
   - Discuss limitations:
     - Use of FPKM (not raw counts) for DE.
     - Any potential batch effects or unmodeled covariates.
     - Limited sample size, multiple testing issues.

5. **Conclusion and Future Work**
   - Concise summary of the most important results:
     - 2–3 sentences on key DE patterns and enriched pathways.
   - Brief statement on biological implications.
   - Optional: suggestions for follow-up analyses (e.g., validation experiments, additional datasets, deeper pathway exploration).

### 7.2 Integrating Outputs from the 4 Tasks into the Report

- **From Task 1 (Data loading & understanding):**
  - Use key descriptive statistics and any QC plots in the *Materials and Methods* (data description) and *Results* (data overview).
  - Report exact numbers: how many genes, how many samples in each group (tumor vs normal), how pairs were defined.
- **From Task 2 (Differential expression):**
  - Present summary counts of DE genes and top DE genes in the *Results* section.
  - In *Methods*, clearly document the design formula and thresholds used.
  - Use DE tables to support interpretations in the *Discussion* (e.g., highlight important DE genes with known roles in cancer).
- **From Task 3 (Volcano & MA plots):**
  - Include selected figures in the *Results* section with informative captions (what the axes represent, thresholds, and main patterns).
  - Refer to these plots in the text when describing the distribution of DE genes (e.g., “As shown in Figure 1…”).
- **From Task 4 (Enrichment / pathway analysis):**
  - Place tables/plots of enriched pathways in the *Results* section.
  - Use these results heavily in the *Discussion* to build a biological narrative (e.g., “Pathways related to cell cycle and DNA repair were enriched among up-regulated genes, indicating…”).
- Ensure cross-referencing:
  - Number tables and figures (e.g., Table 1, Figure 1) and refer to them by number in the text.
  - Make sure all key outputs from the 4 tasks have a clear place and purpose in the narrative.

### 7.3 Adding Study References to the Report

- **Reference types to include:**
  - Primary publications for the dataset (GSE183947 original study, if available).
  - Methodological references:
    - RNA-Seq and DE methodology papers.
    - DESeq2 or other tools used.
  - Biological background sources:
    - Review articles or textbook chapters describing the relevant cancer type, pathways, or biological processes.
- **Referencing style:**
  - Follow the citation style specified in your course/assignment (e.g., numeric [1], author-year (Smith et al., 2020)).
  - Be consistent throughout the report.
- **Where to cite:**
  - Introduction:
    - Cite sources for biological background and dataset description.
  - Methods:
    - Cite software/tools and methods (e.g., DESeq2, enrichment tools, gene set databases).
  - Discussion:
    - Cite studies that support or contrast your findings.
- **Practical tips:**
  - Keep a small `References` or `Bibliography` section at the end with complete citations.
  - As you work on Tasks 2–4, note down any relevant papers or manuals you consult so they can be cited properly later.

---

## 8. Task Checklist for Progress Tracking

Use this checklist to track your progress through the plan. You can duplicate it into your notes and mark items as you complete them.

### 8.1 Global Preparation
- [x] Install/verify R and required packages (`tidyverse`, `DESeq2`, annotation and enrichment packages).
- [x] Set up project structure (`scripts/`, `results/`, `figures/`).
- [x] Review RNA-Seq, FPKM, and basic DE concepts.
- [x] Understand the assignment requirements and download required data files.

### 8.2 Task 1 – Data Loading and Understanding
- [x] Load `data/GSE183947_fpkm.csv` into R and inspect its structure.
- [x] Load `data/GSE183947_series_matrix.txt` and extract sample metadata (tumor/normal status, pairing if available).
- [x] Load `data/Human.GRCh38.p13.annot.tsv` and verify gene ID compatibility with FPKM symbols using the `Symbol` and `Synonyms` columns.
- [x] Construct the FPKM-based expression matrix (genes in rows, samples in columns) and corresponding `colData` (`sample_id`, `condition`, `patient_id` if possible).
- [x] Perform basic EDA (summaries, distributions, PCA/boxplots) using FPKM expression and sample metadata.
- [x] Document assumptions (e.g., how pairing was inferred, handling of missing values, low-expression filtering, ID mapping decisions, and the decision to not use the partial raw counts file for DE).

### 8.3 Task 2 – Differential Gene Expression Analysis (FPKM + limma)
- [x] Decide on the statistical design for DE (unpaired `~ condition` vs paired `~ patient_id + condition`) based on metadata, using FPKM-based expression and the pairing information from the series matrix.
- [x] Ensure the FPKM-based expression matrix is suitably transformed for limma (e.g., compute `log2(FPKM + 1)`) and that low-expression genes are filtered according to a chosen rule (e.g., require minimum log2(FPKM+1) in a minimum number of samples).
- [x] Construct a design matrix using `model.matrix()` with the selected design formula and the `colData` (including `condition` and optionally `patient_id`).
- [x] Fit the limma linear model to the log-transformed FPKM data (`lmFit` + `eBayes`) and verify that the coefficient corresponding to the tumor vs normal contrast (e.g., `conditiontumor`) is correctly identified.
- [x] Extract DE results for the tumor vs normal contrast using `topTable()` (log2FC, p-values, adjusted p-values) and summarize the number of significantly up- and down-regulated genes at selected thresholds.
- [x] Apply significance and log2 fold-change cutoffs (e.g., FDR < 0.05 and |log2FC| > 1) and add a `regulation` label (up/down/not_significant) to the results.
- [x] Save the complete DE results table (all genes) to `results/` (e.g., `de_results_limma_fpkm.csv`) for later visualization and downstream analysis.
- [ ] Document key analysis choices (design formula, thresholds, filtering criteria) and explicitly state that FPKM-based limma was used as the DE method because a complete raw counts matrix was not available.

### 8.4 Task 3 – Visualization: Volcano Plot and MA Plot
- [x] Understand volcano plot axes: x = log2 fold change (tumor vs normal), y = -log10(adjusted p-value).
- [x] Understand MA plot axes: x = mean expression (e.g., average log2(FPKM+1)), y = log2 fold change.
- [x] Decide how to represent significance (color by `regulation`, reference lines for `alpha` and `logfc_cutoff`).
- [x] Use the `res_limma_df` table from Task 2 as the single source of DE results.
- [x] Confirm that `logFC`, `adj.P.Val`, `AveExpr`, and `regulation` in `res_limma_df` are defined based on the limma model on log2(FPKM+1) and thresholds `alpha` and `logfc_cutoff`.
- [x] Derive `-log10(adj.P.Val)` for the volcano plot as `neg_log10_adjP`.
- [x] Use `ggplot2` to create a volcano plot with x = `logFC`, y = `-log10(adj.P.Val)`, colored by `regulation`.
- [x] Add vertical lines at `x = ±logfc_cutoff` and a horizontal line at `y = -log10(alpha)` to show significance thresholds consistent with Task 2.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, add labels for a small number of highly significant genes (e.g., using `ggrepel`) if required by the assignment.
- [x] Save the volcano plot to `figures/` with a clear filename (e.g., `volcano_limma_fpkm_tumor_vs_normal.png`).
- [x] Use `ggplot2` to create an MA plot with x = `AveExpr` (average log2 expression) and y = `logFC`, colored by `regulation`.
- [x] Add a horizontal reference line at `y = 0` to indicate no change in expression.
- [x] Set a descriptive title and axis labels reflecting the tumor vs normal contrast and the limma on log2(FPKM+1) workflow.
- [x] Optionally, highlight or label significantly differentially expressed genes if desired.
- [x] Save the MA plot to `figures/` with a descriptive filename (e.g., `ma_limma_fpkm_tumor_vs_normal.png`).
- [ ] Write short, clear narrative text in the report explaining what the volcano plot reveals about the number, direction, and significance of differentially expressed genes.
- [ ] Write short, clear narrative text explaining what the MA plot reveals about how fold changes depend on average expression and whether any systematic patterns or biases are apparent.
- [ ] Explicitly mention that these visualizations are based on limma DE results using log2(FPKM+1) as input rather than DESeq2 on raw counts, due to the lack of a complete counts' matrix.

---

## 5. Task 4 – Downstream Functional and Pathway Analysis

### 5.1 Study Biological / Bioinformatics Background
- Gene sets and pathways:
  - Understand what GO terms (BP/MF/CC), KEGG pathways, Reactome pathways, and MSigDB collections represent.
  - Distinguish between over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
- Enrichment concepts:
  - ORA: uses a list of significant genes vs a background set.
  - GSEA: uses a ranked list of all genes with a continuous statistic.
  - Key outputs: enrichment score, normalized enrichment score (NES), p-value, FDR/q-value.
- Gene identifier types:
  - Study mapping between Ensembl IDs, Entrez IDs, and gene symbols.
  - Learn why consistent ID types are required for each tool.

### 5.2 Prepare Gene Lists for Enrichment
- From DE results, define:
  - List of significantly up-regulated genes (e.g., `padj < threshold` and log2FC > cutoff).
  - List of significantly down-regulated genes (e.g., `padj < threshold` and log2FC < -cutoff).
- Optionally create a ranked gene list:
  - Rank genes by log2FC or a signed statistic (e.g., sign(log2FC) * -log10(pvalue)).

### 5.3 Map Gene Identifiers
- Determine which ID type the enrichment tools expect:
  - DAVID, KEGG, GO tools, or R packages like `clusterProfiler` may require Entrez IDs or gene symbols.
- Use annotation resources:
  - Use `AnnotationDbi` and `org.Hs.eg.db` (if human) to map from gene IDs in the dataset to required IDs.
  - Track how many genes successfully map and how many are lost.

### 5.4 Choose Enrichment Tools and Strategy
- Option A: Web-based tools (e.g., DAVID, Enrichr, web-based GO/KEGG portals):
  - Prepare text files or copy-paste gene lists of up- and down-regulated genes.
  - Select species (likely Homo sapiens).
  - Choose annotation categories (GO BP/MF/CC, KEGG, Reactome, etc.).
  - Download or record top enriched terms (p-value, FDR, gene counts).
- Option B: R-based tools (preferred for reproducibility if allowed/expected):
  - Use `clusterProfiler` for functions like `enrichGO`, `enrichKEGG`, `enricher`, or `GSEA`.
  - Provide gene lists and universe/background as all expressed genes.
  - Optionally use MSigDB gene sets with `fgsea` or similar packages.

### 5.5 Run Enrichment / Pathway Analysis
- Run ORA for up-regulated genes:
  - Identify significantly enriched GO terms and pathways.
- Run ORA for down-regulated genes:
  - Identify enriched terms/pathways potentially associated with suppressed functions in tumors.
- If performing GSEA:
  - Use ranked gene list and chosen gene sets (e.g., Hallmark, KEGG).
  - Extract top positively and negatively enriched pathways.

### 5.6 Summarize and Visualize Enrichment Results
- Visualization:
  - Create bar plots or dot plots of top enriched GO terms/pathways.
  - If using R tools, consider enrichment maps or network plots (optional, depending on assignment scope).
- Interpretation:
  - Identify major biological themes (e.g., cell cycle, immune response, apoptosis, metabolism).
  - Distinguish processes/pathways enriched among up-regulated vs down-regulated genes.
- Documentation:
  - Record thresholds, databases, and tools used for enrichment.
  - Note limitations (e.g., incomplete mapping, redundancy of GO terms, FPKM-based DE).

---

## 6. Integration, Reporting, and Final Checks

### 6.1 Plan the Structure of the Final Report / Answers
- Organize content around the 4 assignment tasks:
  - Task 1: Data loading and understanding.
  - Task 2: Differential expression analysis.
  - Task 3: Visualization (volcano & MA plots).
  - Task 4: Downstream functional and pathway analysis.
- For each task, plan to include:
  - Brief description of methods used.
  - Key results (tables, counts, or summaries).
  - Representative figures if applicable.
  - Short interpretation statements answering the specific question.

### 6.2 Ensure Reproducibility
- Keep all R commands in one or more scripts or an R Markdown file.
- Save intermediate objects (e.g., DE results) if re-running analyses is expensive.
- Set a random seed where needed (e.g., for GSEA implementations that involve randomness).

### 6.3 Final Sanity Checks
- Verify that all 4 assignment questions are explicitly answered using the planned workflow output.
- Cross-check that:
  - Data loaded and preprocessed correctly.
  - DE analysis steps and thresholds match what is described in the text.
  - Plots are consistent with DE results.
  - Enrichment results logically follow from the DE gene lists.
- Note any caveats and limitations clearly (especially the use of FPKM, not raw counts).

---

## 7. Report Writing Guide (5 Main Sections)

This section describes how to structure the written report around the 4 analysis tasks and how to integrate results, figures, and references.

### 7.1 Suggested Report Structure (5 Sections)

1. **Introduction**
   - Brief biological and experimental context:
     - What dataset GSE183947 represents (tumor vs normal, tissue type, number of samples/patients).
     - High-level goal of the analysis (identify differentially expressed genes and pathways between tumor and normal).
   - Motivation:
     - Why differential expression and pathway analysis are useful in this context.
   - Short overview of the workflow:
     - 1–2 sentences summarizing data loading, DE analysis, visualization, and downstream enrichment.

2. **Materials and Methods**
   - **Data description**:
     - Source: `GSE183947_fpkm.csv` (mention GEO if relevant).
     - Variables: what rows (genes) and columns (samples and any metadata) represent.
   - **Preprocessing**:
     - How data were read into R (tools like `readr`, `data.table`).
     - Any filtering of low-expression genes.
     - Handling of missing values or duplicated IDs.
   - **Differential expression analysis**:
     - Statistical model and design formula (e.g., `~ condition` or `~ patient_id + condition`).
     - Software and versions (e.g., R version, DESeq2 version) if used.
     - Thresholds for calling genes differentially expressed (FDR and log2FC cutoffs).
   - **Visualization**:
     - Methods used for volcano and MA plots (e.g., `ggplot2`, `plotMA` from DESeq2).
   - **Downstream/Enrichment analysis**:
     - Tools (e.g., DAVID, clusterProfiler) and gene sets/databases used (GO, KEGG, MSigDB).
     - How gene lists were defined (up/down lists, ranked lists).

3. **Results**
   - **Task 1 – Data understanding**:
     - Number of genes and samples.
     - Summary statistics of expression (e.g., distribution of FPKM, any outliers).
     - Any major QC findings (e.g., clustering of tumor vs normal).  
   - **Task 2 – Differential expression**:
     - Total number of significantly up- and down-regulated genes at chosen thresholds.
     - Maybe a short table of top N DE genes with log2FC and FDR.
   - **Task 3 – Visualization**:
     - Include volcano plot(s) with a brief description of what they show.
     - Include MA plot with commentary on how fold changes vary with expression level.
   - **Task 4 – Downstream analysis**:
     - Tables or bullet lists of top enriched GO terms/pathways for up-regulated and down-regulated genes.
     - 1–2 plots summarizing enrichment (e.g., bar plot of top GO terms).

4. **Discussion**
   - Interpret the main findings:
     - What biological processes/pathways appear activated in tumor vs normal.
     - How these results relate to general cancer biology (e.g., proliferation, apoptosis, immune response).
   - Reflect on consistency between DE results, plots, and enrichment:
     - Do the enriched pathways match the direction and magnitude seen in the DE results and plots?
   - Discuss limitations:
     - Use of FPKM (not raw counts) for DE.
     - Any potential batch effects or unmodeled covariates.
     - Limited sample size, multiple testing issues.

5. **Conclusion and Future Work**
   - Concise summary of the most important results:
     - 2–3 sentences on key DE patterns and enriched pathways.
   - Brief statement on biological implications.
   - Optional: suggestions for follow-up analyses (e.g., validation experiments, additional datasets, deeper pathway exploration).

### 7.2 Integrating Outputs from the 4 Tasks into the Report

- **From Task 1 (Data loading & understanding):**
  - Use key descriptive statistics and any QC plots in the *Materials and Methods* (data description) and *Results* (data overview).
  - Report exact numbers: how many genes, how many samples in each group (tumor vs normal), how pairs were defined.
- **From Task 2 (Differential expression):**
  - Present summary counts of DE genes and top DE genes in the *Results* section.
  - In *Methods*, clearly document the design formula and thresholds used.
  - Use DE tables to support interpretations in the *Discussion* (e.g., highlight important DE genes with known roles in cancer).
- **From Task 3 (Volcano & MA plots):**
  - Include selected figures in the *Results* section with informative captions (what the axes represent, thresholds, and main patterns).
  - Refer to these plots in the text when describing the distribution of DE genes (e.g., “As shown in Figure 1…”).
- **From Task 4 (Enrichment / pathway analysis):**
  - Place tables/plots of enriched pathways in the *Results* section.
  - Use these results heavily in the *Discussion* to build a biological narrative (e.g., “Pathways related to cell cycle and DNA repair were enriched among up-regulated genes, indicating…”).
- Ensure cross-referencing:
  - Number tables and figures (e.g., Table 1, Figure 1) and refer to them by number in the text.
  - Make sure all key outputs from the 4 tasks have a clear place and purpose in the narrative.

### 7.3 Adding Study References to the Report

- **Reference types to include:**
  - Primary publications for the dataset (GSE183947 original study, if available).
  - Methodological references:
    - RNA-Seq and DE methodology papers.
    - DESeq2 or other tools used.
  - Biological background sources:
    - Review articles or textbook chapters describing the relevant cancer type, pathways, or biological processes.
- **Referencing style:**
  - Follow the citation style specified in your course/assignment (e.g., numeric [1], author-year (Smith et al., 2020)).
  - Be consistent throughout the report.
- **Where to cite:**
  - Introduction:
    - Cite sources for biological background and dataset description.
  - Methods:
    - Cite software/tools and methods (e.g., DESeq2, enrichment tools, gene set databases).
  - Discussion:
    - Cite studies that support or contrast your findings.
- **Practical tips:**
  - Keep a small `References` or `Bibliography` section at the end with complete citations.
  - As you work on Tasks 2–4, note down any relevant papers or manuals you consult so they can be cited properly later.

---
