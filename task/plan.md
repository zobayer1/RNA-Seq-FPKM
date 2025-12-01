# RNA-Seq FPKM Differential Expression & Downstream Analysis Plan

This plan outlines a detailed, step-by-step workflow for solving all 4 problems in the assignment using R and the datasets in the `data/` directory:
- `GSE183947_raw_counts_GRCh38.p13_NCBI.tsv` – raw read counts for DESeq2.
- `GSE183947_fpkm.csv` – FPKM-normalized expression for exploration/visualization.
- `GSE183947_series_matrix.txt` – GEO series matrix with sample metadata (e.g., tumor vs normal, pairing).
- `Human.GRCh38.p13.annot.tsv` – gene annotation for mapping Ensembl IDs to symbols and other attributes.

It focuses on what needs to be studied and done, not on the final solutions.

---

## 1. Global Preparation and Background

### 1.1 Understand the Assignment and Data
- Read `task.md` carefully to identify the 4 specific questions/problems and any constraints (e.g., tools, packages, or plots required).
- Confirm how each file in `data/` is used:
  - `GSE183947_raw_counts_GRCh38.p13_NCBI.tsv` → primary input for DESeq2-based differential expression.
  - `GSE183947_fpkm.csv` → normalized expression matrix used mainly for exploratory plots / QC.
  - `GSE183947_series_matrix.txt` → sample-level metadata, including which samples are tumor vs normal and how they are paired.
  - `Human.GRCh38.p13.annot.tsv` → mapping from gene IDs (e.g., Ensembl) to human-readable gene symbols and possibly Entrez IDs.
- Note that the study involves paired samples (tumor vs normal tissue from the same patients) and plan to use that information in the design when possible.
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

### 2.3 Briefly Inspect Raw Counts TSV
- Load `data/GSE183947_raw_counts_GRCh38.p13_NCBI.tsv` (at least a subset for inspection).
- Examine structure:
  - Confirm that rows correspond to the same gene identifiers used in the FPKM file.
  - Confirm that columns correspond to samples, and compare sample names with FPKM.
- Document intended use:
  - Note that this file will be the primary input for DESeq2 in Task 2.
  - Record any differences in sample naming or gene IDs relative to the FPKM file.

### 2.4 Build Expression Matrix and Sample Metadata from FPKM and Series Matrix
- Construct the expression matrix from FPKM:
  - Set row names to gene identifiers.
  - Ensure all expression columns are numeric.
  - Confirm that rows = genes and columns = samples.
- Load sample metadata from GEO series matrix:
  - Read `data/GSE183947_series_matrix.txt`.
  - Extract sample IDs and relevant characteristics (e.g., tumor vs normal status, any available pairing/patient information).
- Create sample metadata (`colData`):
  - Derive a `sample_id` column matching the column names of the FPKM and counts matrices.
  - Create a `condition` factor with levels `normal` and `tumor` based on the series matrix.
  - If possible, create a `patient_id` to represent paired samples (same patient in normal and tumor).
- Verify consistency:
  - Confirm that the order of samples in `colData` matches the columns in the FPKM expression matrix.
  - Check that sample IDs in the series matrix, FPKM CSV, and counts TSV are consistent or can be reconciled.

### 2.5 Load and Check Gene Annotation
- Load `data/Human.GRCh38.p13.annot.tsv` into R.
- Identify key columns:
  - Gene ID column (e.g., Ensembl ID) expected to match IDs in FPKM/counts.
  - Annotation fields such as gene symbol, gene name/description, and biotype.
- Check compatibility:
  - Verify that most gene IDs in the FPKM/counts files appear in the annotation table.
  - Note any systematic differences in ID formats (e.g., versioned vs unversioned Ensembl IDs).
- Plan use:
  - Decide which annotation columns will be merged into DE results later (e.g., gene symbol, biotype).
  - Document that full annotation merging will occur in Task 2/Task 4, while Task 1 focuses on understanding and basic checks.

### 2.6 Basic Exploratory Data Analysis (EDA)
- Compute simple summaries using the FPKM expression matrix:
  - Number of genes and samples.
  - Distribution of FPKM values for a few random samples.
- Optional quick visualizations for understanding:
  - Boxplots of log-transformed FPKM per sample to assess overall distribution.
  - Sample-sample correlation heatmap or PCA plot for global patterns (e.g., do tumors cluster separately from normals?).
- Cross-check metadata:
  - Summarize the number of `normal` vs `tumor` samples.
  - If `patient_id` is available, verify that each patient has both a normal and a tumor sample.

### 2.7 Document Assumptions and Limitations
- Clearly note:
  - The nature of the data (FPKM vs raw counts) and that DESeq2 will use raw counts from the TSV file in Task 2.
  - Any preprocessing steps applied at this stage (filtering low-expression genes, handling missing values, resolving duplicate IDs).
  - Any assumptions about pairing inferred from sample names or metadata in the series matrix.
  - Any ID-format issues discovered between expression files and annotation and how you plan to handle them.

---

## 3. Task 2 – Differential Gene Expression Analysis

### 3.1 Study Statistical and Methodological Foundations
- DESeq2 fundamentals:
  - Understand that DESeq2 is designed for count data (negative binomial model), not FPKM.
  - If the assignment expects DESeq2, recognize and plan to discuss this limitation.
  - Learn the main steps: size factor estimation, dispersion estimation, model fitting, and hypothesis testing.
- Design formulas:
  - Simple design: `~ condition` (unpaired comparison of tumor vs normal).
  - Paired design: `~ patient_id + condition` to account for matched pairs.
  - Understand how DESeq2 interprets factors and how the reference level of `condition` determines direction of log2FC.
- Multiple testing and thresholds:
  - Review Benjamini–Hochberg FDR control.
  - Decide on significance thresholds (e.g., `padj < 0.05` or more stringent).
  - Decide on a biologically meaningful log2FC cutoff (e.g., |log2FC| > 1).

### 3.2 Prepare Data for DE Analysis
- Ensure the expression matrix is in a suitable format:
  - If raw counts are unavailable and only FPKM are provided, decide how the assignment wants to proceed (e.g., treat FPKM approximately like counts, or apply a transformation) and document this clearly.
- Finalize `colData`:
  - Confirm `condition` is a factor and choose the reference level (e.g., `normal` as reference so positive log2FC means higher in tumor).
  - Include `patient_id` if a paired analysis is required or recommended.
- Filter genes (optional but common):
  - Remove genes with extremely low expression across all samples (e.g., FPKM ~0), as they may not be informative for DE.

### 3.3 Run Differential Expression
- Construct a `DESeqDataSet` using `DESeqDataSetFromMatrix()` (if DESeq2 is used):
  - Provide expression matrix, `colData`, and an appropriate design formula.
- Run the DESeq2 workflow:
  - Call `DESeq()` to perform estimation and fitting.
  - Extract results with `results()` specifying the contrast (tumor vs normal).
- Inspect overall results:
  - Use `summary()` to see the count of significant up- and down-regulated genes.
  - Look at distributions of `log2FoldChange`, `pvalue`, and `padj`.

### 3.4 Post-Processing and Annotation of DE Results
- Apply significance and fold-change thresholds:
  - Create logical flags for significant genes (`padj < threshold`).
  - Create a `regulation` label: `up`, `down`, or `not_significant` based on sign and magnitude of log2FC.
- Annotate genes:
  - If gene identifiers are not human-readable (e.g., Ensembl IDs), map them to gene symbols (e.g., using `org.Hs.eg.db`).
  - Optionally add gene descriptions if available.
- Save DE results:
  - Export a complete results table (including all genes) to `results/` as a CSV.
  - Optionally, export a filtered table of only significant genes.

### 3.5 Quality Assessment of DE Analysis
- Check diagnostic outputs:
  - Plot p-value histogram to ensure it looks reasonable (e.g., uniform under null with a peak near 0 for true positives).
  - Review the number of DE genes; too many or too few may indicate issues.
- Document methodological choices:
  - Design formula used and rationale (paired vs unpaired).
  - Threshold values for `padj` and log2FC.
  - Any filters applied before DE analysis.

---

## 4. Task 3 – Visualization: Volcano Plot and MA Plot

### 4.1 Study Visualization Principles
- Volcano plot:
  - x-axis: log2 fold change (tumor vs normal).
  - y-axis: -log10(p-value) or -log10(adjusted p-value).
  - How to represent significance:
    - Color points by significance and direction (up/down).
    - Add horizontal/vertical lines indicating p-value and log2FC thresholds.
- MA plot:
  - x-axis: mean expression (e.g., log10 of baseMean or average FPKM).
  - y-axis: log2 fold change.
  - Purpose: to visualize how fold changes depend on expression level.
  - Understand built-in DESeq2 `plotMA()` behavior vs custom plots.

### 4.2 Prepare Data for Plots
- From the DE results table, derive plotting variables:
  - `log2FoldChange` for x-axis (volcano) or y-axis (MA plot).
  - `padj` or `pvalue` to compute `-log10(p)`.
  - `baseMean` or mean expression per gene for MA plot x-axis.
  - `regulation` label used for coloring points (up/down/not significant).
- Select thresholds for highlighting genes:
  - Decide which genes to highlight/label (e.g., top 10 by smallest `padj`).

### 4.3 Create Volcano Plot in R
- Use `ggplot2` to create a scatter plot:
  - Points: each gene as a point.
  - Color: by `regulation`.
  - Add reference lines:
    - Vertical lines for chosen log2FC cutoffs.
    - Horizontal line for chosen p-value/FDR cutoff.
- Optionally add labels:
  - Label a small number of highly significant genes (e.g., using `ggrepel`).
- Save the plot:
  - Export to `figures/` with a clear filename (e.g., `volcano_tumor_vs_normal.png`).

### 4.4 Create MA Plot in R
- Option 1: Use DESeq2 built-in `plotMA()` (if DESeq2 is used):
  - Provide the DESeq2 results object.
  - Customize aesthetics if needed.
- Option 2: Custom MA plot with `ggplot2`:
  - x-axis: log-transformed mean expression or `baseMean`.
  - y-axis: `log2FoldChange`.
  - Color points by `regulation`.
  - Optionally overlay smoothing or highlight significant genes.
- Save the plot:
  - Export to `figures/` with descriptive filename (e.g., `ma_plot_tumor_vs_normal.png`).

### 4.5 Interpret and Document Visualizations
- Write short, clear descriptions for the report:
  - What the volcano plot shows about the number and direction of DE genes.
  - What the MA plot reveals about the dependence of fold changes on expression levels.
  - Any notable genes or patterns highlighted in these plots.

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
- [x] Load `data/GSE183947_raw_counts_GRCh38.p13_NCBI.tsv` and inspect its structure.
- [x] Load `data/GSE183947_series_matrix.txt` and extract sample metadata (tumor/normal status, pairing if available).
- [x] Load `data/Human.GRCh38.p13.annot.tsv` and verify gene ID compatibility with expression data.
- [x] Construct FPKM-based expression matrix (genes in rows, samples in columns) and corresponding `colData` (`sample_id`, `condition`, `patient_id` if possible).
- [x] Perform basic EDA (summaries, distributions, optional QC plots) using FPKM expression and sample metadata.
- [x] Document assumptions (e.g., how pairing was inferred, handling of missing values, low-expression filtering, ID mapping decisions).

### 8.3 Task 2 – Differential Gene Expression Analysis
- [ ] Decide on the statistical design for DE (unpaired `~ condition` vs paired `~ patient_id + condition`) based on metadata.
- [ ] Load the full raw counts matrix from `data/GSE183947_raw_counts_GRCh38.p13_NCBI.tsv` into R.
- [ ] Ensure raw counts columns align with `colData$sample_id` (same samples, same order) and reconcile any naming differences.
- [ ] Filter out genes with extremely low counts across all samples (e.g., near-zero expression) according to a chosen rule.
- [ ] Construct a `DESeqDataSet` (or chosen DE object) using the raw counts matrix and `colData`, with the selected design formula.
- [ ] Run the DE pipeline (e.g., `DESeq()`), obtain results for tumor vs normal, and summarize the number of up- and down-regulated genes at selected thresholds.
- [ ] Apply significance and log2 fold-change cutoffs and add a `regulation` label (up/down/not_significant) to the results.
- [ ] Annotate DE results with gene-level information from `Human.GRCh38.p13.annot.tsv` (e.g., gene symbols, biotypes).
- [ ] Export complete and filtered DE results tables to `results/` for later visualization and downstream analysis.
- [ ] Document key analysis choices (design formula, thresholds, filtering criteria) and any notable QC observations from the DE step.

### 8.4 Task 3 – Visualization: Volcano Plot and MA Plot
- [ ] From the DE results table, derive plotting variables (log2FC, p-values/FDR, baseMean/mean expression, `regulation` labels).
- [ ] Implement a volcano plot (using `ggplot2` or equivalent) with log2FC on the x-axis and -log10(p or padj) on the y-axis, coloring points by `regulation`.
- [ ] Add threshold reference lines to the volcano plot (log2FC cutoffs and p-value/FDR cutoff) and optionally highlight/label top genes.
- [ ] Implement an MA plot (either via DESeq2 `plotMA()` or a custom `ggplot2` version) showing mean expression vs log2FC, colored by `regulation`.
- [ ] Save the volcano and MA plots to `figures/` with clear filenames and sufficient resolution for inclusion in the report.
- [ ] Write brief interpretations of what the volcano and MA plots show about the overall DE pattern (for use in the Results/Discussion sections).

### 8.5 Task 4 – Downstream Functional and Pathway Analysis
- [ ] Define gene lists for enrichment: significantly up-regulated and down-regulated genes based on chosen thresholds.
- [ ] Optionally create a ranked gene list (e.g., by log2FC or a signed statistic) if performing GSEA-type analyses.
- [ ] Map DE gene identifiers to the ID type required by chosen enrichment tools (e.g., Entrez ID or gene symbols) using `Human.GRCh38.p13.annot.tsv` and/or `org.Hs.eg.db`.
- [ ] Choose one or more enrichment tools/workflows (e.g., clusterProfiler, DAVID, Enrichr, or web-based GO/KEGG portals) and configure parameters (background universe, databases, cutoffs).
- [ ] Run enrichment analyses separately for up-regulated and down-regulated genes (and optionally GSEA on a ranked list).
- [ ] Summarize top enriched GO terms and pathways (e.g., KEGG, Reactome) with their p-values/FDR and gene counts.
- [ ] Create at least one visualization of enrichment results (e.g., bar plot or dot plot of top terms/pathways) and save to `figures/`.
- [ ] Interpret major biological themes emerging from enrichment (e.g., cell cycle, immune response, metabolism) and relate them to tumor vs normal biology.
- [ ] Document enrichment settings (databases, thresholds, background) and limitations (e.g., incomplete mapping, redundancy of GO terms).

### 8.6 Task 5 – Summary and Report Organization
- [ ] Draft a structured outline for the final report organized by the four main tasks (data loading, DE, visualization, enrichment), following the suggested 5-section format.
- [ ] For each task, list the specific figures, tables, and key numerical results that must appear in the report (with tentative figure/table labels).
- [ ] Collect and record the main conclusions from each task (e.g., number of DE genes, key pathways, QC findings) in concise bullet points.
- [ ] Ensure that all important analysis choices (design formula, thresholds, filters, tools used) are noted for inclusion in the Materials and Methods section.
- [ ] Identify any caveats or limitations (e.g., FPKM vs counts, sample size, unmodeled covariates) that should be discussed explicitly.
- [ ] Compile a preliminary list of references to cite (dataset publication, DESeq2, enrichment tools, key biological background papers).
- [ ] Verify that all assignment questions are explicitly answered by mapping each question to the corresponding results, figures, and discussion points in the report outline.
