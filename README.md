# 🧩 dialipar

Process DIA-NN and Spectronaut LiP‑MS data with ease and clarity.

Welcome to dialipar — a toolkit for analyzing LiP‑MS experiments. Built around the DIA‑LiPA workflow, dialipar helps you process LiP‑MS data searched with DIA‑NN (semi‑tryptic libraries) or Spectronaut, and turn them into clear, high‑quality results.

With a few simple steps you can:

- Clean and preprocess your LiP‑MS data
- Run differential analysis to identify proteolytic changes
- Visualize results with volcano plots, heatmaps, and QC metrics
- Generate polished, shareable reports

Whether you’re new to LiP‑MS or a seasoned proteomics researcher, dialipar makes your workflow smoother, faster, and more reproducible.

---

## 🛠️ Requirements

1. R >= 4.4.0
2. Install Bioconductor:
   - install.packages("BiocManager")
3. Install PhantomJS for web screenshots (used by webshot):
   - webshot::install_phantomjs()
4. Install devtools:
   - install.packages("devtools")
5. Install Quarto (version 1.6.43 or more):
   - Follow the instructions at https://quarto.org/docs/download/

---

## 🛠 How to install the R package

For development:

1. Clone the GitHub repository.
2. Install the package locally using devtools:

```r
devtools::install()
```

To install the package from GitHub for testing:

```r
devtools::install_github("Gevaert-Lab/dialipar")
```

---

## 📂 Quarto templates included

The repository contains Quarto templates for generating reports:

- `Template_.qmd` — Template for standard DIA‑LiPA analysis

Specify the template name via the `template_file` parameter in the `render_dialipa_report` function.

---

## 🚀 How to run an analysis

1. Download DATA_TO_BE_DECIDED.zip (this should include the LiP‑MS experiment files and the experiment design file).
2. Unzip it into a folder, for example `../path/dialipar_test/`.
3. Use the code below (adjust paths as needed). Always use full paths for files and folders.

```r
report_target_folder <- "../path/dialipar_test/DIA-LiPA_result"
template_file <- "Template_.qmd"
output_filename <- "DIA-LiPA_report.html"

params_start <- list()
params_start$description <- "DIA-LiPA from DIA-NN input (TC 1 LiP in the same parquet file)"
params_start$design_file <- "../path/dialipar_test/SampleAnnotation.txt"
params_start$input_file_lip <- "../path/dialipar_test/report.parquet"
params_start$input_file_tc <- ""
params_start$fasta_file <- "../path/dialipar_test/SP_9606_PK.fasta"
params_start$folder_prj <- report_target_folder

# Report metadata
params_start$title <- "DIA-LiPA Report"
params_start$subtitle <- "DIA-LiPA"
params_start$author <- "Your Name"

# Analysis parameters
params_start$formula <- "~ -1 + Condition"
params_start$comparisons <- c("ConditionRapa_LiP - ConditionDMSO_LiP")
params_start$FC_thr <- 1
params_start$adjPval_thr <- 0.1
params_start$comparison_label <- c("Rapa - DMSO")

# Run the report rendering
render_dialipa_report(
  params = params_start,
  template = template_file,
  report_folder = report_target_folder,
  report_filename = output_filename
)
```

---

## ⚙️ Input parameters

Parameters available to customize the DIA‑LiPA analysis and the generated HTML report:

- **title**: Report title
- **subtitle**: Report subtitle
- **author**: Author name
- **description**: Description of the experiment
- **input_file_lip**: Path to the DIA‑NN or Spectronaut report file (parquet)
- **input_file_tc**: Path to the TC (trypsin control) report file (parquet), if separate
- **design_file**: Path to the experiment design file (tab‑separated)
- **fasta_file**: Path to the fasta file
- **folder_prj**: Path to the root output folder for the analysis
- **formula**: Formula used in MSqrob for the linear model (e.g., "~ -1 + Condition")
- **contrast**: Name of the column in the experiment design file used in the model (e.g., "Condition")
- **FC_thr**: log2 fold‑change threshold (default: 1)
- **adjPval_thr**: Adjusted p‑value threshold for significance (default: 0.05)
- **POI**: Proteins of interest, specified using UniProt IDs. If not, POIs will include the proteins with statistal significant precursors.
- **comparisons**: List of comparisons to test, including the variable name specified in the model (e.g., "ConditionA - ConditionB")
- **comparison_label**: Human‑readable labels for comparisons (e.g., "GroupA - GroupB" → "A - B")

In case both lip anff Tc runs are in just one parquet file, *input_file_tc* can be not specified. 
---

## 📝 Experiment Design File (EDF)

The experiment design file is a tab‑separated text file and must include the following columns:

- Run: Raw file name without the file extension (mzML/.d/.raw)
- Pipeline: "TC" for trypsin runs and "LiP" for LiP (semi‑tryptic) runs
- Drug: Drug used in the experiment (e.g., Rapa/DMSO)
- Condition: Experimental groups (used in the model)
- Replicate: Replicate identifier
- CondRep: Combined label consisting of condition and replicate (e.g., DMSO_LiP_1)

Example rows:

| Run                                                                 | Pipeline | Drug | Condition   | Replicate | CondRep         |
|----------------------------------------------------------------------|----------|------|-------------|-----------|------------------|
| F017128_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_1            | LiP      | DMSO | DMSO_LiP    | 1         | DMSO_LiP_1       |
| F017130_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_2            | LiP      | Rapa | Rapa_LiP    | 1         | Rapa_LiP_1       |
| F017144_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_9            | TC       | DMSO | DMSO_TC     | 1         | DMSO_TC_1        |


---
