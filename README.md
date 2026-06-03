# 🧩 dialipar

Process DIA-NN and Spectronaut LiP‑MS data with ease and clarity.

Welcome to dialipar — a toolkit for analyzing LiP‑MS experiments. Built around the DIA‑LiPA workflow, dialipar helps you process LiP‑MS data searched with DIA‑NN (semi‑tryptic libraries) or Spectronaut, and turn them into clear, high‑quality results.

With a few simple steps you can:

- Clean and preprocess your DIA LiP‑MS data
- Run differential analysis to identify proteolytic changes
- Visualize results with volcano plots, heatmaps, and QC metrics
- Generate polished, shareable reports

Whether you’re new to LiP‑MS or a seasoned proteomics researcher, dialipar makes your workflow smoother, faster, and more reproducible.

DiaLipa data analysis is described in details in  [C. Van Leene at. al.](https://pubs.acs.org/doi/full/10.1021/acs.analchem.5c07014)

---

## 🛠️ Requirements

1. **R**>= 4.4.0
2. **Install Bioconductor**:
   - install.packages("BiocManager")
3. **Install devtools**:
   - install.packages("devtools")
4. **Install Quarto** (version 1.6.43 or more):
   - Follow the instructions at https://quarto.org/docs/download/
5. **Install Bioconductor dependencies**   
    ``` r 
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
 
        BiocManager::install(c("QFeatures", "SummarizedExperiment", "MsCoreUtils", "msqrob2","scater"))
   ```
6. **Install CRAN dependencies**   
    ``` r 
      install.packages(c("arrow", "dplyr", "stringr","seqinr","withr",
                     "tidyr",,"tibble","quarto","logger",'fs',
                     'knitr','asserthat','magrittr',"yaml","plotly",
                     "ggrepel","ggplot2","ggsci"))
   ```

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

The package contains Quarto template for generating  the following reports:

- `Standard` — Template for DIA‑LiPA analysis. This includes QC plots, unsupervised analsys (MDS plot ), and differential analysis at precursor level for the comparisons of interest

Specify the template name via the `template` parameter in the `render_dialipa_report` function.

---

##  📦 Data Availability

LipMS experiment processed with DIA-NN and Spectronaut can be found in [MassIVE](<ftp://MSV000099740@massive-ftp.ucsd.edu/>). However, The raw files were searched with the laster version opf the software are available also here for download :

-  [DIA-NN v.2.3 library free]  [LF_DIANN_v2-3](https://cloud.cmb.ugent.be/index.php/s/eRDy32eFm7GsRY2) (predicted spectral library). TC and LiP runs are included in only one parquet file. 

-  [DIA-NN v.2.3 empirical library]  [empirical_DIANN_v2-3](https://cloud.cmb.ugent.be/index.php/s/eRDy32eFm7GsRY2) with empirical DDA based spectral library. TC and Lip runs are saved intwo separate parquest files.
-  [Spectronaut v.20.3 empirical spectral library ]  [empirical_spectronaut_v20-3](https://cloud.cmb.ugent.be/index.php/s/eRDy32eFm7GsRY2) with empirical DDA based spectral library. TC and Lip runs are saved intwo separate parquest files.
-  [Spectronaut v.20.3  library free]  [LF_spectronaut_v20-3](https://cloud.cmb.ugent.be/index.php/s/eRDy32eFm7GsRY2) (predicted spectral library). TC and LiP runs are included in only one parquet file.

In the each folder, you can also find the relative experiment description and the Human fasta files.



## 🚀 How to run an analysis

1. Download the required data from Pride 

- parquet files 
- experiment design file 
- fasta file 

See section [Data Availability](#Data-Availability) for details.

2. Put all the file into a folder, (e.g.`../path/dialipar_test/`).
3. Use the code below (adjust paths as needed)
✅ **Always use full paths for files and folders.**

```r
report_target_folder <- "../path/dialipar_test/DIA-LiPA_result"
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
params_start$paired <- FALSE


# Run the report rendering
render_dialipa_report(
  params = params_start,
  template = 'Standard',
  report_folder = report_target_folder,
  report_filename = output_filename
)
```
💡*Note for Windows users:*: When specifying file paths, use double backslashes `(\\)` or forward slashes `(/)`.
```r 
params_start$design_file <- "C:\\Users\\Name\\dialipar_test\\SampleAnnotation.txt"
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
- **paired**: Boolean flag to select the correction for paired (TRUE) or unpaired (FALSE) experiment


In case both lip and Tc runs are in just one parquet file, *input_file_tc* can be not specified. 
---

## 📝 Experiment Design File (EDF)

The experiment design file is a tab‑separated text file and must include the following columns:

- Run: Raw file name without the file extension (mzML/.d/.raw)
- Pipeline: "TC" for trypsin runs and "LiP" for LiP (semi‑tryptic) runs
- Treatment: Drug used in the experiment (e.g., Rapa/DMSO, used in the model)
- Condition: Experimental groups 
- Replicate: Replicate identifier
- CondRep: Combined label consisting of condition and replicate (e.g., DMSO_LiP_1)

Example rows:

| Run                                                                 | Pipeline | Treatment | Condition   | Replicate | CondRep         |
|----------------------------------------------------------------------|----------|------|-------------|-----------|------------------|
| F017128_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_1            | LiP      | DMSO | DMSO_LiP    | 1         | DMSO_LiP_1       |
| F017130_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_2            | LiP      | Rapa | Rapa_LiP    | 1         | Rapa_LiP_1       |
| F017144_1p_Ih19um_trapPM3_Neo__CMB-1934__Chloe_lipMSDIA_9            | TC       | DMSO | DMSO_TC     | 1         | DMSO_TC_1        |


---
