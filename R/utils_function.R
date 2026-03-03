
#' @title render_quarto_template
#' 
#' @description
#' This internal function manages the heavy lifting of the reporting layer. It 
#' creates a temporary sandbox, bundles all analysis results into a single 
#' RDS "data bag," updates the Quarto parameters, and renders the report. 
#' Finally, it cleans up the temporary directory and moves the output to 
#' the designated folder.
#'
#' @param data_list A named list containing all data objects (tables, stats) 
#' required by the Quarto template.
#' @param template_name Character string. The name of the `.qmd` file located 
#' within the package's `quarto_template` directory.
#' @param report_fld Character string. The final destination directory for 
#' the rendered report and its associated files.
#' @param report_fname Character string. The name of the resulting HTML file 
#' (e.g., "Nterm_Report.html").
#' @param params_report A list of parameters to be passed to Quarto. This 
#' function automatically appends a `data_path` element to this list.
#'
#' @return The function returns `NULL` invisibly. Its primary purpose is 
#' the side effect of file creation and directory management.
#' 
#' @keywords internal
#' 
#' @importFrom quarto quarto_render
#' @importFrom logger log_info log_error
#' @importFrom withr with_dir
#' @importFrom tools file_path_sans_ext
render_quarto_template <- function(data_list, template_name, report_fld, report_fname, params_report) {
  # 1. Setup Temp Dir (Your current logic is good here)
  template_source_folder <- system.file("quarto_template", package = "dialipar")
  if (template_source_folder == "") {
    stop("Template folder not found in the package.")
  }

  temp_work_dir <- file.path(tempdir(), paste0("quarto_temp_", Sys.getpid()))
  dir.create(temp_work_dir, recursive = TRUE, showWarnings = FALSE)
  log_info('Temp folder created : {temp_work_dir}')
  success <- file.copy(from = template_source_folder,
                       to = temp_work_dir,
                       recursive = TRUE)
   if (!success) {
    stop("Failed to copy the template folder to the temporary directory.")
  }
  log_info('Copy Template file ...done')
  
  # 2. THE BIG CHANGE: Save everything into ONE file
  data_rds_path <- file.path(temp_work_dir, "report_data.rds")
  saveRDS(data_list, data_rds_path)
  
  log_info('Copy Rds Results ...done')
  # 3. Add that path to parameters
  params_report$data_path <- data_rds_path

    # Construct the path to the copied template file in the temp directory.
  # Assumes that the template file is directly inside the copied folder.
  temp_template_path <- file.path(temp_work_dir, basename(template_source_folder), template_name)
  if (!file.exists(temp_template_path)) {
    stop("Template file not found in the temporary directory: ", temp_template_path)
  }

  path <- file.path(temp_work_dir, basename(template_source_folder))
  
    tryCatch({
    withr::with_dir(path, {
      quarto_render(
        input = temp_template_path,
        output_format = "html",
        output_file = report_fname,
        execute_params = params_report,
        quarto_args = c( "--no-clean", "--output-dir", path)
      )
    })
  }, error = function(e) {
    log_error("Error in Quarto rendering: {e$message}")
    unlink(temp_work_dir, recursive = TRUE)
    stop(e)
  })

  resource_folder_name <- paste0(tools::file_path_sans_ext(report_fname), "_files")
  rendered_report_path <- file.path(path, report_fname)

  if (!dir.exists(report_fld)) {
    dir.create(report_fld, recursive = TRUE)
  }
    log_info('Copying rendered html report ...')
  # Copy the rendered HTML report to the target folder
  file.copy(from = rendered_report_path, to = file.path(report_fld, report_fname), overwrite = TRUE)
  # If a resource folder was generated, copy it as well
  temp_resource_path <- file.path(temp_work_dir, resource_folder_name)
  if (dir.exists(temp_resource_path)) {
    file.copy(from = temp_resource_path,
              to = file.path(report_fld, resource_folder_name),
              recursive = TRUE, overwrite = TRUE)
  }
  log_info('Cleaning temp folder ...')
  # Optionally, remove the temporary working directory to clean up
  unlink(temp_work_dir, recursive = TRUE)
  return(invisible(NULL))
  
}


process_dialipa_data <- function (params_report, analysis_type = "unpaired" ){
      
  fastaproc <- read_fasta_ann(params_report$fasta_file )
      input_data <- parse_input(params_report$input_file_tc, 
                          params_report$input_file_lip,  
                          dual = FALSE, 
                          params_report$design_file)

        qf_base <- create_qfeat_(input_data$design, input_data$lip, input_data$tc)
        if (qf_base$status == 1) stop(qf_base$error)

        qf_norm <- normalization_scaling_factor(qf_base$result)
        if (qf_norm$status == 1) stop(qf_norm$error) 

        qf_final <- qf_norm$result 

        if (analysis_type == 'paired') {
            usage_res <- compute_usage(qf_norm$result) 
            if (usage_res$status == 1) stop(usage_res$error)
            qf_final <- usage_res$result
        }
        qc_ann <- qc_precursor_annotation(qf_final,fastaproc$result, type= analysis_type  )
        if (qc_ann$status == 1) stop(qc_ann$error)
        #debug 
        #saveRDS(qc_ann$result, './TEST_exp.Rds')
        #stop('Debug')
        ## qc_ann$result go to data bag.
        if (analysis_type == 'paired') {
              log_info('Paired branch ...')
              res_de_norm <-  msqrob_model(pe = qf_final, params = params_report, layer = 'precursors_lip_norm' )
              if (res_de_norm$status == 1) stop(res_de_norm$error)
              res_de_usage <-  msqrob_model(pe = res_de_norm$q_feat, params = params_report, layer = 'precursors_lip_usage' )
              if (res_de_usage$status == 1) stop(res_de_usage$error)
              df_ann <-  as.data.frame(rowData(res_de_usage$q_feat[["precursors_lip_norm"]])[c("Precursor.Id", "Protein.Group", "Genes", "Proteotypic", "Stripped.Sequence")])
              test_ <-  lapply(params_report$comparisons, 
                    build_df_result,
                    data= res_de_usage$q_feat  ,
                    df_anno = df_ann  ,
                    mapping_df = fastaproc$result,
                    layer=  'precursors_lip_norm' ,
                    layer_ = 'precursors_lip_usage',
                  params = params_report)  
              assays_to_keep <- c("precursors_lip_norm", "precursors_tc_norm",'precursors_lip_usage' )
                # Subset to only these
              qf_final <- res_de_usage$q_feat[, , assays_to_keep]
        }else{
           # upaired 
              log_info('Unpaired branch ...')
                a <-  msqrob_model(pe = qf_final, params = params_report, layer = 'precursors_lip_norm' )
                if (a$status == 1) stop(a$error)
                b  <-  msqrob_model(pe = a$q_feat, params = params_report, layer = 'proteins_tc' )
                if (b$status == 1) stop(b$error)  
                qf_unpair <- calculate_lip_usage(b$q_feat, i_lip = "precursors_lip_norm", 
                                          i_tc = "proteins_tc",
                                           contrasts = colnames(b$contr_exp))
                if (qf_unpair$status == 1) stop(qf_unpair$error)  
                df_ann <-  as.data.frame(rowData(qf_unpair$qf[["precursors_lip_norm"]])[c("Precursor.Id", "Protein.Group", "Genes", "Proteotypic", "Stripped.Sequence")])
                test_ <-  lapply(params_report$comparisons, 
                    build_df_result,
                    data= qf_unpair$qf  ,
                    df_anno = df_ann  ,
                    mapping_df = fastaproc$result,
                    layer= 'precursors_lip_norm' ,
                    layer_ = NULL,
                  params = params_report)
                names(test_)<- params_report$comparison_label
                #qf_unpair$qf   

                assays_to_keep <- c("precursors_lip_norm", "precursors_tc_norm" )
                # Subset to only these
                qf_final <- qf_unpair$qf[, , assays_to_keep]
        }
        ##  qc -> data for QC
        ##  pe -> Qfeat subseted
       
        quarto_bag <- list( qc_data = qc_ann$result,
                            pe = qf_final,
                            res_DE= test_  )

        return( list(error= '', status= 0,result =quarto_bag ))
 
  
}

#' Build Result Data Frame for a Single Contrast
#'
#' Internal helper function to extract, format, and standardize differential expression 
#' results for a specific contrast. It handles both "paired" (usage in a separate assay) 
#' and "unpaired" (usage calculated within the same assay) workflows by unifying column names.
#'
#' @param label \code{character(1)}. The name of the contrast to extract (e.g., "TreatmentA - TreatmentB").
#' @param data A \code{QFeatures} object containing the statistical results in its \code{rowData}.
#' @param layer \code{character(1)}. The name of the main assay (usually the normalized LiP assay) 
#'   containing the protein-level or peptide-level fold changes.
#' @param layer_ \code{character(1)} or \code{NULL}. The name of the secondary assay (usually the 
#'   Usage assay). 
#'   \itemize{
#'     \item If provided (Paired analysis), results are fetched from this layer and joined.
#'     \item If \code{NULL} (Unpaired analysis), the function looks for usage statistics 
#'           (e.g., \code{pval_usage}) inside \code{layer} and renames them to match 
#'           the standard format (\code{usage_pval}).
#'   }
#' @param df_anno \code{data.frame}. A data frame containing annotation metadata. Must contain 
#'   a column named \code{"Precursor.Id"} for joining.
#'
#' @return A \code{list} containing a single element \code{df}, which is the combined 
#'   data frame of statistics and annotations for the requested contrast.
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr rename rename_with left_join mutate full_join relocate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove fixed
#' @keywords internal


build_df_result <- function (label, data , layer , layer_ = NULL , df_anno, mapping_df, params){
 # --- 1. Process Assay A (e.g., Lip normalized) ---
      res_layer <-  rowData(data[[layer]])[[label]]
      res_layer_df <- as.data.frame(res_layer, check.names = FALSE) %>% 
      rownames_to_column(var = "Precursor.Id") 
# --- 2. Process Assay B (Optional, e.g., Protein/Norm) ---
  if (!is.null(layer_)) {
     ## usage for paired 
      res_layer_ <- rowData(data[[layer_]])[[label]]
      res_layer__df <- as.data.frame(res_layer_, check.names = FALSE) %>% 
      rownames_to_column(var = "Precursor.Id") %>% 
      rename_with(~paste0("usage_", str_remove(.x, fixed(paste0(label, ".")))), .cols = -Precursor.Id)
      browser()
    # Merge A and B
    res_combined <- full_join(res_layer_df, res_layer__df, by = "Precursor.Id")  
    }else {
    # === UNPAIRED BRANCH ===
    # Usage is already in res_layer_df, but with suffixes (_usage).
    # We rename them to match the paired prefixes (usage_).
    res_combined <- res_layer_df %>%
      dplyr::rename(
        usage_logFC   = usage,            # 'usage' is the logFC
        usage_pval    = pval_usage,
        usage_adjPval = adjPval_usage,
        usage_se      = se_usage,
        usage_t       = t_usage,
        usage_df      = df_usage
      )
  }
  # --- 3. Final Join with Annotation ---
  final_df <- res_combined %>% 
    left_join(df_anno, by = "Precursor.Id") %>%
    mutate(contrast = label) # Good for downstream filtering
  
  # --4 join with sequence 
  final_df <- final_df %>% 
    left_join(mapping_df, by = join_by ( Protein.Group == Accession)) %>%
    pep_char()
   
  # top table corrected and not 
  full_toptable <- final_df %>% relocate("Precursor.Id", contains("usage"))
  
  # first filter the data              
  #vol_volcano_make(final_df)  
  # volcano = volcano_out,  POI_l = POI_Plot 

    res_usage  <- make_volcano_plot (
      df = full_toptable,
      params = params,
      title = paste0("Volcano ", 'TEST' ),
      to_save= FALSE
    )
  
    res_ncorrect   <- make_volcano_plot (
      df = full_toptable,
      params = params,
      title = paste0("Volcano ", 'TEST' ),
      usage = FALSE,
      to_save= FALSE
    )
  prot_seq_ <-  mapping_df %>% filter(Accession %in% res_usage$POI_l)
  
  se <- joinAssays(data,i=c("precursors_lip_norm","precursors_tc_norm"),  fcol = "Precursor.Id") %>%  
          getWithColData("joinedAssay")
  
  
   barplot_id <- plot_barcode (res_usage$POI_l, se , 
                  DE_result = full_toptable,   prot_seq = prot_seq_ ,  
                group_column = 'Treatment'  , params= params )
  
  return(  list( full_toptable =    full_toptable, 
                 plotvolcano_usage=  res_usage$volcano , 
                 plotvolcano_ncorr= res_ncorrect$volcano ,
                 barplot_ =  barplot_id$gg ))

}


#' Calculate LiP-MS Usage (Corrected logFC)
#'
#' This function corrects LiP peptide log-fold changes (logFC) by subtracting the 
#' total protein (TC) logFC. It performs error propagation for standard errors 
#' and re-calculates p-values for the "usage" (the peptide change relative to 
#' the protein change).
#'
#' @param qf A \code{QFeatures} object.
#' @param i_lip \code{character(1)}: The name of the assay containing LiP-level models.
#' @param i_tc \code{character(1)}: The name of the assay containing TC-level models.
#' @param contrasts \code{character()}: A vector of contrast names (e.g., "A - B") 
#'   that exist in the \code{rowData} of both assays.
#' @param fcol_lip \code{character(1)}: Column in \code{rowData(qf[[i_lip]])} 
#'   used for matching to proteins. Default is "Protein.Group".
#' @param fcol_tc \code{character(1)}: Column in \code{rowData(qf[[i_tc]])} 
#'   identifying proteins. Default is "Protein.Group".
#' @param satterthwaite \code{logical(1)}: Whether to use the Satterthwaite 
#'   approximation for degrees of freedom. Default is \code{FALSE}.
#'
#' @return A \code{QFeatures} object with updated \code{rowData} in the \code{i_lip} assay.
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom stats pt p.adjust

calculate_lip_usage <- function(qf, i_lip, i_tc, contrasts, 
                               fcol_lip = "Protein.Group", 
                               fcol_tc = "Protein.Group",  
                               satterthwaite = FALSE) {
  
   tryCatch({
    log_info('Unparied compute usage ...')
  # 1. Match precursors to proteins
  # This creates an index mapping each LiP row to its corresponding TC row
    match_idx <- match(SummarizedExperiment::rowData(qf[[i_lip]])[[fcol_lip]],
                     SummarizedExperiment::rowData(qf[[i_tc]])[[fcol_tc]])
  
  for (contrast in contrasts) {
    # 2.1. Extract TC results for this contrast
      log_info('Extract TC results  ...')
    res_tc <- SummarizedExperiment::rowData(qf[[i_tc]])[[contrast]]
    
    if (is.null(res_tc)) {
      warning(paste("Contrast", contrast, "not found in TC assay. Skipping."))
      next
    }
    
    # Rename TC columns to avoid collisions
    colnames(res_tc) <- paste0(colnames(res_tc), "_tc")
    
    # 2.2. Align TC results to the LiP precursors
    # We subset the TC results using the match index
    log_info('Align TC results to the LiP precursors  ...')
    aligned_tc <- res_tc[match_idx, , drop = FALSE]
    
    # 2.3. Combine and calculate usage
    # We convert to a standard data frame temporarily for easier calculation
    res_lip <- SummarizedExperiment::rowData(qf[[i_lip]])[[contrast]]
    combined <- cbind(res_lip, aligned_tc)
    
    # Calculation logic
     log_info('Compute usage  logFC, se, df ...')
    combined$usage <- combined$logFC - combined$logFC_tc
    combined$se_usage <- sqrt(combined$se^2 + combined$se_tc^2)
    combined$t_usage <- combined$usage / combined$se_usage
    
    # Degrees of freedom calculation
    if (satterthwaite) {
      combined$df_usage <- (combined$se^2 + combined$se_tc^2)^2 / 
        (combined$se^4/combined$df + combined$se_tc^4/combined$df_tc)
    } else {
      combined$df_usage <- combined$df
    }
    
    # P-value and FDR
    combined$pval_usage <- stats::pt(abs(combined$t_usage), 
                                     df = combined$df_usage, 
                                     lower.tail = FALSE) * 2
    combined$adjPval_usage <- stats::p.adjust(combined$pval_usage, method = "fdr")
    
    # Update the rowData of the original object
    log_info('Update rowData with usage corrected  ...')

    SummarizedExperiment::rowData(qf[[i_lip]])[[contrast]] <- combined
  }

  return( list(error= '', status= 0, qf =qf ))
  }, error = function(e) {
        print(paste("Unpaired calculate_lip_usage  :  ",err))
        return( list(error= err, status= 1, qf =NULL ))

  })
}

##-----------------

#' @author Andrea Argentini
#' @title check_design_requirement
#' @description This function performs the following checks on the design file:
#'  1) Check required columns exist
#'  2) Check Run column has no file extensions
#'  3) Check Pipeline column values
#' @param df Input data frame where features need to be checked
#' @param required_cols Character vector of required column names
#' @return A list with elements:
#'   \item{status}{integer; 0 if no error, 1 if an error was found}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#' @importFrom logger log_info


check_design_requirement <- function(df, required_cols) {
  # Default result = no error
  result <- list(error = "", status = 0)
  
  # 1. Check required columns exist
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    result$error <- paste("Missing required columns:", paste(missing_cols, collapse = ", "))
    result$status <- 1
    return(result)
  }
  
  # 2. Check Run column has no file extensions
  if ("Run" %in% names(df)) {
    bad_runs <- grep("\\.[A-Za-z0-9]+$", df$Run, value = TRUE)
    if (length(bad_runs) > 0) {
      result$error <- paste(
        "Column 'Run' contains values with file extensions:",
        paste(unique(bad_runs), collapse = ", ")
      )
      result$status <- 1
      return(result)
    }
  }
  
  # 3. Check Pipeline column values
  if ("Pipeline" %in% names(df)) {
    allowed <- c("LiP", "TC")
    bad_vals <- setdiff(unique(df$Pipeline), allowed)
    if (length(bad_vals) > 0) {
      result$error <- paste(
        "Column 'Pipeline' contains invalid values:",
        paste(bad_vals, collapse = ", "),
        "| Allowed values are:", paste(allowed, collapse = ", ")
      )
      result$status <- 1
      return(result)
    }
  }
  # no error detected
  return (result)
}


#' Create QFeatures Object from Proteomics Reports
#'
#' @title create_qfeat_
#' 
#' @description 
#' Converts raw data frames into a QFeatures object, calculates group-specific 
#' detection heuristics, and performs initial quality filtering.
#'
#' @param annotation_df A data frame containing sample metadata (colData).
#' @param report_file1 A data frame containing quantitative proteomics data.
#' @param report_file2 Optional; secondary report file (default is NULL).
#'
#' @return A list with status, error message, and the resulting QFeatures object.
#'
#' @import QFeatures
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom MultiAssayExperiment getWithColData
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats model.matrix
#' @importFrom matrixStats rowMins
#' @importFrom logger log_info
#'
create_qfeat_ <- function(annotation_df, report_file1, report_file2 = NULL) {

 calculatedSmallestGroupSize <- function(qf, i, facName){
    se <- MultiAssayExperiment::getWithColData(qf, i = i) #1.
    fac <- SummarizedExperiment::colData(se)[[facName]]
    SummarizedExperiment::rowData(qf[[i]])$smallestGroupSize <- 
      (se %>%
         SummarizedExperiment::assay() %>% #2
         Negate(is.na)()  #3
      ) %*% stats::model.matrix(~-1 + fac) %>% #4 
      matrixStats::rowMins() #5
    return(qf)
  }

  tryCatch(expr = {
    log_info('Importing data into QFeatures ...')

    # 1. Initial Import and subsetting
    # Filter using .data to avoid global variable warnings
    input_data <- report_file1 %>%
      dplyr::filter(.data$Precursor.Quantity > 4)

    qf <- QFeatures::readQFeatures(
      assayData = input_data,
      colData = annotation_df,
      quantCols = "Precursor.Quantity",
      runCol = "Run",
      fnames = "Precursor.Id"
    )

    # 2. Add Proteotypic info
    for (i in seq_along(qf)) {
      rd <- SummarizedExperiment::rowData(qf[[i]])
      rd$Proteotypic <- ifelse(grepl(";", rd$Protein.Group), 0, 1)
      SummarizedExperiment::rowData(qf[[i]]) <- rd
    }

    log_info('Filtering Qvalue/PG.Qvalue < 0.01 ...')

    # filterFeatures uses a formula interface; it handles its own variable scope
    qf <- QFeatures::filterFeatures(qf, ~ Q.Value <= 0.01 & 
                                      PG.Q.Value <= 0.01 & 
                                      Lib.Q.Value <= 0.01 & 
                                      Precursor.Id != "" & 
                                      Decoy == 0)

    # 3. Join assays based on Pipeline metadata
    lip_samples <- which(SummarizedExperiment::colData(qf)$Pipeline == "LiP")
    tc_samples <- which(SummarizedExperiment::colData(qf)$Pipeline == "TC")

    log_info('Joining Tc and Lip samples ...')
    qf <- QFeatures::joinAssays(x = qf, i = lip_samples, fcol = "Precursor.Id", name = "precursors_lip")
    qf <- QFeatures::joinAssays(x = qf, i = tc_samples,  fcol = "Precursor.Id", name = "precursors_tc")

    log_info('Filtering TC based on smallest group heuristic ...')
    qf <- calculatedSmallestGroupSize(qf, i = "precursors_lip", facName = "Treatment")
    qf <- calculatedSmallestGroupSize(qf, i = "precursors_tc", facName = "Treatment")

    # Filter for precursors present in at least 2 replicates per group
    qf <- QFeatures::filterFeatures(qf, ~ smallestGroupSize >= 2, i = "precursors_tc")

    log_info('Filtering Lip with precursor > 80% detection ...')
    nObs <- 2
    n <- ncol(qf[["precursors_lip"]])
    qf <- QFeatures::filterNA(qf, i = "precursors_lip", pNA = (n - nObs) / n)

    return(list(error = '', status = 0, result = qf))
    
  }, error = function(err) {
    msg <- conditionMessage(err)
    message(paste("Creating Qfeat Error: ", msg))
    return(list(error = msg, status = 1, result = NULL))
  })
}



#' Normalize and Aggregate Proteomics Data
#'
#' @title Log-Transformation, Median Normalization, and Protein Aggregation
#' 
#' @description 
#' Performs log2 transformation on LiP and TC assays, calculates sample-based 
#' normalization factors using common features, and aggregates TC precursors 
#' to the protein level.
#' 
#' The normalization scaling factor is calculated by:
#' \enumerate{
#'   \item Extracting the assay data.
#'   \item Removing features with missing values.
#'   \item Calculating column-wise medians to obtain log2 scale normalization factors.
#'   \item Zero-centering the normalization factors.
#'   \item Subtracting factors from intensities using \code{QFeatures::sweep}.
#' }
#'
#' @param q_feat A \code{QFeatures} object containing "precursors_lip" and "precursors_tc" assays.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{error}: Character string containing error messages, if any.
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{result}: The updated \code{QFeatures} object.
#' }
#'
#' @import QFeatures
#' @importFrom SummarizedExperiment assay
#' @importFrom matrixStats colMedians
#' @importFrom MsCoreUtils medianPolish
#' @importFrom logger log_info
#' @importFrom rlang .data
#' @importFrom stats median na.exclude
#'

normalization_scaling_factor <- function(q_feat ){

  # log_transform

  # 1. Extracts the assay data 
  # 2. Removes features with missing values
  # 3. Subsequently takes the column wise median to obtain the sample based normalisation factor on the log2 scale.
  # 4. Zero center the normalisation factors
  # 5. Subtract these log2-norm factors from the intensities of each corresponding column of the assay data and store the result in the new assay peptides_norm. (We adopt the sweep function to the peptides_log assay of the spikein QFeatures object with as statistic the log2 normfactor STATS=nf the default function FUN = "-", MARGIN = 2 to substract the column wise log2 norm factor from each entry of the corresponding assay data)

  # Internal helper to calculate median normalization factors
  medianNormCommonFeatures <- function(qf, i, name = "normAssay") {
    # 1. Extract assay, 2. Remove missing, 3. Calculate column medians
    m <- SummarizedExperiment::assay(qf[[i]])
    m_complete <- stats::na.exclude(m)
    
    norm_factors <- matrixStats::colMedians(m_complete)
    
    # 4. Zero-center the normalization factors
    norm_factors <- norm_factors - stats::median(norm_factors)
    
    # 5. Sweep out the factors
    qf <- QFeatures::sweep(
      qf, 
      MARGIN = 2, 
      STATS = norm_factors, 
      i = i, 
      name = name,
      FUN = "-"
    )
    return(qf)
  }

  tryCatch(
    expr = {
      log_info('Log2 transformation TC and Lip ...')
      q_feat <- QFeatures::logTransform(q_feat, i = "precursors_lip", name = "precursors_lip_log")
      q_feat <- QFeatures::logTransform(q_feat, i = "precursors_tc", name = "precursors_tc_log")

      log_info('Normalization features based on median common features ...')
      q_feat <- medianNormCommonFeatures(q_feat, i = "precursors_lip_log", name = "precursors_lip_norm")
      q_feat <- medianNormCommonFeatures(q_feat, i = "precursors_tc_log", name = "precursors_tc_norm")

      log_info('TC aggregation protein level ...')
      # Aggregates precursors to protein level using Median Polish
      q_feat <- QFeatures::aggregateFeatures(
        q_feat,
        i = "precursors_tc_norm",
        fcol = "Protein.Group",
        name = "proteins_tc",
        fun = MsCoreUtils::medianPolish,
        na.rm = TRUE
      )

      return(list(error = '', status = 0, result = q_feat))
    },
    error = function(err) {
      msg <- conditionMessage(err)
      message(paste("normalization_scaling_factor Error: ", msg))
      # Return the object in its current state even if error occurs
      return(list(error = msg, status = 1, result = q_feat))
    }
  )
}

####----

#' Compute LiP-MS Usage (Accessibility)
#'
#' @title compute_usage
#' 
#' @description 
#' This function calculates the relative peptide usage by subtracting the 
#' corresponding protein-level abundance (Total Control) from normalized 
#' peptide intensities. This step adjusts for changes in total protein expression 
#' to isolate changes in protein conformation or accessibility.
#'
#' @param q_feat A \code{QFeatures} object.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{error}: Character string containing error messages, if any.
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{result}: The updated \code{QFeatures} object containing the "precursors_lip_usage" assay.
#' }
#'
#' @import QFeatures
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom rlang .data
#' @importFrom logger log_info
#' @importFrom methods as
#'
compute_usage <- function(q_feat ){

calculate_usage_paired <- function(qf, i_lip = "precursors_lip_norm", 
                                    i_tc = "proteins_tc", 
                                    fcol_lip = "Protein.Group", 
                                    fcol_tc = "Protein.Group", 
                                    name = "precursors_lip_usage", 
                                    match_cols_lip_to_tc = colnames(qf[[i_lip]]))
{
    qf <-  QFeatures::addAssay(qf, qf[[i_lip]], name = name) #1.
    qf <- QFeatures::filterFeatures(qf, 
                         i = name,
                         filter = formula(paste0("~",fcol_tc, "%in% rowData(qf[[i_tc]])[[fcol_tc]]") ))

    match_precursor_to_protein_FC <- match(SummarizedExperiment::rowData(qf[[name]])[[fcol_lip]],
                       SummarizedExperiment::rowData(qf[[i_tc]])[[fcol_tc]]) 
    
    SummarizedExperiment::assay(qf[[name]]) <- 
      SummarizedExperiment::assay(qf[[name]]) - 
      SummarizedExperiment::assay(qf[[i_tc]])[match_precursor_to_protein_FC, match_cols_lip_to_tc]
    return(qf) 
} 

   tryCatch( expr = {
     log_info('Computing usage in paired design ...')
     q_feat <- calculate_usage_paired(q_feat, match_cols_lip_to_tc = 1:ncol(q_feat[["precursors_lip_norm"]]))
      return( list(error= '', status= 0,result =q_feat ))
   },error = function(err){
        print(paste("compute usage/paired design :  ",err))
        return( list(error= err, status= 1,result = NULL ))
  } )
 
}

#' Classify Peptide Trypticity
#' @title classify_trypticity
#' Categorizes peptides as Tryptic, Semi-Tryptic, or Non-Tryptic based on the 
#' presence of Lysine (K) or Arginine (R) at the cleavage sites, while 
#' accounting for protein termini and N-terminal Methionine excision.
#'
#' @param peptide A character vector of peptide sequences.
#' @param protein A character vector of the parent protein sequences.
#' @param start_pos A numeric vector indicating the starting position of the 
#' peptide within the protein (1-based indexing).
#'
#' @return A character vector of the same length as `peptide` containing 
#' "Tryptic", "Semi-Tryptic", "Non-Tryptic", or "Ambiguous".
#' 
#' @details 
#' The classification rules are:
#' \itemize{
#'   \item \strong{Tryptic}: Both ends follow tryptic rules (preceded by K/R or at N-term; ends in K/R or at C-term).
#'   \item \strong{Semi-Tryptic}: Only one end follows tryptic rules.
#'   \item \strong{Non-Tryptic}: Neither end follows tryptic rules.
#' }
#' Special case: If a peptide starts at position 2 and the first amino acid 
#' of the protein is Methionine (M), the N-terminus is considered tryptic 
#' due to common N-terminal Methionine excision.
#'
  
classify_trypticity <- function(peptide, protein, start_pos) {
 
  pep_len <- sapply(peptide, nchar)
  end_pos <- start_pos + pep_len - 1
 
    # Get surrounding amino acids
  AA_before <- ifelse(start_pos > 1, 
                    substr(protein, start_pos - 1, start_pos - 1), 
                    "")
  AA_after <- ifelse(end_pos < nchar(protein), 
                    substr(protein, end_pos + 1, end_pos + 1), 
                    "")
  AA_first <- substr(peptide, 1, 1)
  AA_last  <- substr(peptide, pep_len, pep_len)

 
  # Classification
  ambiguous <- is.na(protein) | is.na(start_pos)

  nterm_tryptic <- (AA_before %in% c("K", "R", "") | (AA_before == "M" & start_pos==2))
  cterm_tryptic <- (AA_last %in% c("K", "R") | AA_after=="")
  return(
    ifelse(ambiguous, 
         "Ambiguous",
         ifelse(nterm_tryptic & cterm_tryptic,
                "Tryptic",
                ifelse(nterm_tryptic | cterm_tryptic,
                       "Semi-Tryptic",
                       "Non-Tryptic")
                )
    )
  )
}
  
#' Calculate Protein Sequence Coverage
#' @title calculate_coverage
#' @author Andrea Argentini
#' This function calculates the fraction of a protein sequence covered by a set 
#' of peptides or fragments. It accounts for overlapping regions by treating 
#' the segments as genomic-style ranges.
#'
#' @param start A numeric vector of start positions for the peptides.
#' @param end A numeric vector of end positions for the peptides.
#' @param protein_length A single numeric value representing the total length 
#' of the protein sequence.
#'
#' @return A numeric value (0 to 1) representing the percentage of the protein 
#' sequence covered. Returns `NA` if `protein_length` is `NA` or if no valid 
#' start/end pairs are provided.
#'
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors coverage
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' 

 calculate_coverage <- function(start, end, protein_length) {
  #Checks 
  if (is.na(protein_length)) return(NA)
  startEnd <- cbind(start,end) %>% na.omit()
  if (nrow(startEnd) <1) return(NA)
  covered <- IRanges(startEnd[,1], startEnd[,2]) %>%
    coverage() %>%
    as.logical() %>%
    sum()
  return(covered/protein_length)
} 

#' Characterize Peptide Properties and Mapping
#'
#' @title pep_char
#' 
#' @description 
#' This function annotates a data frame of peptides with biochemical and 
#' positional properties. Finally, we define a function for adding and 
#' adjusting data for the features to:
#' \enumerate{
#'   \item Define missed_cleavages
#'   \item Calculate how many times precursor is repeated in protein
#'   \item Calculate start and end position in a protein
#'   \item Define peptide type
#'   \item Extract last amino acid
#' }
#'
#' @param table A data frame or tibble containing peptide and protein sequences.
#' @param prot_seq A character string specifying the column name for the 
#' full protein sequence. Default is `"Protein.Sequence"`.
#' @param pep_seq A character string specifying the column name for the 
#' stripped peptide sequence. Default is `"Stripped.Sequence"`.
#'
#' @return A data frame with the following additional columns:
#' \itemize{
#'   \item \code{missed_cleavages}: Count of internal [RK] not followed by P.
#'   \item \code{total_repeats}: Number of times the peptide occurs in the protein.
#'   \item \code{start}: Starting position of the first occurrence.
#'   \item \code{end}: Ending position of the first occurrence.
#'   \item \code{pep_type}: Trypticity classification (Tryptic, Semi, Non, or Ambiguous).
#'   \item \code{AA_last}: The C-terminal amino acid of the peptide.
#' }
#'
#' @importFrom dplyr mutate select
#' @importFrom stringr str_count str_locate str_sub
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#'
pep_char <- function(table, prot_seq="Protein.Sequence", pep_seq="Stripped.Sequence")
{
  table <- table |>
    mutate(missed_cleavages = get(pep_seq) %>%
             str_count("[RK](?!(P|$))"), #1.
         total_repeats = str_count(get(prot_seq), get(pep_seq)), #2.
         tmp = str_locate(Protein.Sequence, Stripped.Sequence),
         start = tmp[,1],
         end = tmp[,2], #3.
         pep_type = ifelse(
           (total_repeats > 1) | is.na(total_repeats),
           "Ambiguous",
           classify_trypticity(
             peptide = Stripped.Sequence, 
             protein = Protein.Sequence, 
             start_pos = start)), #4.
         AA_last = str_sub(Stripped.Sequence, -1, -1) #5.
           ) %>%
    select(-tmp)
} 


#' Annotate Precursors for QC Visualization
#'
#' @title Precursor Annotation for Quality Control Plots
#' 
#' @description 
#' Converts QFeatures assays into a long-format data frame and performs 
#' comprehensive annotation including protein mapping, peptide characterization, 
#' and condition formatting. This prepared data is typically used for QC 
#' plots like MDS or intensity distributions.
#'
#' @param q_feat A \code{QFeatures} object containing "precursors_lip_norm", 
#' "precursors_tc_norm", and "precursors_lip_usage" assays.
#' @param mapping A data frame used for joining protein accessions to metadata. 
#' Must contain an \code{Accession} column.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{error}: Character string containing error messages, if any.
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{result}: An annotated data frame in long format.
#' }
#'
#' @importFrom QFeatures longForm
#' @importFrom dplyr left_join mutate case_match join_by
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom logger log_info
#'


qc_precursor_annotation <- function(q_feat, mapping, type) {
  tryCatch(expr = { 
    log_info('Annotate precursor for QC plot ...')
    
    # 1. Define layers
    layer <- if (type == 'paired') {
      c("precursors_lip_norm", "precursors_tc_norm", "precursors_lip_usage")
    } else {
      c("precursors_lip_norm", "precursors_tc_norm")
    }
    
    # 2. Extract and Join
    qcObj <- q_feat[,,layer] %>%
      QFeatures::longForm(
        colvars = c("Condition", "CondRep", "Treatment", "Replicate", "Pipeline"), 
        rowvars = c("Precursor.Id", "Protein.Group", "Stripped.Sequence", 
                    "Precursor.Charge", "Genes")
      ) %>%
      as.data.frame() %>%
      # Fix: Removed x$ and y$ references for join_by compatibility
      dplyr::left_join(mapping, by = dplyr::join_by(Protein.Group == Accession)) %>%
      pep_char() %>%
      # Fix: Moved this OUTSIDE the 'if' so all types get clean labels
      dplyr::mutate(
        assay = dplyr::case_match(
          .data$assay,
          "precursors_lip_norm" ~ "LiP", 
          "precursors_lip_usage" ~ "usage", 
          "precursors_tc_norm" ~ "TC",
          .default = .data$assay
        ),
        CondRep = paste(.data$Condition, .data$assay, .data$Replicate, sep = "_"),
        Condition = paste(.data$Condition, .data$assay, sep = "_")
      )

    return(list(error = '', status = 0, result = qcObj))
    
  }, error = function(err) {
    msg <- conditionMessage(err)
    message(paste("Annotation precursor error: ", msg))
    return(list(error = msg, status = 1, result = NULL))
  })
}




#' @author Andrea Argentini
#' @title annotate_diann_OLD
#' @description Annotate and standardize DIA-NN reports: bind reports, join sample annotations,
#' compute peptide-level metrics (missed cleavages, proteotypic), expand protein accessions,
#' join FASTA-derived sequence information, compute peptide start/end positions and flanking AAs,
#' and classify peptides as Tryptic / SemiTryptic / NonTryptic. Returns an annotated data.frame inside a result list.
#' @param annotation_df Data frame with sample annotation (must contain "Run" for joining)
#' @param report_file1 Data frame or tibble with DIA-NN report data
#' @param report_file2 Optional second DIA-NN report data frame or tibble (default NULL)
#' @param fasta_ann Data frame with FASTA-derived annotations (must contain "Accession" and "Protein.Sequence")
#' @return A list with elements:
#'   \item{status}{integer; 0 if no error, 1 if an error was found}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{result}{data.frame; annotated DIA-NN data (NULL if an error occurred)}
#' @importFrom arrow read_parquet read_tsv_arrow
#' @importFrom tidyr separate_rows
#' @importFrom dplyr rename select bind_rows left_join group_by mutate ungroup filter slice n
#' @importFrom stringr str_split str_count str_locate_all str_sub str_trim str_detect str_remove_all
#' @importFrom magrittr %>%
#' @importFrom logger log_info

 annotate_diann_OLD <- function(annotation_df, report_file1, report_file2 = NULL,fasta_ann){

     tryCatch( expr = {
    log_info('Annotation DIA-NN standardize column name ...')

    app <- bind_rows(report_file1, report_file2) %>%
            select(Run,
              Precursor.Id,
              Modified.Sequence,
              Stripped.Sequence,
              Precursor.Charge,
              Protein.Group,
              Protein.Names,
              Genes,
              Precursor.Quantity,
              Ms1.Area,
              # Q.Value, #If no MBR is done, all precursors pass the threshold in set in DIA-NN, there's no need to filter further
              Lib.Q.Value) %>%
  left_join(annotation_df, by = "Run") %>% filter(Precursor.Quantity>4,
         Lib.Q.Value <= 0.01)
       
  app_a <- app %>% 
    group_by(Precursor.Id, Condition) %>%
    filter(n() >= 2)   #Keep only precursors which were found in at least 2 samples per condition


app_b <- app_a %>%
  mutate(
    Genes = str_remove_all(Genes, "_.*?(?=;|$)"),
    missed_cleavages = str_count(Stripped.Sequence, "[RK](?!(P|$))"),
    Proteotypic     = ifelse(str_detect(Protein.Group, ";"), 0, 1),
    Accession       = Protein.Group         # create column to split so original stays
  ) %>%
  separate_rows(Accession, sep = ";") %>%
  mutate(Accession = str_trim(Accession)) %>%  # trim spaces if any
  ungroup()
  ## splited in two parts
  log_info('Annotation spectronaut adding fasta info ...')
  
  app_c <- app_b %>%  left_join(fasta_ann, by = "Accession") %>%   # option if fasta file is used
  mutate(total_repeats = str_count(Protein.Sequence, Stripped.Sequence)) %>%
  group_by(Run, Stripped.Sequence, Accession) %>%
  dplyr::slice(rep(1, total_repeats[1])) %>%
  mutate(repeat_nr = 1:max(total_repeats)) %>%
  group_by(Stripped.Sequence, Accession) %>%
  mutate(start = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 1],
         end = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 2],
         AA_before = str_sub(Protein.Sequence[1], start-1, start-1),
         AA_last = str_sub(Stripped.Sequence[1], -1, -1),
         AA_after = str_sub(Protein.Sequence[1], end+1, end+1),
         pep_type = ifelse(AA_before %in% c("K", "R") & (AA_last %in% c("K", "R")) | #internal tryptic peptides
                          (AA_before %in% c("K", "R") & AA_after == "") | #C-terminal tryptic peptides
                          ((AA_before == "" | (AA_before == "M" & start==2)) & AA_last %in% c("K", "R")), #N-terminal tryptic peptides
                           "Tryptic",
                          ifelse(AA_before %in% c("K", "R") | (AA_last %in% c("K", "R")),
                                 "SemiTryptic",
                                 "NonTryptic"))) %>% ungroup()

  return( list(error= '', status= 0, result =app_c ))
     },error = function(err){
    print(paste("Annotation on Spectrounaut :  ",err))
    return( list(error= err, status= 1,result =NULL ))
  } )

}

#' @author Andrea Argentini
#' @title parse_input
#' @description Read and parse input parquet file(s) and the experiment design file. Detects whether reports are in DIA‑NN format, reads TC and/or LiP reports (from one or two parquet files depending on 'dual'), validates the design using check_design_requirement, and returns parsed reports along with the design and a diann_flag.
#' @param input_parquet_tc Path to the TC parquet file (used when dual = TRUE)
#' @param input_parquet_lip Path to the LiP parquet file (or the single combined parquet when dual = FALSE)
#' @param dual Logical; TRUE if TC and LiP are in separate parquet files, FALSE if both are in one file
#' @param input_design Path to the experiment design file (tsv/arrow readable)
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{lip}{data.frame or tibble; LiP report data (NULL on error)}
#'   \item{tc}{data.frame or tibble; TC report data or NULL if not provided}
#'   \item{design}{data.frame or tibble; parsed design file}
#'   \item{diann_flag}{logical; TRUE if the report looks like DIA‑NN (contains a "Run" column)}
#' @importFrom arrow read_parquet read_tsv_arrow
#' @importFrom dplyr rename
#' @importFrom rlang .data
#' @importFrom logger log_info
parse_input <- function(input_parquet_tc, input_parquet_lip, dual, input_design) {
  
  # Internal helper to check for DIA-NN format
  is_diann <- function(df) {
    "Run" %in% colnames(df)
  }
  
  tryCatch(
    expr = {
      # 1. Loading Reports
      if (isTRUE(dual)) {
        log_info('Reading Tc and Lip from SEPARATE parquet files ...')
        TC_report <- arrow::read_parquet(input_parquet_tc)
        LiP_report <- arrow::read_parquet(input_parquet_lip)
        # Note: 'stop' will be caught by the error block below
        stop('Dual mode logic is currently being fixed.') 
      } else {
        log_info('Reading both LiP and TC from ONE parquet file ...')
        LiP_report <- arrow::read_parquet(input_parquet_lip)
        TC_report <- NULL
      }
      
      diann_flag <- is_diann(LiP_report)
      
      # 2. Loading and Validating Design
      log_info('Reading experiment Design file ...')
      design <- arrow::read_tsv_arrow(input_design)
      
      col_design_required <- c('Run', 'Pipeline', 'Treatment', 'Condition', 'Replicate', 'CondRep')
      
      # Assuming check_design_requirement is an internal package function
      checkdesign <- check_design_requirement(design, col_design_required)
      
      # Use .data$Run to avoid "no visible binding for global variable" warning
      design <- design %>% 
        dplyr::rename(runCol = .data$Run)
      
      if (checkdesign$status == 1) {
        return(list(error = checkdesign$error, status = 1, lip = NULL))
      } else {
        return(list(
          error      = '', 
          status     = 0,
          lip        = LiP_report,
          tc         = TC_report,
          design     = design,
          diann_flag = diann_flag
        ))
      }
    },
    error = function(err) {
      # Extracting the error message specifically
      msg <- conditionMessage(err)
      message(paste("Input Parquet Error: ", msg))
      return(list(error = msg, status = 1, lip = NULL))
    }
  )
}


#' @author Andrea Argentini
#' @title annotate_spectronaut
#' @description Annotate and standardize Spectronaut reports: bind reports, rename/standardize columns,
#' join sample annotations, filter and keep precursors found in at least two samples per condition,
#' compute peptide-level metrics (missed cleavages, proteotypic), expand protein accessions,
#' join FASTA-derived sequence information, compute peptide start/end positions and flanking amino acids,
#' and classify peptides as Tryptic / SemiTryptic / NonTryptic. Returns an annotated data.frame inside a result list.
#' @param annotation_df Data frame with sample annotation (must contain "Run" for joining)
#' @param report_file1 Data frame or tibble with Spectronaut report data
#' @param report_file2 Optional second Spectronaut report data frame or tibble (default NULL)
#' @param fasta_ann Data frame with FASTA-derived annotations (must contain "Accession" and "Protein.Sequence")
#' @return A list with elements:
#'   \item{status}{integer; 0 if no error, 1 if an error was found}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{result}{data.frame; annotated Spectronaut data (NULL if an error occurred)}
#' @importFrom arrow read_parquet read_tsv_arrow
#' @importFrom tidyr separate_rows
#' @importFrom dplyr rename select bind_rows left_join group_by mutate ungroup filter slice n
#' @importFrom stringr str_split str_count str_locate_all str_sub str_trim str_detect
#' @importFrom magrittr %>%

#' @importFrom logger log_info
annotate_spectronaut <- function (annotation_df, report_file1, report_file2 = NULL,fasta_ann){

  tryCatch( expr = {
    log_info('Annotation spectronaut standardize column name ...')

    app <- bind_rows(report_file1, report_file2) %>%
          dplyr::rename(  "Run" = "R_FileName" ,
                 "Precursor.Id" = "PEP_GroupingKey" ,
                  "Modified.Sequence" = "FG_LabeledSequence" ,
                  "Stripped.Sequence" = "PEP_StrippedSequence" ,
                 "Precursor.Charge" = "FG_Charge" ,
                 "Protein.Group" = "PG_ProteinGroups",
                 "Genes" =  "PG_Genes",
                 "Protein.Names" = "PG_ProteinNames",
                 "Precursor.Quantity"= "FG_MS2RawQuantity"  ,
                  "Ms1.Area" = "FG_MS1RawQuantity",
                  "Lib.Q.Value" = "FG_Qvalue" ) %>%
        select(Run,
              Precursor.Id,
              Modified.Sequence,
              Stripped.Sequence,
              Precursor.Charge,
              Protein.Group,
              Protein.Names,
              Genes,
              Precursor.Quantity,
              Ms1.Area,
              # Q.Value, #If no MBR is done, all precursors pass the threshold in set in DIA-NN, there's no need to filter further
              Lib.Q.Value) %>%
  left_join(annotation_df, by = "Run")
app_a <- app %>%  filter(Precursor.Quantity> 4,
         Lib.Q.Value <= 0.01) %>% #MBR was used, so filter on Lib.Q.Value
  group_by(Precursor.Id, Condition) %>%
  filter(n() >= 2)   #Keep only precursors which were found in at least 2 samples per condition

# emin original code
#app_b <- app_a %>% mutate(missed_cleavages = (Stripped.Sequence %>% str_count("[RK](?!(P|$))"))) %>%
#          group_by(Run, Precursor.Id) %>%
#          mutate(Proteotypic = ifelse(str_detect(Protein.Group, ";"), 0, 1)) %>%
#          slice(rep(1, str_count(Protein.Group, ";")+1)) %>%
#          mutate(Accession = unlist(str_split(Protein.Group[1], ";"))) %>%
#          ungroup()

app_b <- app_a %>%
  mutate(
    missed_cleavages = str_count(Stripped.Sequence, "[RK](?!(P|$))"),
    Proteotypic     = ifelse(str_detect(Protein.Group, ";"), 0, 1),
    Accession       = Protein.Group         # create column to split so original stays
  ) %>%
  separate_rows(Accession, sep = ";") %>%
  mutate(Accession = str_trim(Accession)) %>%  # trim spaces if any
  ungroup()
  ## splited in two parts
  log_info('Annotation spectronaut adding fasta info ...')
  
  app_c <- app_b %>%  left_join(fasta_ann, by = "Accession") %>%   # option if fasta file is used
  mutate(total_repeats = str_count(Protein.Sequence, Stripped.Sequence)) %>%
  group_by(Run, Stripped.Sequence, Accession) %>%
  dplyr::slice(rep(1, total_repeats[1])) %>%
  mutate(repeat_nr = 1:max(total_repeats)) %>%
  group_by(Stripped.Sequence, Accession) %>%
  mutate(start = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 1],
         end = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 2],
         AA_before = str_sub(Protein.Sequence[1], start-1, start-1),
         AA_last = str_sub(Stripped.Sequence[1], -1, -1),
         AA_after = str_sub(Protein.Sequence[1], end+1, end+1),
         pep_type = ifelse(AA_before %in% c("K", "R") & (AA_last %in% c("K", "R")) | #internal tryptic peptides
                          (AA_before %in% c("K", "R") & AA_after == "") | #C-terminal tryptic peptides
                          ((AA_before == "" | (AA_before == "M" & start==2)) & AA_last %in% c("K", "R")), #N-terminal tryptic peptides
                           "Tryptic",
                          ifelse(AA_before %in% c("K", "R") | (AA_last %in% c("K", "R")),
                                 "SemiTryptic",
                                 "NonTryptic"))) %>% ungroup()

  return( list(error= '', status= 0, result =app_c ))
  },error = function(err){
    print(paste("Annotation on Spectrounaut :  ",err))
    return( list(error= err, status= 1,result =NULL ))
  } )

}


#' @author Andrea Argentini
#' @title consensus_normalisation_OLD
#' @description Function to calculate scaling factors for median normalisation based on precursors that are identified in every sample. Computes per-sample medians on shared precursors, derives scaling factors, and applies them to produce normalized precursor intensities.
#' @param report Data frame or tibble containing at least the columns "CondRep", "Precursor.Id", and "Precursor.Quantity"
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{normalized}{data.frame; original report augmented with column normPQ (normalized log2 precursor quantity)}
#'   \item{scale_factor}{data.frame; per-sample medians and computed scaling factors}
#' @importFrom arrow read_parquet read_tsv_arrow
#' @importFrom dplyr distinct select group_by mutate ungroup n_distinct summarise
#' @importFrom logger log_info



consensus_normalisation_OLD <- function(report){

  tryCatch( expr = {

    log_info('Consensus normalization ...')

   report_long <- report %>%
      distinct(CondRep, Precursor.Id, .keep_all = TRUE) %>%
      select(CondRep, Pipeline, Precursor.Id, Precursor.Quantity) %>%
      group_by(Pipeline) %>% 
         mutate(samples = n_distinct(CondRep)) %>% 
      group_by(Precursor.Id) %>% 
      filter(n_distinct(CondRep)==samples) %>% 
      group_by(CondRep) %>% 
      mutate(log_Precursor.Quantity = log2(Precursor.Quantity))

  
    # Step 2: compute sample medians safely
    log_info('Computing scaling factor ...')
    scaling_factors <- report_long %>%
      
      summarise(sample_shared_median = median(log_Precursor.Quantity), .groups = "drop") %>%
      ungroup() %>% 
      mutate(sample_scaling_factor = median(sample_shared_median) - sample_shared_median)

    # Step 3: apply scaling factors to original report
    log_info('Correcting original intensity ...')

    report_norm <- report %>%
      mutate(normPQ = log2(Precursor.Quantity) +
               scaling_factors$sample_scaling_factor[match(CondRep, scaling_factors$CondRep)])

    return(list(error = '', status = 0, normalized = report_norm, scale_factor = scaling_factors))

  },error = function(err){
      print(paste(" Consensus normalization  :  ",err))
      return( list(error= err, status= 1,normalized =NULL ))
  } )
}



#' @author Andrea Argentini
#' @title calculate_coverages
#' @description Calculate sequence coverage per protein accession across the dataset.
#' For each accession the function reduces peptide ranges (start/end) and computes the
#' fraction of covered amino-acid positions relative to the reported protein length, then
#' joins the coverage value back to the original report.
#' @param report Data frame or tibble containing at minimum the columns:
#'   "Accession", "length" (protein length), "start" and "end" (peptide coordinates)
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{result}{data.frame; original report augmented with a numeric column "coverage" (per-accession coverage)}
#' @importFrom dplyr distinct left_join join_by group_by summarise
#' @importFrom IRanges IRanges reduce width
#' @importFrom logger log_info



calculate_coverages <- function(report) {

  tryCatch( expr = {
        log_info('Computing Sequence Coverage ...')

        coverages <-
          report %>%
          distinct(Accession, length, start, end) %>%
          group_by(Accession) %>%
          summarise(
            coverage = {
              ranges <- IRanges(start = start, end = end)
              covered_positions <- sum(width(reduce(ranges)))
              covered_positions / length[1]},
            .groups = "drop")
      report_ <- report %>%   left_join(coverages, join_by(Accession))
     return( list(error= '', status= 1,result =report_ ))
  },error = function(err){
    print(paste(" Coverage computation :  ",err))
    return( list(error= err, status= 1,result =NULL ))
  } )

}


#' @author Andrea Argentini
#' @title read_fasta_ann
#' @description Read a FASTA file and return a mapping of accession to protein sequence and sequence length. Extracts the accession (second field when splitting the FASTA header by "|") and computes protein length.
#' @param input_fasta Path to the FASTA file to read (character)
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{result}{data.frame; mapping with columns Accession, Protein.Sequence, length (NULL on error)}
#' @importFrom seqinr read.fasta
#' @importFrom stringr str_split_i
#' @importFrom logger log_info



read_fasta_ann <- function(input_fasta ){

tryCatch( expr = {
## code here
     log_info('Reading Fasta file ...')
    fasta <- read.fasta(input_fasta , seqtype = "AA", as.string = T)
  # read in human fasta file for protein sequence & length -> to calculate the coverage
    mapping <- data.frame(Accession = fasta %>% names() %>% str_split_i("\\|", 2),
                        Protein.Sequence = unlist(fasta))
    mapping$length <- nchar(mapping$Protein.Sequence)  #add length of the protein sequence
      return( list(error= '', status= 0,result =mapping ))
  ## good exit
},error = function(err){
    print(paste(" Reading Fasta  :  ",err))
    return( list(error= err, status= 1,result =NULL ))
  } )

}





#' @author Andrea Argentini
#' @title input_qf
#' @description Create a QFeatures object from a precursor-level input data.frame. The function pivots quantification to wide format, selects the appropriate quantification column depending on flag_tc (normPQ for TC, adjPQ for LiP), filters the design for the chosen pipeline, and returns the QFeatures precursor assay and the filtered design.
#' @param df_input Data frame or tibble containing precursor-level data (must include Run, Precursor.Id, Protein.Group, repeat_nr and the quantification columns like normPQ or adjPQ)
#' @param design Data frame with experiment design (must include Run, Pipeline, CondRep, etc.)
#' @param columns_not_wide Character vector of identifier/feature columns to retain prior to pivoting to wide format
#' @param flag_tc Logical; TRUE to process TC (uses normPQ), FALSE to process LiP (uses adjPQ). Default TRUE.
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{qf_pe}{QFeatures object; precursor-level QFeatures (NULL on error)}
#'   \item{design_filt}{data.frame; filtered design used as colData for the QFeatures object (NULL on error)}
#' @importFrom dplyr distinct
#' @importFrom tidyr drop_na
#' @importFrom stringr str_detect
#' @importFrom QFeatures readQFeatures
#' @importFrom logger log_info


input_qf <- function(df_input  , design , columns_not_wide, flag_tc = TRUE   ){

tryCatch( expr = {
  if (flag_tc  == TRUE){
      quantCol = "normPQ"
      filt_pipeline = 'TC'
      }else{
        quantCol = "adjPQ"
        filt_pipeline = 'LiP'
      }

  log_info(paste0('Creating Qfeature obj for ', filt_pipeline )) 

 input_wide <- df_input %>%
    filter( repeat_nr==1) %>%
    distinct(Run, Precursor.Id, Protein.Group, .keep_all = T) %>%
    dfToWideMsqrob(precursorquan = quantCol,
      wide_colums = columns_not_wide )  

  if (flag_tc  == TRUE){
    # for TC  no NAN allowed
     input_wide <- input_wide %>% drop_na()
 
  }

  
  if (flag_tc  == TRUE){
     design_filt <- design %>% filter(Pipeline  == filt_pipeline ) %>% rename(quantCols = Run)
      }else{
       ## lip 
      design_filt <- design %>% filter(Pipeline  == filt_pipeline ) %>% 
        rename(quantCols = Run) %>%
        rename(SampleName = CondRep  ) %>%  mutate(Condition = as.factor(Condition))
    }

  ## create qfeat obj
  qf_ <- readQFeatures( input_wide ,
            fnames = "Precursor.Id",
            quantCols =  str_detect(names(input_wide), paste( columns_not_wide , collapse = "|"), negate=TRUE) ,
            colData = design_filt,
            name = "precursor")
  
  

  return ( list( error= '', status= 0, qf_pe = qf_,  design_filt = design_filt ))
},error = function(err){
    print(paste(" QF Creation :  ",err))
    return( list(error= err, status= 1, qf_pe = NULL,  design_filt = NULL ))
  } )

}






#' @author Andrea Argentini
#' @title processing_tc_qfeat
#' @description Aggregate precursor-level QFeatures to proteins, filter proteins with more than one peptide,
#' compute per-protein median deviations and drug-specific abundance adjustment factors to support normalization.
#' The function returns a data.frame of adjustment factors that can be applied downstream.
#' @param input_pe QFeatures object (precursor-level) to aggregate (expects a "precursor" assay)
#' @param design Data frame with experiment design (must include a "Run" column and "Drug" column for grouping)
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{adj_scaling_df}{data.frame; per-protein and per-drug abundance adjustment factors (NULL on error)}
#' @importFrom dplyr distinct
#' @importFrom stringr str_detect
#' @importFrom QFeatures readQFeatures aggregateFeatures filterFeatures longFormat
#' @importFrom logger log_info
#' @importFrom stats median


processing_tc_qfeat <- function(input_pe  , design    ){

  tryCatch( expr = {

          log_info('Aggregating to protein level ...')
          input_pe_ <- aggregateFeatures(input_pe,
                                i = "precursor",
                                fcol = "Protein.Group",
                                name = "proteins",
                                fun = MsCoreUtils::robustSummary,
                                na.rm = TRUE)
          ## work aroung for the filtering
         

          prot <- input_pe_[["proteins"]]
          prot <- prot[rowData(prot)$.n > 1, ]
          input_pe_[["proteins"]] <- prot    
    
          log_info('Filtering more than one petides  ...')
          #input_pe_ <- filterFeatures(input_pe_ ,
          #                      ~.n >1 ,
          #                      "proteins",
          #                      keep = F)
          log_info('Computing scaling factor  ...')
           ## TODO  hard coded stuff 
          abundance_adjustment_factors <-
            longFormat(input_pe_[["proteins"]]) %>%
            rename( "Run" = "colname",
                     "Protein.Group"= "rowname")%>%
            left_join(design, join_by(Run)) %>%
            group_by(Protein.Group) %>%
            mutate(deviation = value - median(value, na.rm = TRUE)) %>%
            group_by(Protein.Group,  Condition, Treatment) %>%
            summarise(abundance_adjustment = median(deviation, na.rm = TRUE)) %>%
            mutate(ID = paste(Protein.Group, Treatment)) %>%
            ungroup()
    
      return ( list( error= '', status= 0 ,  adj_scaling_df = abundance_adjustment_factors  ))
  },error = function(err){
    print(paste(" TC computing scaling  :  ",err))
    return( list(error= err, status= 1,adj_scaling_df =NULL ))
  } )


  
}


#' @author Andrea Argentini
#' @title msqrob_model
#' @description Fit msqrob2 models and perform hypothesis testing for differential expression.
#' Given a QFeatures object, a parameter list (including a model formula and comparisons),
#' and the assay layer name, this function fits msqrob models, constructs contrasts from
#' the provided comparisons, runs hypothesis tests, and returns the QFeatures object with results.
#' @param pe QFeatures object containing the assay to model
#' @param params List of parameters; must include at least:
#'   \describe{
#'     \item{formula}{character or formula string used for msqrob modeling}
#'     \item{comparisons}{character vector of coefficient names to test (e.g. c("condB-condA"))}
#'   }
#' @param layer Character; name/index of the assay layer in the QFeatures object to analyze
#' @return A list with elements:
#'   \item{status}{integer; 0 if successful, 1 if an error occurred}
#'   \item{error}{character; error message when status is 1, otherwise an empty string}
#'   \item{q_feat}{QFeatures object; the input QFeatures augmented with fitted models and hypothesis test results (NULL on error)}
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom msqrob2 msqrob getCoef makeContrast hypothesisTest
#' @importFrom dplyr left_join select group_by summarise distinct n
#' @importFrom stats as.formula
#' @importFrom logger log_info


msqrob_model <- function(pe, params, layer  ){

  tryCatch( expr = {

    
    log_info('Msqrob model ...')
    
    pe <- msqrob(object = pe, i = layer,
          formula = as.formula(params$formula),

          ridge = FALSE,
          overwrite = T)

    contrast_list <- paste0(params$comparisons, "=0")

    coef <- rowData(pe[[layer]])$msqrobModels[[1]] %>% getCoef %>% names
    model_type <- sapply(rowData(pe[[layer]])$msqrobModels, function(x) x@type)
    success_model <- which(model_type != "fitError")
    if (length(success_model) == 0)
      return(list(
          error = 'Model is not able to converge', status = 1,
          q_feat = NULL, de_comp = NULL
      ))

    getCoef(rowData(pe[[layer]])$msqrobModels[[1]])

    log_info('Making contrast & testing ...')
    L <- makeContrast(contrast_list, parameterNames = coef)
    pe <- hypothesisTest(object = pe, i = layer, contrast = L , overwrite=TRUE)

    return (list(error= '', status= 0,q_feat = pe , contr_exp = L   ))


  },error = function(err){
    print(paste("Msqrob modeling :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

}


#' @author Andrea Argentini
#' @title make_volcano_plot
#' @description Generates a volcano plot from differential expression results using ggplot2.
#' Supports optional tooltip composition from annotation fields for interactive displays,
#' draws fold-change and significance guide lines, colours points by interest, and can
#' either return an interactive-style plot (with text tooltips prepared) or a static plot
#' with labeled significant proteins.
#' @param df Data.frame containing differential expression results (must include columns like logFC, pval, adjPval, interest, Protein.Group, significance)
#' @param params List of parameters; expected entries include numeric FC_thr (fold-change threshold) and adjpval_thr (adjusted p-value threshold)
#' @param title Character; title for the plot
#' @param annotation_fields Character vector; column names to use in tooltip annotations (e.g. c("Protein.Names", "Genes"))
#' @param poi_vis Character or numeric vector; points-of-interest used to determine colour palette length (e.g. vector of proteins/genes of interest)
#' @param to_save Logical; if FALSE prepare interactive-style tooltip text (for saving as interactive plot), if TRUE produce a static plot with labeled significant points
#' @return A ggplot object representing the volcano plot (interactive tooltip text is placed in a column named `tooltip_text` when to_save is FALSE)
#' @importFrom ggplot2 ggplot  scale_alpha_manual scale_shape_manual aes theme_bw geom_point geom_vline geom_hline scale_colour_manual labs
#' @importFrom ggsci pal_npg
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate select all_of filter distinct
#' @importFrom plotly plot_ly add_lines layout


make_volcano_plot <- function(df, params, title, annotation_fields = NULL,   to_save = FALSE, usage= TRUE  ) {
  # Compose tooltip string based on requested annotation fields
  if (to_save== F ){
    df <- df %>%
      mutate(
        tooltip_text = apply(
          select(., all_of(annotation_fields)),
          1,
          function(row) {
            paste(paste0(annotation_fields, ": ", row), collapse = "<br>")
          })
         )
    }
    if (usage){ 
      df <- df %>% select(.data$usage_adjPval,.data$usage_pval, .data$usage_logFC, .data$Protein.Group, .data$Genes ,.data$pep_type) %>%  
      mutate(adjPval = usage_adjPval,
           pval = usage_pval,
           logFC = usage_logFC) 
    }else{
           df <- df %>% select(.data$adjPval,.data$pval, .data$logFC, .data$Protein.Group, .data$Genes ,.data$pep_type) 
    }
    # plotly native code 
   sigPG <- df |>
          filter(.data$adjPval <= params$adjPval_thr) %>%
          arrange(.data$pval) %>% 
          pull(.data$Protein.Group)
   
  #  POIs <-df %>%
  #         dplyr::filter(! is.na(.data$adjPval), .data$adjPval <= params$adjPval_thr) %>%
  #         dplyr::pull(.data$Accession) %>%
  #         unique()
    nPOI <- function(x) sapply(seq_along(x), function(i,x){length(unique(x[1:i]))},x=x)
    POI_Plot <- sigPG[nPOI(sigPG) <=10] |> unique()

    POI_Plot_Genes <- df %>%
      select(.data$Protein.Group, .data$Genes) %>%
      filter(Protein.Group %in% POI_Plot) %>%
      pull(.data$Genes) %>% 
        unique() 
    
    adjAlpha <-  params$adjPval_thr * mean(df$adjPval<= params$adjPval_thr, na.rm = TRUE)
  
  # #
  #  ggrepel::geom_text_repel(data = . %>% filter(interest != "Not Relevant"),
  #                              aes(label = Genes,
  #                                  colour = interest),
  #                              size = 2,
  #                              segment.size = 0.25,
  #                              show.legend = F) +
  
   # dataset 


  levels_interest <- c(POI_Plot_Genes, "Relevant", "Not Relevant")

  # 2. Build the Color Palette Vector
  # We use named vectors: this ensures "Relevant" is ALWAYS red
  poi_cols <- pal_npg()(length(POI_Plot))
  names(poi_cols) <- POI_Plot_Genes

  # Combine into a single named vector
  my_colors <- c(poi_cols, "Relevant" = "red", "Not Relevant" = "black")
  # 3. Build the Alpha Vector
  my_alphas <- c(setNames(rep(1, length(POI_Plot)), POI_Plot_Genes), 
                "Relevant" = 0.2, 
                "Not Relevant" = 0.1)
  df <- df %>%  filter(!is.na(.data$adjPval)) %>% 
  mutate(relevance = ifelse(.data$adjPval<= params$adjPval_thr & abs(.data$logFC)>=params$FC_thr, "Relevant", "Not Relevant"),
         interest = ifelse(.data$Protein.Group %in% POI_Plot, .data$Genes, .data$relevance) %>%
           factor(levels=c(POI_Plot_Genes,"Relevant","Not Relevant"))) %>%
    arrange(desc(.data$interest)) 

  volcano_out <- ggplot(df, 
          aes(x = logFC,
              y = -log10(pval))) +
      geom_hline(yintercept = -log10(adjAlpha)) +
      geom_vline(xintercept = c(-1, 1)) +
      geom_point(size = 1,
                 aes(colour = interest,
                     alpha = interest,
                     shape = pep_type)) +
      # Use a named vector or direct values instead of the .data subset
      scale_colour_manual(values = my_colors) + 
      scale_alpha_manual(values = my_alphas) +
      scale_shape_manual(values = c(15:18),
                     guide = guide_legend(override.aes = list(size = 3))) + 
      scale_x_continuous(breaks = seq(-10, 10, by = 2)) + # seq is safer than -100:100*2
      labs(title = paste0("Differentially ", ifelse(usage, "used", "abundant"), " LiP precursors"),
           x = expression(Log[2](Fold~change)),
           y = expression(-Log[10](P~value)),
           colour = "Interest",
           alpha = "Interest",
           shape = "Peptide type") +
      theme_bw()
  
  return(list ( volcano = volcano_out,  POI_l = POI_Plot ))

}



#' @author Andrea Argentini
#' @title dep_volcano_barcode
#' @description Generates volcano plots and returns both raw and annotated differential expression results for a given contrast. Assumes model results are stored in the specified layer (e.g., "proteinRS", "peptideNorm"); selects proteins of interest (POIs), annotates significance and interest, prepares interactive and static volcano plots, and produces barcode plots for POIs.
#' @param label Character. Contrast label identifying the result column in rowData (e.g., "condB-condA")
#' @param data QFeatures object containing proteomics results and rowData with DE results
#' @param params List. Analysis parameters; expected entries include numeric FC_thr, adjPval_thr and a character vector poi (optional)
#' @param layer Character. Name of the assay/layer in the QFeatures object where rowData contains the result table
#' @param df_anno Data.frame. Annotated Lip/experiment dataframe used to produce barcode plots
#' @return A list with the following elements:
#'   \item{toptable}{data.frame; differential expression results filtered to non-NA adjPval}
#'   \item{volcano}{ggplot object; volcano plot prepared for interactive display (tooltip text added when to_save = FALSE)}
#'   \item{volcano2file}{ggplot object; static/annotated volcano plot suitable for saving/export (to_save = TRUE)}
#'   \item{barcode_plot}{list; barcode plot objects (one per POI)}
#'   \item{POI}{character vector; proteins of interest used for plotting (may be NULL)}
#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate case_when filter left_join arrange desc pull select case_match slice
#' @importFrom stringr str_subset


dep_volcano_barcode <- function ( label, data  ,params,layer , df_anno ){
  cmp = label

  ## get data
  res <- rowData(data[[layer]])[[label]] %>% rownames_to_column(var = "precursor.Id")
  data_df <- as.data.frame(rowData(data[[layer]])) %>% rownames_to_column("precursor.Id")
  temp <- data_df %>% dplyr::select(precursor.Id,  Protein.Group, Accession, Genes,pep_type , total_repeats, start, end )
  all_res <- res %>% left_join(temp, by = "precursor.Id")
  
  all_res__ <- all_res %>% 
      group_by(precursor.Id) %>% 
      dplyr::slice(rep(1, str_count(Protein.Group, ";")+1)) %>%
      mutate(Accession = unlist(str_split(Protein.Group[1], ";"))) %>% 
      ungroup()
  
  
  #  mutate(Proteotypic = ifelse(str_detect(Protein.Group, ";"), 0, 1)) %>%  I think we do not use it 
  #" select POI
  if  (    all(!(params$poi == '')  &  (length( params$poi) >= 1)) ){
    POIs <- params$poi
    if ( !all(params$poi %in% all_res$Accession)){ 
      stop('Uniprot IDs not recognized or not detected in your result ! ')
    } 

  }else{
    POIs <-all_res %>%
          dplyr::filter(!is.na(.data$adjPval), .data$adjPval <= params$adjPval_thr) %>%
          dplyr::pull(.data$Accession) %>%
          unique()
  }
 

  if(length(POIs) > 10){POIs <- NULL}
  # add annotation
 
  all_res_f <- all_res__ %>%
    filter(!is.na(adjPval)) %>%
   mutate(significance = ifelse(adjPval<=0.05 & abs(logFC)>=1, "Significant", "Not Significant"),
         interest = ifelse(Accession %in% POIs, Genes, significance) %>% 
           factor(levels = c(all_res__$Genes[match(POIs, all_res__$Accession)], "Significant", "Not Significant")),
         pep_type = case_match(pep_type,
                              "SemiTryptic" ~ "Semi-Tryptic",
                              "Tryptic" ~ "Tryptic",
                              "NonTryptic" ~ "Non-Tryptic") %>% 
           factor(levels = c("Tryptic", "Semi-Tryptic", "Non-Tryptic"))) %>% 
    arrange(desc(interest))
  

  DEall <- all_res_f


  volcano <- make_volcano_plot (
      df = all_res_f,
      params = params,
      title = paste0("Volcano ", cmp ),
      annotation_fields = c("Genes",'pep_type', 'significance', 'start', 'end' ),
      poi_vis= POIs,
      to_save= FALSE
    )
    #perc_field <- rowData(data[['proteinRS']]) %>% colnames() %>%  stringr::str_subset('perc')
   # export table


  ## volcano annotate with gene name


# log_info(paste0(cmp,' preparing annotated volcano plot ...'))

p_toFile <-  make_volcano_plot (
      df = all_res_f,
      params = params,
      title = paste0("Volcano ", cmp ),
      annotation_fields = c("Accession" ),
      poi_vis= POIs,
      to_save= TRUE
    )
  
  #LiP_annotated --> missing 
 ## barcode plot  
barcode <- lapply(POIs, plot_barcode, df_anno, "Treatment", DEall, T )  
names(barcode) <- POIs
barcode_plotly_list <- lapply(barcode, function(el) {
   #el$gg is the ggplot, el$colour_mapping_significance is the mapping
  plotly_from_ggplot(el$gg, el$colour_mapping_significance, title_center = 0.0, top_margin = 80, tooltip = "text")
})
  
return ( list( toptable =DEall , 
              volcano = volcano, 
            volcano2file = p_toFile , 
            barcode_gg = barcode , 
            barcode_plty = barcode_plotly_list,  POI = POIs ) )

}

#' @author Andrea Argentini
#' @title plot_barcode
#' @description Generate a barcode-style plot for one or more proteins of interest (POIs). The function summarizes precursor-level signals, annotates significance and peptide types, computes coverage and directionality, and delegates plotting to make_barplot.
#' @param POI Character vector of protein accessions (Uniprot IDs) to plot
#' @param input Data frame or tibble with precursor-level annotations and quantitative columns (must include Accession, Precursor.Id, Drug, normPQ, Proteotypic, total_repeats, coverage, Stripped.Sequence, repeat_nr)
#' @param group_column Unquoted column name (tidy evaluated) indicating the grouping variable in input (e.g., Drug)
#' @param DE_result Data frame of differential expression results (must contain precursor.Id, logFC, pval, adjPval)
#' @param indicate_direction Logical; if TRUE annotate significance with directionality (default FALSE)
#' @return A plotting object (the result returned by make_barplot), typically a ggplot or list containing plot elements for the barcode visualization
#' @importFrom dplyr filter distinct group_by mutate left_join ungroup arrange pull
#' @importFrom tidyr uncount
#' @importFrom stats setNames



plot_barcode <- function(POI, se, group_column= 'Treatment', DE_result, 
                              indicate_direction = F, 
                              prot_seq,
                            usage=TRUE, directionality = FALSE, expand = FALSE, params){
  

grouping <- colData(se)[[group_column]]
groups <- unique(grouping)

signif_names <- c(
  "Not Significant",
  paste(groups[1], "Missing"), # verplaatst zodat significante bovenop komen
  paste(groups[2], "Missing"),
  paste(groups[1], "Up"),
  paste(groups[2], "Up"),
  "Missing",
  "Significant"
)

signif_colours <- c(
  "grey40",
  "yellow2", # ook mee verplaatst door hierboven
  "lightskyblue",
  "lightslateblue",
  "orange2",
  "blue",
  "orange"
)  
 
colour_mapping_significance <- setNames(signif_colours, signif_names)

colour_mapping_type <- c("Non-proteotypic, internally repeating" = "green", 
                    "Non-proteotypic" = "red", 
                    "Internally repeating" = "purple", 
                    "Proteotypic" = "grey80")  

  if (usage) 
  { DE_result <- DE_result %>% select ( .data$Precursor.Id, .data$usage_adjPval,.data$usage_pval, .data$usage_logFC, .data$Protein.Group, .data$Genes, .data$length ,.data$pep_type,.data$Stripped.Sequence, .data$Proteotypic)  %>% mutate(pval = usage_pval, 
         adjPval = usage_adjPval,
         effectSize = usage)
  } else {
    DE_result <-  DE_result <- DE_result %>% select(.data$adjPval,.data$pval, .data$logFC, .data$Protein.Group, .data$Genes ,.data$pep_type) %>% mutate(effectSize = logFC)
  }
## to be checked 
 DE_result <- DE_result %>%
    filter(grepl(POI, Protein.Group)) %>%
    select(.data$pval, .data$adjPval,.data$effectSize,.data$Stripped.Sequence, .data$Protein.Group, .data$Proteotypic, .data$Precursor.Id, .data$length, .data$Genes) %>%
   left_join(prot_seq %>% select(Accession, Protein.Sequence), join_by(Protein.Group == Accession) ) %>% 
    mutate( total_repeats = str_count(Protein.Sequence, Stripped.Sequence)) %>%
    uncount(total_repeats, .remove=FALSE) %>%
    group_by(Stripped.Sequence)%>%
    mutate(repeat_nr = 1:max(total_repeats),
           start = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 1],
           end = str_locate_all(Protein.Sequence[1], Stripped.Sequence[1])[[1]][repeat_nr, 2],
           pep_type = classify_trypticity(
             peptide = Stripped.Sequence[1], 
             protein = Protein.Sequence[1], 
             start_pos = start)) %>%
    mutate(tier = match(Precursor.Id, sort(unique(Precursor.Id))), 
         max_tier = max(tier)) %>%
    ungroup() %>%
    mutate(
    type = ifelse(Proteotypic==0 & total_repeats>1, 
                  "Non-proteotypic, internally repeating",
                  ifelse(total_repeats>1, 
                         "Internally repeating",
                         ifelse(Proteotypic==0, 
                                "Non-proteotypic", 
                                "Proteotypic"))) %>%
      factor(levels = c("Non-proteotypic, internally repeating", 
                        "Non-proteotypic",
                        "Internally repeating",
                        "Proteotypic"))
    ) # 
  
completeness <- assay(se)[DE_result$Precursor.Id,] %>%
  is.na() %*% model.matrix(~0+grouping) 
DE_result <- DE_result %>% 
  mutate(
    directionality = directionality,
    missing = 
      (completeness == matrix(
        table(grouping), 
        byrow=TRUE,
        nrow=nrow(completeness),
        ncol=ncol(completeness))) |> 
      apply(1, function(x){
        h <- which(x) |> names() |> paste0() |> gsub(pattern="grouping",replacement="")
        return(ifelse(sum(x)==0, NA, h))
        }), 
  significance = 
    ifelse(
    !is.na(missing),
    ifelse(directionality,
           paste0(missing," Missing"),
           "Missing"),
    ifelse(
      adjPval < params$adjPval_thr & !is.na(adjPval),
      ifelse(directionality,
             ifelse(effectSize < 0 & !is.na(effectSize),
                    paste0(groups[1]," Up"),
                    paste0(groups[2]," Up")),
             "Significant"),
      "Not Significant")) |> factor(levels=signif_names),
    section = ifelse(significance == "Not Significant",
                     "Not Significant",
                          ifelse(!is.na(missing),
                                 "Missing",
                                 "Significant")) |> 
  factor(levels = c("Not Significant", "Significant", "Missing"))
  ) %>% 
  filter(!(significance=="Not Significant"&is.na(adjPval))) %>%
 arrange(significance, type)  

  barplot_obj <- make_barplot( df = DE_result, POI_ = POI, colour_mapping_significance, colour_mapping_type, expand)
 
  return (list(
        gg = barplot_obj
        ##colour_mapping_significance = colour_mapping_significance,
        #colour_mapping_type = colour_mapping_type
      ) )
}

#' @author Andrea Argentini
#' @title make_barplot
#' @description ggplot2 code for the barcode plot. Draws peptide rectangles along protein sequence coordinates, colours by significance and proteotypicity, and returns a ggplot object ready for display or saving.
#' @param df Data frame prepared for plotting (must contain start, end, tier, max_tier, length, coverage, significance, type)
#' @param POI_ Character; protein accession or identifier used in the plot title
#' @param colour_mapping_significance Named character vector mapping significance categories to fill colours
#' @param colour_mapping_type Named character vector mapping peptide type categories to outline colours
#' @return A ggplot object representing the barcode plot
#' @importFrom ggplot2 ggplot aes geom_rect scale_x_continuous scale_y_continuous labs scale_fill_manual scale_colour_manual theme_bw theme element_blank element_rect element_text guide_legend

# colour = type
make_barplot <- function ( df, POI_ , colour_mapping_significance, colour_mapping_type,expand ){
 
   plot <- 
  ggplot(df, aes(x=start, y = 1)) +
    geom_rect(aes(xmin = start,
                  xmax = end,
                  ymin = (1/max_tier)*(tier-1),
                  ymax = (1/max_tier)*(tier),
                  fill = significance,
                  colour = type
                  ),
              linewidth = 0.2) +
    scale_x_continuous(breaks = c(0:1000*10^(floor(log10(df$length[1]-1)))), 
                       limits = c(0,df$length[1]),
                       expand = c(0,0)
                       ) +
    scale_y_continuous(expand = c(0,0)) +
    labs(title = "Significant changes by precursor",
         subtitle = paste(df$Genes[1],"/",POI_,": ",
                          round(calculate_coverage(
                            df$start,
                            df$end,
                            df$length[1])*100,
                            1), "% coverage", sep = ""),
         x = "Residue",
         y = "Precursors",
         fill = "Significance",
         colour = "Proteotypicity") +
    scale_fill_manual(values = colour_mapping_significance) +
    scale_colour_manual(values = colour_mapping_type, 
                        guide = guide_legend(override.aes = list(fill = "transparent"))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          # axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.title.y = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))

if (expand) return(plot + facet_grid(pep_type~section)) else return(plot)
}


#' @author Andrea Argentini
#' @title Convert ggplot to interactive plotly with grouped legends
#' @description
#' Convert a ggplot2 object to an interactive plotly object (via ggplotly),
#' remap trace names to significance group labels (based on the names of
#' colour_mapping_significance) and assign legend groups so traces that
#' represent the same significance category appear as a single legend entry.
#' Only the first trace in each legend group is shown in the legend. The
#' function also preserves ggplot title and subtitle in the plotly layout and
#' allows adjusting title horizontal alignment and top margin. The ggplotly
#' tooltip argument can be controlled with the `tooltip` parameter.
#'
#' @param ggp A ggplot object to convert (created with ggplot2).
#' @param colour_mapping_significance A named character vector whose names are
#'   significance category labels to match against ggplotly trace names. Only
#'   the names are used for matching/legend grouping; values (colours) are not
#'   used by this function but are commonly provided for downstream styling.
#' @param title_center Numeric between 0 and 1 specifying the horizontal
#'   position of the title in the plotly layout (0 = left, 0.5 = center,
#'   1 = right). Default is 0.0 (left).
#' @param top_margin Integer top margin in pixels to apply in the plotly
#'   layout. Default is 70.
#' @param tooltip Character scalar passed to ggplotly's tooltip argument
#'   (e.g. "text", "x", "y", or a vector). Default is "text".
#'
#' @return A plotly object (list-like) produced by ggplotly with adjusted
#'   trace legend groups and layout title/margin.
#' @importFrom plotly ggplotly
#' @importFrom magrittr %>%
plotly_from_ggplot <- function(ggp, colour_mapping_significance, title_center = 0.0, top_margin = 70, tooltip = "text") {
  pp <- ggplotly(ggp, tooltip = tooltip)

  sig_names <- names(colour_mapping_significance)

  # Map ggplotly trace names to significance labels (if matched) and set legendgroup
  for (i in seq_along(pp$x$data)) {
    tr <- pp$x$data[[i]]
    nm <- if (!is.null(tr$name)) tr$name else ""
    matched <- sig_names[vapply(sig_names, function(s) grepl(s, nm, fixed = TRUE), logical(1))]
    if (length(matched) >= 1) {
      matched <- matched[1]
      pp$x$data[[i]]$legendgroup <- matched
      pp$x$data[[i]]$name <- matched
    }
  }

  # show only the first trace per legendgroup
  seen <- character(0)
  for (i in seq_along(pp$x$data)) {
    lg <- pp$x$data[[i]]$legendgroup
    if (!is.null(lg)) {
      if (lg %in% seen) {
        pp$x$data[[i]]$showlegend <- FALSE
      } else {
        pp$x$data[[i]]$showlegend <- TRUE
        seen <- c(seen, lg)
      }
    }
  }

  # ensure layout title + margin; use ggplot labels if present
  gtitle <- if (!is.null(ggp$labels$title)) ggp$labels$title else NULL
  gsubtitle <- if (!is.null(ggp$labels$subtitle)) ggp$labels$subtitle else NULL
  title_text <- if (!is.null(gtitle) || !is.null(gsubtitle)) {
    paste0(
      if (!is.null(gtitle)) paste0("<b>", gtitle, "</b>") else "",
      if (!is.null(gsubtitle)) paste0("<br>", gsubtitle) else ""
    )
  } else NULL

  pp <- pp %>% layout(
    title = list(text = if (!is.null(title_text)) title_text else "", x = title_center),
    margin = list(t = top_margin)
  )

  pp
}

#' @author Andrea Argentini
#' @title check_dependencies
#' @description Check for required R packages and install them if missing. Loads each package after ensuring installation.
#' @param required_packages Character vector of package names to check and (if needed) install
#' @return Invisibly returns NULL. Side effects: installs and loads requested packages.
#' @importFrom utils install.packages
#' @importFrom BiocManager install
#' @export
check_dependencies = function(required_packages = required_packages){
  suppressPackageStartupMessages(
    for(i in required_packages){
      # require returns TRUE invisibly if it was able to load package
      if(! require(i, character.only = TRUE, quietly = TRUE)){
        #  If package was not able to be loaded then re-install
        tryCatch(install.packages(i , dependencies = TRUE), error = function(e) { NULL })
        tryCatch(BiocManager::install(i), error = function(e) { NULL })
        require(i, character.only = TRUE, quietly = TRUE)
      }
    }
  )

}


#' @author Andrea Argentini
#' @title render_child
#' @description Render an Rmd/Quarto child template into the main document. Reads the template file, evaluates it with a small environment containing provided objects, and writes the rendered text to stdout.
#' @param data Data object to expose to the child template (e.g., DE result for a contrast)
#' @param path Character; path where results or outputs should be stored (exposed to the template)
#' @param pe QFeatures or other object that may be required by the template
#' @param label Character; layer/label name passed to the template
#' @param template Path to the template file to render (character)
#' @return Invisibly returns NULL. Side effects: writes rendered child output to stdout.
#' @export
render_child <- function(data, path, pe, label ,  template) {
  
    # _templateBArPlot.Rmd _templateContrast _templatePval
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe,  label = label,  path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  
}


#' @author Andrea Argentini
#' @title flag_complex_formula
#' @description Inspect a linear model formula and detect disallowed complex terms.
#' Complex terms include interaction operators (":" or "*") and function/transformation
#' calls (e.g., log(A), I(A^2), poly(A,2), scale(A)). Only simple main-effect
#' variable names (letters, numbers, underscores, dots) are considered valid.
#' @param formula A model formula (object of class "formula") to check, e.g. ~ Group + Age
#' @return A list with elements:
#'   \item{flag}{logical; TRUE if disallowed (complex) terms are present, FALSE otherwise}
#'   \item{problematic_terms}{character; vector of the disallowed term labels found, or NULL if none}
#' @importFrom stats terms

flag_complex_formula <- function(formula) {
  # Extract terms
  terms_obj <- terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  
  # Allowed = simple variable names (letters, numbers, underscores, dots)
  allowed_pattern <- "^[A-Za-z0-9_.]+$"
  
  # Flag terms that are not simple names
  bad_terms <- term_labels[!grepl(allowed_pattern, term_labels)]
  
  if (length(bad_terms) > 0) {
    return(list(flag = TRUE, problematic_terms = bad_terms))
  } else {
    return(list(flag = FALSE, problematic_terms = NULL))
  }
}