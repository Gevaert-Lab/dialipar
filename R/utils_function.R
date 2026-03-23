
#' @title Create MDS Plot Grob
#' @description 
#' Performs Multi-Dimensional Scaling (MDS) on joined assays within a QFeatures 
#' object, calculates variance explained, and returns a rendered ggplot grob.
#' 
#' @param pe A \code{QFeatures} object containing \code{precursors_lip_norm} 
#' and \code{precursors_tc_norm} assays.
#' 
#' @return A \code{gtable} (grob) object representing the MDS plot.
#' 
#' @importFrom QFeatures joinAssays
#' @importFrom MultiAssayExperiment getWithColData
#' @importFrom scater runMDS plotMDS
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw ggplotGrob
#' @importFrom methods as
#' @importFrom rlang .data

make_mds <- function(pe) {
  # 1. Join assays and extract as SummarizedExperiment
  se <- QFeatures::joinAssays(pe, i = c("precursors_lip_norm", "precursors_tc_norm"), fcol = "Precursor.Id") |> 
    MultiAssayExperiment::getWithColData("joinedAssay")

  # 2. Convert to SingleCellExperiment (required by scater) and run MDS
  sce <- as(se, "SingleCellExperiment")
  sce <- scater::runMDS(sce, exprs_values = 1)

  # 3. Extract Eigenvalues for axis labels
  eig <- attr(SingleCellExperiment::reducedDim(sce, "MDS"), "eig")

  if (!is.null(eig)) {
    percent_var <- (eig / sum(eig)) * 100
    xlbl <- paste0("Component 1 (", round(percent_var[1], 1), "%)")
    ylbl <- paste0("Component 2 (", round(percent_var[2], 1), "%)")
  } else {
    xlbl <- "MDS1"
    ylbl <- "MDS2"
  }

  # 4. Internal Helper to strip data environments
  as_lean_grob <- function(p) {
    if (is.null(p)) return(NULL)
    return(ggplot2::ggplotGrob(p))
  }

  # 5. Create Plot
  # Note: scater::plotMDS creates the base ggplot

  mds_plot <- scater::plotMDS(sce, colour_by = "Condition")
  
# 2. Add the 'text' aesthetic to the existing mapping for plotly
  # This avoids adding a 2nd layer and keeps the scater internal data
  mds_plot$layers[[1]]$mapping <- ggplot2::aes(
    color = colData(se)$Condition,
    text = paste0(
      "Sample: ", colData(se)$CondRep, "<br>",
      "Condition: ", colData(se)$Condition, "<br>",
      "Replicate: ", colData(se)$Replicate
    )
  )

  mds_plot <- mds_plot + 
    ggplot2::labs(x = xlbl, y = ylbl) +
    ggplot2::theme_bw()

  return(as_lean_grob(mds_plot))
}



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
#' @importFrom utils object.size
#' @importFrom rlang .data
#' @importFrom lobstr obj_size

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
  log_info(sprintf("Size INSIDE RENDER  : %.2f MB", as.numeric(obj_size(data_list)) / 1024^2))
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

#' @title process_dialipa_data
#' 
#' @description
#' Main processing pipeline for DIA-LiPA data. Orchestrates input parsing, 
#' normalization, precursor annotation, and differential expression modeling 
#' for both paired and unpaired experimental designs.
#'
#' @param params_report A list containing report parameters including file paths, 
#' comparisons, and thresholds.
#' @param analysis_type Character string. Either "unpaired" (default) or "paired".
#'
#' @return A list with three elements: `error` (string), `status` (integer, 0 for success), 
#' and `result` (a list containing QC data, MDS plots, and DE results).
#' 
#' @export
#' 
#' @importFrom logger log_info
#' @importFrom rlang .data
#' @importFrom utils object.size
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom S4Vectors colnames
#' @importFrom lobstr obj_size

process_dialipa_data <- function (params_report, analysis_type = "unpaired" ){
  fastaproc <- read_fasta_ann(params_report$fasta_file )
      input_data <- parse_input(params_report$input_file_tc, 
                                params_report$input_file_lip,   
                                params_report$design_file)

        qf_base <- create_qfeat_(input_data$design, input_data$lip, input_data$tc, input_data$diann_flag)
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
        ## qc_ann$result go to data bag.
        if (analysis_type == 'paired') {
              log_info('Paired branch ...')
              res_de_norm <-  msqrob_model(pe = qf_final, params = params_report, layer = 'precursors_lip_norm' )
              if (res_de_norm$status == 1) stop(res_de_norm$error)
              res_de_usage <-  msqrob_model(pe = res_de_norm$q_feat, params = params_report, layer = 'precursors_lip_usage' )
              if (res_de_usage$status == 1) stop(res_de_usage$error)
              df_ann <-  as.data.frame(SummarizedExperiment::rowData(res_de_usage$q_feat[["precursors_lip_norm"]])[c("Precursor.Id", "Protein.Group", "Genes", "Proteotypic", "Stripped.Sequence")])
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
                                           contrasts = base::colnames(b$contr_exp))
                if (qf_unpair$status == 1) stop(qf_unpair$error)  
                ## remove all :
                assays_to_keep <- c("precursors_lip_norm", "precursors_tc_norm" )
                # Subset to only these
                qf_final <- qf_unpair$qf[, , assays_to_keep]         
                df_ann <-  as.data.frame(rowData(qf_final[["precursors_lip_norm"]])[c("Precursor.Id", "Protein.Group", "Genes", "Proteotypic", "Stripped.Sequence")])
                
          log_info(sprintf("Size INPUT qf_unpair: %.2f MB", as.numeric(obj_size(qf_final)) / 1024^2))
                log_info(sprintf("Size INPUT df_ann: %.2f MB", as.numeric(obj_size(df_ann)) / 1024^2))
                test_ <-  lapply(params_report$comparisons, 
                    build_df_result,
                    data= qf_final  ,
                    df_anno = df_ann  ,
                    mapping_df = fastaproc$result,
                    layer= 'precursors_lip_norm' ,
                    layer_ = NULL,
                  params = params_report)
                names(test_)<- params_report$comparison_label
                # log_info('Size from the returned list  ')
                # log_info(sprintf("Size ggplot: %.2f MB", as.numeric(obj_size(test_$`Rapa - Dmso`$plotvolcano_usage)) / 1024^2))
                #qf_unpair$qf   
                
                mds_plot <- make_mds( pe = qf_final  )
          
                log_info('Size mds  ')
                log_info(sprintf("Size ggplot: %.2f MB", as.numeric(obj_size(mds_plot)) / 1024^2))
                

        }
        ##  qc -> data for QC
        ##  pe -> Qfeat subseted
       
        quarto_bag <- list( qc_data = qc_ann$result,
                            mds = mds_plot,
                            res_DE= test_  )

        return( list(error= '', status= 0,result =quarto_bag ))
 
  
}

#' @title Build Result Data Frame for a Single Contrast
#'
#' @description
#' Internal helper function to extract, format, and standardize differential expression 
#' results for a specific contrast. It handles both "paired" and "unpaired" 
#' workflows by unifying column names and generating associated plots.
#'
#' @param label Character. The name of the contrast (e.g., "TreatmentA - TreatmentB").
#' @param data A QFeatures object.
#' @param layer Character. Main assay name (normalized LiP).
#' @param layer_ Character or NULL. Secondary assay name (Usage).
#' @param df_anno data.frame. Precursor metadata for joining.
#' @param mapping_df data.frame. FASTA/Protein mapping metadata.
#' @param params list. Analysis parameters.
#'
#' @return A list containing the results table and serialized/lean plot objects.
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom QFeatures joinAssays
#' @importFrom dplyr rename rename_with left_join mutate full_join relocate select filter join_by contains
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove fixed
#' @importFrom ggplot2 ggplotGrob
#' @importFrom rlang .data
#' @keywords internal

build_df_result <- function (label, data , layer , layer_ = NULL , df_anno, mapping_df, params){
 # --- 1. Process Assay A (e.g., Lip normalized) ---
      res_layer <-  rowData(data[[layer]])[[label]]
      res_layer_df <- as.data.frame(res_layer, check.names = FALSE) %>% 
      rownames_to_column(var = "Precursor.Id") 
      #res_layer_df <- res_layer_df[, !sapply(res_layer_df, is.list)]
# --- 2. Process Assay B (Optional, e.g., Protein/Norm) ---
  if (!is.null(layer_)) {
     ## usage for paired 
      res_layer_ <- rowData(data[[layer_]])[[label]]
    
      res_layer__df <- as.data.frame(res_layer_, check.names = FALSE) 
    
      #res_layer__df <- res_layer__df[, !sapply(res_layer__df, is.list)] # Drop models
      res_layer__df <-  res_layer__df %>%   rownames_to_column(var = "Precursor.Id") %>% 
                                      rename_with(~paste0("usage_", str_remove(.x, fixed(paste0(label, ".")))), .cols = -.data$Precursor.Id)
    # Merge A and B
    res_combined <- full_join(res_layer_df, res_layer__df, by = "Precursor.Id")  
  }else {
    # === UNPAIRED BRANCH ===
    # Usage is already in res_layer_df, but with suffixes (_usage).
    # We rename them to match the paired prefixes (usage_).
    res_combined <- res_layer_df %>%
      dplyr::rename(
        usage_logFC   = .data$usage,            # 'usage' is the logFC
        usage_pval    = .data$pval_usage,
        usage_adjPval = .data$adjPval_usage,
        usage_se      = .data$se_usage,
        usage_t       = .data$t_usage,
        usage_df      = .data$df_usage
      )
  }
  # --- 3. Final Join with Annotation ---
  final_df <- res_combined %>% 
    left_join(df_anno, by = "Precursor.Id") %>%
    mutate(contrast = .env$label) # Good for downstream filtering
  
  # --4 join with sequence 
  final_df <- final_df %>% 
    left_join(mapping_df, by = join_by ( Protein.Group == Accession)) %>%
    pep_char()
   
  # top table corrected and not 
  full_toptable <- final_df %>% relocate("Precursor.Id", contains("usage"))
  
  ## filtering null adj pval 

  full_toptable <- full_toptable %>% dplyr::filter(!is.na(.data$usage_adjPval) & (!is.na(.data$adjPval)))
  full_toptable <- full_toptable %>% 
                                  dplyr::select(
                                    .data$Precursor.Id, .data$usage_adjPval, .data$usage_df, .data$usage_logFC, 
                                    .data$usage_pval, .data$usage_se, .data$usage_t, .data$Protein.Group, 
                                    .data$Genes, .data$Proteotypic, .data$Stripped.Sequence, .data$contrast, 
                                    .data$Protein.Sequence, .data$length, .data$missed_cleavages, 
                                    .data$total_repeats, .data$start, .data$end, .data$pep_type, .data$AA_last
                                  ) %>% 
                                  dplyr::filter(!is.na(.data$usage_adjPval))

   # Define a NEW "Grob Stripper"
  as_lean_grob <- function(p) {
    if (is.null(p)) return(NULL)
    
    # 1. Convert ggplot to a gtable/grob (this 'renders' the data into coordinates)
    g <- ggplot2::ggplotGrob(p)
    
    # 2. This object is now just 'points and lines'—the 160MB dataframe is GONE.
    return(g)
  }


  clean_plotly <- function(p) {
  # 1. Remove the original ggplot object hidden inside
  p$x$visdat <- NULL 
  
 # 2. Scrub environments from any remaining attributes
  # This targets internal plotly metadata that might still have pointers
  if (!is.null(p$x$attrs)) {
    p$x$attrs <- lapply(p$x$attrs, function(x) {
      attr(x, ".Environment") <- NULL
      return(x)
    })
  }

  # 3. Handle the data list
  # Plotly stores the processed data here; we just want the raw values, not the env
  p$x$data <- lapply(p$x$data, function(d) {
    attr(d, ".Environment") <- NULL
    return(d)
  })
  
  # 3. Nuclear Strip (Serialization Reboot)
  p <- unserialize(serialize(p, NULL))
  
  return(p)
}  
    res_usage  <- make_volcano_plot_plotly (
      df = full_toptable,
      params = params,
      title = paste0("Volcano ", 'TEST' ),
      to_save= FALSE
    )
      res_usage_plt <-  res_usage  # clean_plotly(res_usage$volcano)
     
     #res_usage_grob <- as_lean_grob(res_usage$volcano)
    
     res_correct_file   <- make_volcano_plot (
       df = full_toptable,
       params = params,
       title = paste0("Volcano ", 'TEST' ),
       usage = TRUE,
      to_save= TRUE
     )
  
    res_ncorrect_grob <- as_lean_grob(res_correct_file$volcano)

  prot_seq_ <-  mapping_df %>% dplyr::filter(.data$Accession %in% res_usage$POI_l)
  
  se <- joinAssays(data,i=c("precursors_lip_norm","precursors_tc_norm"),  fcol = "Precursor.Id") %>%  
          MultiAssayExperiment::getWithColData("joinedAssay")
  
  
   barplot_id <- plot_barcode (res_usage$POI_l, se , 
                  DE_result = full_toptable,   prot_seq = prot_seq_ ,  
                group_column = 'Treatment'  , params= params )
  
  barplot_id  <-  as_lean_grob( barplot_id$gg)

  #plotvolcano_ncorr = res_ncorrect_grob,
 

    return(  list( full_toptable    =    full_toptable,
                   plotvolcano_usage  =  res_usage_plt,
                   plotvolcano_file  =  res_ncorrect_grob,
                   barplot_= barplot_id
                  )
                )

}


#' @title Calculate LiP-MS Usage (Corrected logFC)
#'
#' @description
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
#' @return A \code{list} with \code{status}, \code{error}, and the updated \code{QFeatures} object.
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom stats pt p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom logger log_info

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
    base::colnames(res_tc) <- paste0(base::colnames(res_tc), "_tc")
    
    # 2.2. Align TC results to the LiP precursors
    # We subset the TC results using the match index
    log_info('Align TC results to the LiP precursors  ...')
    aligned_tc <- res_tc[match_idx, , drop = FALSE]
    
    # 2.3. Combine and calculate usage
    # We convert to a standard data frame temporarily for easier calculation
    res_lip <- SummarizedExperiment::rowData(qf[[i_lip]])[[contrast]]
    combined <- BiocGenerics::cbind(res_lip, aligned_tc)
    
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
  }, error = function(err) {
        msg <- conditionMessage(err) # <--- Using 'err' here
        logger::log_error("LiP usage error: {msg}")
        return( list(error= err, status= 1, qf =NULL ))

  })
}

##-----------------

#' @title Check Design File Requirements
#' 
#' @description 
#' Performs critical validation checks on the experimental design data frame:
#' \enumerate{
#'   \item Verifies all required columns are present.
#'   \item Ensures the 'Run' column does not contain file extensions (e.g., .raw, .mzML).
#'   \item Validates that the 'Pipeline' column only contains "LiP" or "TC".
#' }
#'
#' @param df A \code{data.frame} containing the experimental design.
#' @param required_cols A \code{character} vector of column names that must exist in \code{df}.
#' 
#' @return A \code{list} with two elements:
#' \itemize{
#'   \item \code{status}: Integer; 0 for success, 1 if a validation error was found.
#'   \item \code{error}: Character; a descriptive error message if \code{status} is 1.
#' }
#' 
#' @author Andrea Argentini
#' 
#' 
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
  if ("Run" %in% base::colnames(df)) {
    # Added trimws to handle accidental spaces before checking extensions
    run_vals <- trimws(as.character(df$Run))
    bad_runs <- grep("\\.[A-Za-z0-9]+$", run_vals, value = TRUE)
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
    # unique(df$Pipeline) handles factors or characters correctly
    bad_vals <- base::setdiff(unique(as.character(df$Pipeline)), allowed)
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


#' @title create_qfeat_
#' 
#' @description 
#' Converts raw data frames into a QFeatures object, calculates group-specific 
#' detection heuristics, and performs initial quality filtering.
#'
#' @param annotation_df A data frame containing sample metadata (colData).
#' @param report_file1 A data frame containing quantitative proteomics data.
#' @param report_file2 Optional; secondary report file (default is NULL).
#' @param diann_flag Logical; TRUE if data is from DIA-NN, FALSE for Spectronaut.
#' 
#' @return A list with status, error message, and the resulting QFeatures object.
#'
#' @export
#'
#' @import QFeatures
#' @importFrom SummarizedExperiment assay colData rowData "rowData<-"
#' @importFrom MultiAssayExperiment getWithColData
#' @importFrom dplyr filter bind_rows
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats model.matrix
#' @importFrom matrixStats rowMins
#' @importFrom logger log_info
#' @importFrom BiocGenerics ncol
#'
create_qfeat_ <- function(annotation_df, report_file1, report_file2 = NULL, diann_flag) {

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
    if (! is.null(report_file2) ){
        input_lip <- report_file1 %>% dplyr::filter(.data$Precursor.Quantity > 4)
        input_tc  <- report_file2 %>% dplyr::filter(.data$Precursor.Quantity > 4)
        input_data <- dplyr::bind_rows(input_lip, input_tc)
    }else{
      input_data <- report_file1 %>% dplyr::filter(.data$Precursor.Quantity > 4)
      }
    
    
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
     if (diann_flag ){
    qf <- QFeatures::filterFeatures(qf, ~ Q.Value <= 0.01 & 
                                      PG.Q.Value <= 0.01 & 
                                      Lib.Q.Value <= 0.01 & 
                                      Precursor.Id != "" & 
                                      Decoy == 0)
     }else {
     qf <- QFeatures::filterFeatures(qf, ~ Q.Value <= 0.01 & 
                                       PG.Q.Value <= 0.01 & 
                                       Precursor.Id != "" )
     }

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


#' @title Log-Transformation, Median Normalization, and Protein Aggregation
#' 
#' @description 
#' Performs log2 transformation on LiP and TC assays, calculates sample-based 
#' normalization factors using common features, and aggregates TC precursors 
#' to the protein level.
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
#'
#' @import QFeatures
#' @importFrom SummarizedExperiment assay
#' @importFrom matrixStats colMedians
#' @importFrom MsCoreUtils medianPolish
#' @importFrom logger log_info
#' @importFrom rlang .data
#' @importFrom QFeatures logTransform sweep
#' @importFrom stats median na.exclude

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
    #m_complete <- stats::na.exclude(m)
    m_complete <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
    if (nrow(m_complete) == 0) {
             stop(paste("No common features found in assay", i, "to calculate normalization factors."))
    }
    norm_factors <- matrixStats::colMedians(m_complete)
    
    # 4. Zero-center the normalization factors
    norm_factors <- norm_factors - stats::median(norm_factors)
    
    # 5. Sweep out the factors
    qf <- sweep(
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
      logger::log_info("normalization_scaling_factor Error: {msg}")
      # Return the object in its current state even if error occurs
      return(list(error = msg, status = 1, result = q_feat))
    }
  )
}

####----

#' Compute LiP-MS Usage 
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
#' @importFrom stats formula 
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
        msg <- conditionMessage(err)
        logger::log_info("compute_usage error: {msg}")
        return( list(error= err, status= 1,result = NULL ))
  } )
 
}

#' @title Classify Peptide Trypticity
#' 
#' @description 
#' Categorizes peptides as Tryptic, Semi-Tryptic, or Non-Tryptic based on 
#' Protease (Trypsin) cleavage rules, accounting for protein termini 
#' and N-terminal Methionine excision.
#'
#' @param peptide Character vector of peptide sequences.
#' @param protein Character vector of the parent protein sequences.
#' @param start_pos Numeric vector of the 1-based start position.
#'
#' @return A character vector: "Tryptic", "Semi-Tryptic", "Non-Tryptic", or "Ambiguous".
#'
#' #' @importFrom rlang .data
  
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
  
#' @title Calculate Protein Sequence Coverage
#' 
#' @description 
#' Calculates the fraction of a protein sequence covered by a set 
#' of peptides or fragments. It accounts for overlapping regions by treating 
#' the segments as range-based intervals.
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
#' @author Andrea Argentini
#' @export
#'
#' @importFrom IRanges IRanges coverage
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @importFrom methods as
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

#' @title Characterize Peptide Properties and Mapping
#' 
#' @description 
#' Annotates a data frame of peptides with biochemical and positional properties 
#' such as missed cleavages, mapping positions, and trypticity.
#'
#' @param table A data frame or tibble.
#' @param prot_seq Character. Column name for protein sequence (default "Protein.Sequence").
#' @param pep_seq Character. Column name for stripped peptide sequence (default "Stripped.Sequence").
#' 
#' @return A data frame with additional columns for sequence properties.
#' 

#' @importFrom dplyr mutate select .data
#' @importFrom stringr str_count str_locate str_sub
#' @importFrom magrittr %>%
#'
pep_char <- function(table, prot_seq = "Protein.Sequence", pep_seq = "Stripped.Sequence") {
  
  table <- table %>%
    dplyr::mutate(
      # 1. Missed Cleavages (using .data instead of get() for stability)
      missed_cleavages = stringr::str_count(.data[[pep_seq]], "[RK](?!(P|$))"),
      
      # 2. Total Repeats
      total_repeats = stringr::str_count(.data[[prot_seq]], .data[[pep_seq]]),
      
      # 3. Positional Mapping
      # We extract the matrix columns immediately to avoid keeping the matrix 'tmp'
      start = stringr::str_locate(.data[[prot_seq]], .data[[pep_seq]])[, 1],
      end   = stringr::str_locate(.data[[prot_seq]], .data[[pep_seq]])[, 2],
      
      # 4. Trypticity Classification
      pep_type = ifelse(
        (.data$total_repeats > 1) | is.na(.data$total_repeats),
        "Ambiguous",
        classify_trypticity(
          peptide = .data[[pep_seq]], 
          protein = .data[[prot_seq]], 
          start_pos = .data$start
        )
      ),
      
      # 5. C-terminal Amino Acid
      AA_last = stringr::str_sub(.data[[pep_seq]], -1, -1)
    )
  
  return(table)
}


#' @title Precursor Annotation for Quality Control Plots
#' 
#' @description 
#' Converts QFeatures assays into a long-format data frame and performs 
#' comprehensive annotation including protein mapping, peptide characterization, 
#' and condition formatting.
#'
#' @param q_feat A \code{QFeatures} object.
#' @param mapping A data frame used for joining protein accessions to metadata.
#' @param type Character. Either "paired" or "unpaired".
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{error}: Character string containing error messages.
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{result}: An annotated data frame in long format.
#' }
#'
#' @export
#' @importFrom QFeatures longForm
#' @importFrom dplyr left_join mutate case_match join_by .data
#' @importFrom magrittr %>%
#' @importFrom logger log_info


qc_precursor_annotation <- function(q_feat, mapping, type) {
  tryCatch(expr = { 
    log_info('Annotate precursor for QC plot ...')
    
    # 1. Define layers based on analysis type
    # Using 'intersect' is safer in case an assay name was modified elsewhere
    available_layers <- names(q_feat)
    target_layers <- if (type == 'paired') {
      c("precursors_lip_norm", "precursors_tc_norm", "precursors_lip_usage")
    } else {
      c("precursors_lip_norm", "precursors_tc_norm")
    }
    layer <- intersect(target_layers, available_layers)
    
    # 2. Extract and Join
    qcObj <- q_feat[,,layer] %>%
      QFeatures::longForm(
        colvars = c("Condition", "CondRep", "Treatment", "Replicate", "Pipeline"), 
        rowvars = c("Precursor.Id", "Protein.Group", "Stripped.Sequence", 
                    "Precursor.Charge", "Genes")
      ) %>%
      as.data.frame() %>%
      # Fix: Removed x$ and y$ references for join_by compatibility
dplyr::left_join(mapping, by = dplyr::join_by(Protein.Group == Accession)) %>%      pep_char() %>%
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



#' @title Parse Proteomics Input and Design Files
#' 
#' @description 
#' Reads and parses input parquet file(s) and the experiment design file. 
#' Detects whether reports are in DIA-NN or Spectronaut format, performs 
#' column renaming to a unified standard, and validates the design file.
#'
#' @param input_parquet_tc Path to the TC parquet file (optional).
#' @param input_parquet_lip Path to the LiP parquet file.
#' @param input_design Path to the experiment design file (tsv).
#' 
#' @return A list containing status, error message, parsed reports (lip, tc), 
#'   the design data frame, and a diann_flag.
#' 
#' @author Andrea Argentini
#' 
#' @importFrom arrow read_parquet read_tsv_arrow
#' @importFrom dplyr rename .data
#' @importFrom logger log_info
#' @importFrom magrittr %>%
parse_input <- function(input_parquet_tc, input_parquet_lip, input_design) {
  
  # Helper to check for DIA-NN format
  is_diann <- function(df) { "Run" %in% base::colnames(df) }
  
  # Unified renaming map for Spectronaut (New_Name = Old_Name)
  spec_rename_map <- c(
    Run = "R_FileName",
    Genes = "PG_Genes",
    Protein.Group = "PG_ProteinGroups",
    Protein.Names = "PG_ProteinNames",
    PG.Q.Value = "PG_Qvalue",
    Precursor.Id = "PEP_GroupingKey",
    Stripped.Sequence = "PEP_StrippedSequence",
    Decoy = "EG_IsDecoy",
    Precursor.Charge = "FG_Charge",
    Q.Value = "FG_Qvalue",
    Precursor.Quantity = "FG_MS2RawQuantity"
  )

  tryCatch(
    expr = {
      # 1. Loading Reports
      if (!is.null(input_parquet_tc) && input_parquet_tc != "") {
        logger::log_info('Reading TC and LiP from SEPARATE parquet files ...')
        TC_report <- arrow::read_parquet(input_parquet_tc)
        LiP_report <- arrow::read_parquet(input_parquet_lip)
      } else {
        logger::log_info('Reading both LiP and TC from ONE parquet file ...')
        LiP_report <- arrow::read_parquet(input_parquet_lip)
        TC_report <- NULL
      }
      
      # 2. Detect Format and Standardize
      diann_flag <- is_diann(LiP_report)
      
      if (!diann_flag) {
        logger::log_info('Standardizing Spectronaut format to DIA-NN style...')
        # Use any_of to safely rename only columns that exist
        LiP_report <- LiP_report %>% dplyr::rename(dplyr::any_of(spec_rename_map))
        
        if (!is.null(TC_report)) {
          TC_report <- TC_report %>% dplyr::rename(dplyr::any_of(spec_rename_map))
        }
      }

      # 3. Loading and Validating Design
      logger::log_info('Reading experiment Design file ...')
      design <- arrow::read_tsv_arrow(input_design)
      
      col_design_required <- c('Run', 'Pipeline', 'Treatment', 'Condition', 'Replicate', 'CondRep')
      checkdesign <- check_design_requirement(design, col_design_required)
      
      if (checkdesign$status == 1) {
        return(list(error = checkdesign$error, status = 1, lip = NULL))
      }

      # Rename Run to runCol for QFeatures compatibility
      design <- design %>% dplyr::rename(runCol = .data$Run)
      
      return(list(
        error      = '', 
        status     = 0,
        lip        = LiP_report,
        tc         = TC_report,
        design     = design,
        diann_flag = diann_flag
      ))
    },
    error = function(err) {
      msg <- conditionMessage(err)
      logger::log_error("Input Parquet Error: {msg}")
      return(list(error = msg, status = 1, lip = NULL))
    }
  )
}



#' @title Calculate Sequence Coverage per Protein
#' 
#' @description 
#' Calculates the fraction of each protein sequence covered by the detected 
#' peptides. Overlapping peptide regions are reduced to unique amino acid 
#' positions to ensure accurate coverage mapping.
#'
#' @param report A \code{data.frame} or \code{tibble} containing "Accession", 
#'   "length", "start", and "end".
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{error}: Character string containing error messages.
#'   \item \code{result}: The original data frame with an added "coverage" column.
#' }
#' 
#' @author Andrea Argentini
#'
#' @importFrom dplyr distinct left_join group_by summarise .data
#' @importFrom IRanges IRanges reduce width
#' @importFrom logger log_info
#' @importFrom magrittr %>%

calculate_coverages <- function(report) {

  tryCatch(expr = {
    logger::log_info('Computing Sequence Coverage ...')

    # 1. Calculate unique coverage per Accession
    coverages <- report %>%
      dplyr::distinct(.data$Accession, .data$length, .data$start, .data$end) %>%
      dplyr::group_by(.data$Accession) %>%
      dplyr::summarise(
        coverage = {
          # Use stats::na.omit to prevent IRanges from crashing on missing coords
          s <- stats::na.omit(.data$start)
          e <- stats::na.omit(.data$end)
          
          if (length(s) == 0) {
            0
          } else {
            ranges <- IRanges::IRanges(start = s, end = e)
            # reduce() merges overlapping intervals into single contiguous ranges
            covered_positions <- sum(IRanges::width(IRanges::reduce(ranges)))
            # Accessing the first element of length for this group
            covered_positions / .data$length[1]
          }
        },
        .groups = "drop"
      )

    # 2. Join back to the original report
    report_annotated <- report %>% 
      dplyr::left_join(coverages, by = "Accession")

    return(list(error = '', status = 0, result = report_annotated))

  }, error = function(err) {
    msg <- conditionMessage(err)
    logger::log_error("Coverage computation error: {msg}")
    return(list(error = msg, status = 1, result = NULL))
  })
}




#' @title Read and Annotate FASTA File
#' 
#' @description 
#' Reads a FASTA file and extracts protein accessions, full sequences, and 
#' sequence lengths. It assumes a standard UniProt-style header where the 
#' accession is the second field when split by a pipe character.
#'
#' @param input_fasta Path to the FASTA file (character).
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{error}: Character string containing error messages.
#'   \item \code{result}: A \code{data.frame} with columns: \code{Accession}, 
#'     \code{Protein.Sequence}, and \code{length}.
#' }
#' 
#' @author Andrea Argentini
#'
#' @importFrom seqinr read.fasta
#' @importFrom stringr str_extract word
#' @importFrom stringr str_split_i
#' @importFrom logger log_info
#' @importFrom magrittr %>%

read_fasta_ann <- function(input_fasta) {

  tryCatch(expr = {
    logger::log_info('Reading FASTA file ...')
    
    # Read FASTA; as.string = TRUE returns sequences as single strings
    # seqtype = "AA" specifies Amino Acids
    fasta_list <- seqinr::read.fasta(file = input_fasta, 
                                     seqtype = "AA", 
                                     as.string = TRUE, 
                                     forceDNAtolower = FALSE)
    
    headers <- unlist(lapply(fasta_list, function(x) attr(x, "Annot")))
    
    # Clean the '>' from the start of the annotation if it exists
    headers <- stringr::str_remove(headers, "^>")
    # --- ADD THE NEW LOGIC HERE ---
    # 1. Try to extract a valid UniProt Accession using Regex
    accessions <- stringr::str_extract(headers, "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    
    # 3. Fallback for headers without a standard UniProt ID pattern
    if (any(is.na(accessions))) {
      na_idx <- is.na(accessions)
      
      # Check if header contains pipes
      # If YES: split by pipe and take the 2nd element
      # If NO: take the 2nd word (to skip 'sp')
      accessions[na_idx] <- ifelse(
        grepl("\\|", headers[na_idx]),
        stringr::str_split_i(headers[na_idx], "\\|", 2),
        stringr::word(headers[na_idx], 2) # <--- This fixes your 'sp' issue
      )
    }
    
    # 4. Final fallback: If still NA, just use the whole header
    accessions <- ifelse(is.na(accessions), headers, accessions)
    # --- END OF NEW LOGIC ---

    # Now create the data frame using the prepared 'accessions' vector
    mapping <- data.frame(
      Accession = accessions,
      Protein.Sequence = as.character(unlist(fasta_list, use.names = FALSE)),
      stringsAsFactors = FALSE
    )
    
    # 3. Add sequence length
    mapping$length <- nchar(mapping$Protein.Sequence)
    
    return(list(error = '', status = 0, result = mapping))

  }, error = function(err) {
    msg <- conditionMessage(err)
    logger::log_error("Reading FASTA Error: {msg}")
    return(list(error = msg, status = 1, result = NULL))
  })
}



#' @title Fit msqrob2 Models and Test Contrasts
#' 
#' @description 
#' Fits robust linear models to proteomics data using the msqrob2 framework. 
#' It performs hypothesis testing based on a user-defined formula and 
#' contrast matrix, handling peptide-level or protein-level quantification.
#'
#' @param pe A \code{QFeatures} object.
#' @param params A list containing \code{formula} (string) and \code{comparisons} (character vector).
#' @param layer Character. The name of the assay in \code{pe} to model.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item \code{status}: Integer (0 for success, 1 for error).
#'   \item \code{error}: Character string containing error messages.
#'   \item \code{q_feat}: The updated \code{QFeatures} object with models and test results.
#'   \item \code{contr_exp}: The contrast matrix used for testing.
#' }
#' 
#' @author Andrea Argentini
#' @importFrom SummarizedExperiment rowData
#' @importFrom msqrob2 msqrob getCoef makeContrast hypothesisTest
#' @importFrom stats as.formula
#' @importFrom logger log_info
#' @importFrom methods is


msqrob_model <- function(pe, params, layer  ){

  tryCatch( expr = {

    
   logger::log_info(paste('Fitting msqrob2 model on layer:', layer, '...'))    
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

   logger::log_info('Generating contrast matrix and performing hypothesis testing...')
    L <- makeContrast(contrast_list, parameterNames = coef)
    pe <- hypothesisTest(object = pe, i = layer, contrast = L , overwrite=TRUE)

    return (list(error= '', status= 0,q_feat = pe , contr_exp = L   ))


  },error = function(err){
    print(paste("Msqrob modeling :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

}

#' @title Generate Interactive Volcano Plot with Plotly
#' 
#' @description 
#' Creates a highly customized, interactive volcano plot for LiP-MS data. 
#' The plot features multi-layered legends (Interest and Peptide Type), 
#' RGBA transparency for non-significant points, and automated labeling 
#' of the top 10 most significant protein groups.
#'
#' @param df A data frame containing statistical results (logFC, p-values, adjPval).
#' @param params A list containing threshold parameters:
#'   \itemize{
#'     \item \code{adjPval_thr}: Adjusted P-value threshold (e.g., 0.05).
#'     \item \code{FC_thr}: Log2 Fold Change threshold (e.g., 1).
#'   }
#' @param title Character. The title of the plot.
#' @param annotation_fields Optional character vector of column names for tooltips.
#' @param to_save Logical. If TRUE, skips tooltip generation to reduce object size.
#' @param usage Logical. If TRUE, uses "usage" specific columns (e.g., usage_logFC).
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{volcano}: A plotly object.
#'   \item \code{POI_l}: A character vector of Protein Groups of Interest.
#' }
#' 
#' @export
#'
#' @import plotly
#' @importFrom dplyr mutate select filter arrange pull sample_frac bind_rows .data any_of
#' @importFrom ggsci pal_npg
#' @importFrom stats setNames
#' @importFrom grDevices col2rgb
#' @importFrom rlang .data
#'
make_volcano_plot_plotly <- function(df, params, title, annotation_fields = NULL, to_save = FALSE, usage = TRUE) {
  
  # Compose tooltip string based on requested annotation fields
  if (to_save == FALSE) {
    df <- df %>%
      dplyr::mutate(tooltip_text = paste0(
        "Gene: ", .data$Genes, "<br>",
        "logFC: ", round(.data$usage_logFC, 2), "<br>",
        "adjPval: ", formatC(.data$usage_adjPval, format = "e", digits = 2)
      ))
  }
  
  if (usage) { 
    df <- df %>% 
      dplyr::select(.data$usage_adjPval, .data$usage_pval, .data$usage_logFC, 
                    .data$Protein.Group, .data$Genes, .data$pep_type, 
                    dplyr::any_of("tooltip_text")) %>%  
      dplyr::mutate(adjPval = .data$usage_adjPval,
                    pval = .data$usage_pval,
                    logFC = .data$usage_logFC) 
  } else {
    df <- df %>% 
      dplyr::select(.data$adjPval, .data$pval, .data$logFC, .data$Protein.Group, 
                    .data$Genes, .data$pep_type, dplyr::any_of("tooltip_text")) 
  }
  
  # Filter significant Protein Groups
  sigPG <- df %>%
    dplyr::filter(!is.na(.data$adjPval), .data$adjPval <= params$adjPval_thr) %>%
    dplyr::arrange(.data$pval) %>% 
    dplyr::pull(.data$Protein.Group)
  
  nPOI <- function(x) sapply(seq_along(x), function(i, x) { length(unique(x[1:i])) }, x = x)
  POI_Plot <- sigPG[nPOI(sigPG) <= 10] %>% unique()
  
  POI_Plot_Genes <- df %>%
    dplyr::select(.data$Protein.Group, .data$Genes) %>%
    dplyr::filter(.data$Protein.Group %in% POI_Plot) %>%
    dplyr::pull(.data$Genes) %>% 
    unique() 
  
  adjAlpha <- params$adjPval_thr * mean(df$adjPval <= params$adjPval_thr, na.rm = TRUE)
  
  # Define Interest Levels
  levels_interest <- c(POI_Plot_Genes, "Relevant", "Not Relevant")
  
  # 1. Build the Base Color Vector
  poi_cols <- ggsci::pal_npg()(length(POI_Plot)) 
  names(poi_cols) <- POI_Plot_Genes
  my_colors <- c(poi_cols, "Relevant" = "red", "Not Relevant" = "black")
  
  # 2. Build the Alpha Vector
  my_alphas <- c(stats::setNames(rep(1, length(POI_Plot)), POI_Plot_Genes), 
                 "Relevant" = 0.2, 
                 "Not Relevant" = 0.1)
  
  # 3. Combine Colors and Alphas into Plotly-friendly RGBA strings
  rgba_colors <- mapply(function(color, alpha) {
    rgb_vals <- grDevices::col2rgb(color)[, 1]
    sprintf("rgba(%d, %d, %d, %f)", rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha)
  }, my_colors[levels_interest], my_alphas[levels_interest])
  names(rgba_colors) <- levels_interest
  
  # Mapping R pch to Plotly symbols
  pep_types <- unique(df$pep_type)
  plotly_symbols <- c("square", "circle", "triangle-up", "diamond")[1:length(pep_types)]
  names(plotly_symbols) <- pep_types
  
  # Assign relevance and interest
  df <- df %>% 
    dplyr::filter(!is.na(.data$adjPval), !is.na(.data$logFC)) %>%
    dplyr::mutate(
      relevance = ifelse(.data$adjPval <= params$adjPval_thr & abs(.data$logFC) >= params$FC_thr, 
                         "Relevant", "Not Relevant"),
      interest = ifelse(.data$Protein.Group %in% POI_Plot, 
                        as.character(.data$Genes), 
                        as.character(.data$relevance)) %>%
        factor(levels = c(POI_Plot_Genes, "Relevant", "Not Relevant")),
      
      # PRE-CALCULATE RGBA COLOR 
      point_color = mapply(function(col, alph) {
        rgb_vals <- grDevices::col2rgb(col)[, 1]
        sprintf("rgba(%d, %d, %d, %f)", rgb_vals[1], rgb_vals[2], rgb_vals[3], alph)
      }, my_colors[as.character(.data$interest)], my_alphas[as.character(.data$interest)]),
      
      # PRE-CALCULATE SYMBOL 
      point_symbol = plotly_symbols[as.character(.data$pep_type)]
    ) %>%
    dplyr::arrange(desc(.data$interest))
  
  # Subsample for performance
  df_sig <- df %>% dplyr::filter(.data$interest != "Not Relevant")
  df_noise <- df %>% dplyr::filter(.data$interest == "Not Relevant") %>% dplyr::sample_frac(0.1)
  df_for_plot <- dplyr::bind_rows(df_sig, df_noise)
  
  # ---------------------------------------------------------
  # GENERATE PLOTLY WITH SPLIT LEGENDS
  # ---------------------------------------------------------
  hx <- -100
  hy <- -100
  volcano_out <- plotly::plot_ly()
  
volcano_out <- volcano_out %>% plotly::add_trace(
    data = df_for_plot,         # 1. Provide the specific data frame here
    x = ~logFC,                 # 2. Use ~ to reference columns inside df_for_plot
    y = ~-log10(pval),          # 3. Reference 'pval' inside df_for_plot
    type = 'scatter',
    mode = 'markers',
    marker = list(
      color = ~point_color,     # 4. Use ~ for colors
      symbol = ~point_symbol,   # 5. Use ~ for symbols
      size = 7, 
      line = list(width = 0)
    ),
    text = ~tooltip_text,       # 6. Use ~ for the tooltip
    hoverinfo = "text",
    showlegend = FALSE 
  )
  
  volcano_out <- volcano_out %>% plotly::add_trace(
    x = hx, y = hy, 
    type = 'scatter', mode = 'lines',
    name = '<b>Interest</b>',
    line = list(width = 0, color = 'rgba(0,0,0,0)'),
    hoverinfo = "none",
    showlegend = TRUE
  )
  
  present_interests <- intersect(levels_interest, unique(as.character(df$interest)))
  for (int_name in present_interests) {
    volcano_out <- volcano_out %>% plotly::add_trace(
      x = hx, y = hy, 
      type = 'scatter', mode = 'markers',
      name = paste0("&nbsp;&nbsp;&nbsp;", int_name),
      marker = list(color = unname(rgba_colors[int_name]), size = 9, symbol = "circle"),
      hoverinfo = "none",
      showlegend = TRUE
    )
  }

  volcano_out <- volcano_out %>% plotly::add_trace(
    x = hx, y = hy, 
    type = 'scatter', mode = 'lines',
    name = '<br><b>Peptide type</b>',
    line = list(width = 0, color = 'rgba(0,0,0,0)'),
    hoverinfo = "none",
    showlegend = TRUE
  )
  
  for (pt_name in pep_types) {
    volcano_out <- volcano_out %>% plotly::add_trace(
      x = hx, y = hy, 
      type = 'scatter', mode = 'markers',
      name = paste0("&nbsp;&nbsp;&nbsp;", pt_name),
      marker = list(color = "grey50", size = 9, symbol = unname(plotly_symbols[pt_name])),
      hoverinfo = "none",
      showlegend = TRUE
    )
  }
  
  x_min <- min(df$logFC, na.rm = TRUE) - 0.5
  x_max <- max(df$logFC, na.rm = TRUE) + 0.5
  y_max <- max(-log10(df$pval), na.rm = TRUE) + 0.5

  volcano_out <- volcano_out %>%
    plotly::layout(
      showlegend = TRUE, 
      legend = list(
        title = list(text = ""),
        orientation = "v",      
        x = 1.02,               
        y = 1,
        xanchor = "left"
      ),
      xaxis = list(
        title = "Log<sub>2</sub>(Fold change)", 
        range = c(x_min, x_max)
      ),
      yaxis = list(
        title = "-Log<sub>10</sub>(P value)", 
        range = c(0, y_max)
      ),
      plot_bgcolor = 'white',
      paper_bgcolor = 'white',
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(adjAlpha), y1 = -log10(adjAlpha), line = list(color = "black", dash = "dash")),
        list(type = "line", x0 = -params$FC_thr, x1 = -params$FC_thr, y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1)),
        list(type = "line", x0 = params$FC_thr, x1 = params$FC_thr, y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", width = 1))
      )
    )
  
  # Environment scrubbing for smaller object size/portability
  volcano_out$x$attrs <- lapply(volcano_out$x$attrs, function(trace) {
    lapply(trace, function(attr_val) {
      if (inherits(attr_val, "formula")) {
        environment(attr_val) <- .GlobalEnv
      }
      return(attr_val)
    })
  })

  if (!is.null(volcano_out$x$data)) {
    volcano_out$x$data <- lapply(volcano_out$x$data, function(d) {
      attr(d, ".Environment") <- NULL
      return(d)
    })
  }

  volcano_out <- unserialize(serialize(volcano_out, NULL))
  
  return(list(volcano = volcano_out, POI_l = POI_Plot))
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
#' @param to_save Logical; if FALSE prepare interactive-style tooltip text (for saving as interactive plot), if TRUE produce a static plot with labeled significant points
#' @param usage Logical; if TRUE, uses "usage_" prefixed columns (default TRUE)
#' @return A list containing:
#' \itemize{
#'   \item \code{volcano}: A ggplot object representing the volcano plot.
#'   \item \code{POI_l}: A character vector of the top 10 significant Protein Groups.
#' }
#' @importFrom ggplot2 ggplot scale_alpha_manual scale_shape_manual aes theme_bw geom_point geom_vline geom_hline scale_colour_manual labs scale_x_continuous guide_legend
#' @importFrom ggsci pal_npg
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate select all_of filter distinct arrange pull desc
#' @importFrom stats setNames
#' @importFrom rlang .data

make_volcano_plot <- function(df, params, title, annotation_fields = NULL, to_save = FALSE, usage = TRUE) {
  
  # 1. Compose tooltip string
  if (!to_save && !is.null(annotation_fields)) {
    # Using a safer way to build tooltips without the '.' binding issue
    df$tooltip_text <- apply(df[, annotation_fields, drop = FALSE], 1, function(row) {
      paste(paste0(annotation_fields, ": ", row), collapse = "<br>")
    })
  }

  # 2. Handle Usage vs Abundance
  if (usage) {
    df <- df %>% 
      dplyr::select(.data$usage_adjPval, .data$usage_pval, .data$usage_logFC, 
                    .data$Protein.Group, .data$Genes, .data$pep_type) %>%
      dplyr::mutate(adjPval = .data$usage_adjPval,
                    pval = .data$usage_pval,
                    logFC = .data$usage_logFC)
  } else {
    df <- df %>% 
      dplyr::select(.data$adjPval, .data$pval, .data$logFC, 
                    .data$Protein.Group, .data$Genes, .data$pep_type)
  }

  # 3. Logic for POIs
  sigPG <- df %>%
    dplyr::filter(.data$adjPval <= params$adjPval_thr) %>%
    dplyr::arrange(.data$pval) %>%
    dplyr::pull(.data$Protein.Group)

  nPOI <- function(x) sapply(seq_along(x), function(i, x) { length(unique(x[1:i])) }, x = x)
  POI_Plot <- unique(sigPG[nPOI(sigPG) <= 10])

  POI_Plot_Genes <- df %>%
    dplyr::select(.data$Protein.Group, .data$Genes) %>%
    dplyr::filter(.data$Protein.Group %in% POI_Plot) %>%
    dplyr::pull(.data$Genes) %>%
    unique()

  adjAlpha <- params$adjPval_thr * mean(df$adjPval <= params$adjPval_thr, na.rm = TRUE)

  # 4. Color and Alpha Palettes
  poi_cols <- ggsci::pal_npg()(length(POI_Plot_Genes))
  names(poi_cols) <- POI_Plot_Genes

  my_colors <- c(poi_cols, "Relevant" = "red", "Not Relevant" = "black")
  
  my_alphas <- c(stats::setNames(rep(1, length(POI_Plot_Genes)), POI_Plot_Genes), 
                 "Relevant" = 0.2, 
                 "Not Relevant" = 0.1)

  # 5. Build Plot Data
  df <- df %>%
    dplyr::filter(!is.na(.data$adjPval)) %>%
    dplyr::mutate(
      relevance = ifelse(.data$adjPval <= params$adjPval_thr & abs(.data$logFC) >= params$FC_thr, 
                         "Relevant", "Not Relevant"),
      interest = ifelse(.data$Protein.Group %in% POI_Plot, .data$Genes, .data$relevance) %>%
                 factor(levels = c(POI_Plot_Genes, "Relevant", "Not Relevant"))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$interest))

  # 6. Generate Plot
  volcano_out <- ggplot2::ggplot(df, ggplot2::aes(x = .data$logFC, y = -log10(.data$pval))) +
    ggplot2::geom_hline(yintercept = -log10(adjAlpha)) +
    ggplot2::geom_vline(xintercept = c(-1, 1)) +
    ggplot2::geom_point(
      size = 1,
      ggplot2::aes(
        colour = .data$interest,
        alpha = .data$interest,
        shape = .data$pep_type
      )
    ) +
    ggplot2::scale_colour_manual(values = my_colors) +
    ggplot2::scale_alpha_manual(values = my_alphas) +
    ggplot2::scale_shape_manual(
      values = c(15:18),
      guide = ggplot2::guide_legend(override.aes = list(size = 3))
    ) +
    ggplot2::scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
    ggplot2::labs(
      title = paste0("Differentially ", ifelse(usage, "used", "abundant"), " LiP precursors"),
      x = expression(Log[2](Fold ~ change)),
      y = expression(-Log[10](P ~ value)),
      colour = "Interest",
      alpha = "Interest",
      shape = "Peptide type"
    ) +
    ggplot2::theme_bw()

  volcano_out$plot_env <- emptyenv()

  return(list(volcano = volcano_out, POI_l = POI_Plot))
}


#' @author Andrea Argentini
#' @title Plot Peptide Barcode for a Protein of Interest
#' @description 
#' Generates a barcode-style visualization for a specific protein. It processes 
#' differential expression results, calculates peptide coordinates (including 
#' handling of repeating sequences), determines significance based on both 
#' p-values and data completeness (missing values), and calls \code{make_barplot}.
#' 
#' @param POI Character. The protein accession (e.g., Uniprot ID) to visualize.
#' @param se A \code{SummarizedExperiment} or \code{QFeatures} object.
#' @param group_column Character. The column name in \code{colData(se)} used for grouping.
#' @param DE_result Data frame of differential expression results.
#' @param indicate_direction Logical. Currently a placeholder for directionality logic.
#' @param prot_seq Data frame containing \code{Accession} and \code{Protein.Sequence}.
#' @param usage Logical. If TRUE, uses \code{usage_} columns from DE results.
#' @param directionality Logical. If TRUE, labels significance with "Up" or "Down" strings.
#' @param expand Logical. Passed to \code{make_barplot} for faceting.
#' @param params List. Contains \code{adjPval_thr} and other thresholds.
#' 
#' @return A list containing the ggplot object (\code{gg}).
#' 
#' @importFrom dplyr filter distinct group_by mutate left_join ungroup arrange pull select join_by
#' @importFrom tidyr uncount
#' @importFrom stats setNames model.matrix
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stringr str_count str_locate_all
#' @importFrom rlang .data


plot_barcode <- function(POI, se, group_column = 'Treatment', DE_result, 
                         indicate_direction = FALSE, 
                         prot_seq,
                         usage = TRUE, directionality = FALSE, expand = FALSE, params) {

  grouping <- SummarizedExperiment::colData(se)[[group_column]]
  groups <- unique(grouping)

  # 1. Define Color Mappings
  signif_names <- c("Not Significant", paste(groups[1], "Missing"), paste(groups[2], "Missing"),
                    paste(groups[1], "Up"), paste(groups[2], "Up"), "Missing", "Significant")

  signif_colours <- c("grey40", "yellow2", "lightskyblue", "lightslateblue", "orange2", "blue", "orange")
  colour_mapping_significance <- stats::setNames(signif_colours, signif_names)

  colour_mapping_type <- c("Non-proteotypic, internally repeating" = "green", 
                           "Non-proteotypic" = "red", 
                           "Internally repeating" = "purple", 
                           "Proteotypic" = "grey80")

  # 2. Select and Rename Columns
  if (usage) {
    DE_result <- DE_result %>% 
      dplyr::select(.data$Precursor.Id, .data$usage_adjPval, .data$usage_pval, .data$usage_logFC, 
                    .data$Protein.Group, .data$Genes, .data$length, .data$pep_type, 
                    .data$Stripped.Sequence, .data$Proteotypic) %>% 
      dplyr::mutate(pval = .data$usage_pval, adjPval = .data$usage_adjPval, effectSize = .data$usage_logFC)
  } else {
    DE_result <- DE_result %>% 
      dplyr::select(.data$adjPval, .data$pval, .data$logFC, .data$Protein.Group, .data$Genes, .data$pep_type,
                    .data$Precursor.Id, .data$length, .data$Stripped.Sequence, .data$Proteotypic) %>% 
      dplyr::mutate(effectSize = .data$logFC)
  }

  # 3. Coordinate Calculation Logic
  DE_result <- DE_result %>%
    dplyr::filter(grepl(POI, .data$Protein.Group)) %>%
    dplyr::left_join(prot_seq %>% dplyr::select(.data$Accession, .data$Protein.Sequence), 
                     by = dplyr::join_by( Protein.Group == Accession)) %>% 
    dplyr::mutate(total_repeats = stringr::str_count(.data$Protein.Sequence, .data$Stripped.Sequence)) %>%
    tidyr::uncount(.data$total_repeats, .remove = FALSE) %>%
    dplyr::group_by(.data$Precursor.Id) %>%
    dplyr::mutate(
      repeat_nr = 1:dplyr::n(),
      start = stringr::str_locate_all(.data$Protein.Sequence[1], .data$Stripped.Sequence[1])[[1]][.data$repeat_nr, 1],
      end = stringr::str_locate_all(.data$Protein.Sequence[1], .data$Stripped.Sequence[1])[[1]][.data$repeat_nr, 2]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      type = dplyr::case_when(
        .data$Proteotypic == 0 & .data$total_repeats > 1 ~ "Non-proteotypic, internally repeating",
        .data$total_repeats > 1 ~ "Internally repeating",
        .data$Proteotypic == 0 ~ "Non-proteotypic",
        TRUE ~ "Proteotypic"
      ) %>% factor(levels = names(colour_mapping_type)),
      tier = match(.data$Precursor.Id, sort(unique(.data$Precursor.Id))),
      max_tier = max(.data$tier)
    )

  # 4. Completeness / Missing Values Logic
  # Use %*% for matrix multiplication and model.matrix for group counts
  comp_data <- SummarizedExperiment::assay(se)[DE_result$Precursor.Id, , drop=FALSE]
  completeness <- is.na(comp_data) %*% stats::model.matrix(~ 0 + grouping)
  
  group_counts <- as.numeric(table(grouping))
  
  DE_result <- DE_result %>% 
    dplyr::mutate(
      missing = apply(completeness == matrix(group_counts, byrow = TRUE, 
                                            nrow = nrow(completeness), 
                                            ncol = ncol(completeness)), 1, function(x) {
        h <- names(which(x)) %>% gsub(pattern = "grouping", replacement = "")
        return(if(length(h) == 0) NA else h[1])
      }),
      significance = dplyr::case_when(
        !is.na(.data$missing) ~ if(directionality) paste0(.data$missing, " Missing") else "Missing",
        .data$adjPval < params$adjPval_thr & !is.na(.data$adjPval) ~ 
           if(directionality) ifelse(.data$effectSize < 0, paste0(groups[1], " Up"), paste0(groups[2], " Up")) else "Significant",
        TRUE ~ "Not Significant"
      ) %>% factor(levels = signif_names),
      section = dplyr::case_when(
        .data$significance == "Not Significant" ~ "Not Significant",
        !is.na(.data$missing) ~ "Missing",
        TRUE ~ "Significant"
      ) %>% factor(levels = c("Not Significant", "Significant", "Missing"))
    ) %>% 
    dplyr::filter(!(.data$significance == "Not Significant" & is.na(.data$adjPval))) %>%
    dplyr::arrange(.data$significance, .data$type)

  # 5. Delegate to make_barplot
  barplot_obj <- make_barplot(df = DE_result, POI_ = POI, colour_mapping_significance, colour_mapping_type, expand)
  
  return(list(gg = barplot_obj))
}

#' @author Andrea Argentini
#' @title Create Peptide Barcode Plot
#' @description 
#' Generates a barcode-style plot using ggplot2. Draws peptide rectangles along 
#' protein sequence coordinates, colors by significance and proteotypicity, 
#' and handles faceting by peptide type and section.
#' 
#' @param df Data frame prepared for plotting (must contain start, end, tier, max_tier, length, Genes, significance, type).
#' @param POI_ Character; protein accession or identifier used in the plot title.
#' @param colour_mapping_significance Named character vector mapping significance categories to fill colors.
#' @param colour_mapping_type Named character vector mapping peptide type categories to outline colors.
#' @param expand Logical. If TRUE, facets the plot by \code{pep_type} and \code{section}.
#' 
#' @return A \code{ggplot} object.
#' 
#' @importFrom ggplot2 ggplot aes geom_rect scale_x_continuous scale_y_continuous labs scale_fill_manual scale_colour_manual theme_bw theme element_blank element_rect element_text guide_legend facet_grid
#' @importFrom rlang .data

# colour = type
make_barplot <- function(df, POI_, colour_mapping_significance, colour_mapping_type, expand) {
  
  # Note: ensure calculate_coverage is available in your package namespace
  
  plot <- ggplot(df, aes(x = .data$start, y = 1)) +
    geom_rect(
      aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = (1 / .data$max_tier) * (.data$tier - 1),
        ymax = (1 / .data$max_tier) * (.data$tier),
        fill = .data$significance,
        colour = .data$type
      ),
      linewidth = 0.2
    ) +
    scale_x_continuous(
      breaks = seq(0, 1000 * 10^(floor(log10(df$length[1] - 1))), by = 100), 
      limits = c(0, df$length[1]),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = "Significant changes by precursor",
      subtitle = paste0(
        df$Genes[1], "/", POI_, ": ",
        round(calculate_coverage(
          df$start,
          df$end,
          df$length[1]
        ) * 100, 1), 
        "% coverage"
      ),
      x = "Residue",
      y = "Precursors",
      fill = "Significance",
      colour = "Proteotypicity"
    ) +
    scale_fill_manual(values = colour_mapping_significance) +
    scale_colour_manual(
      values = colour_mapping_type, 
      guide = guide_legend(override.aes = list(fill = "transparent"))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )

  if (expand) {
    # Using explicit namespacing for facet_grid to ensure it's found
    return(plot + ggplot2::facet_grid(.data$pep_type ~ .data$section))
  } else {
    return(plot)
  }
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
        tryCatch(utils::install.packages(i , dependencies = TRUE), error = function(e) { NULL })
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