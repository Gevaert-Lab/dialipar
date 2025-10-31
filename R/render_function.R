

#' @author andrea Argentini
#' @title Validate template parameter
#' @description Validate that the provided template is a single string and belongs to the allowed set of templates. Stops with an informative error if validation fails; returns TRUE when valid.
#' @param template Character scalar specifying the template filename to validate
#' @return Logical TRUE if the template is valid; otherwise the function stops with an error message
#' @importFrom assertthat assert_that is.string

validate_template <- function(template) {
  # Define the list of valid templates
  valid_templates <- c( "Template_.qmd")

  # Check if template is a string
  assertthat::assert_that(assertthat::is.string(template), msg = "template must be a string.")

  # Check if template belongs to the list of valid templates
  if (!template %in% valid_templates) {
    stop("Invalid template. The template must be one of the following: ", paste(valid_templates, collapse = ", "))
  }

  TRUE
}


#' @author andrea argentini
#' @title Validate file name
#' @description Validate a proposed file name. Ensures the input is a single string, is shorter than 40 characters, and does not contain invalid filesystem characters (<>:"/\\|?*). The function stops with an informative error if validation fails and returns TRUE when valid.
#' @param filename Character scalar; the file name to validate
#' @return Logical TRUE if the file name is valid; otherwise the function stops with an error message
#' @importFrom assertthat assert_that is.string

validate_filename <- function(filename) {
  # Define invalid characters for file names
  invalid_chars <- "[<>:\"/\\|?*]"

  # Check if filename is a string
  assertthat::assert_that(assertthat::is.string(filename), msg = "filename must be a string.")

  # Check if filename length is less than 40 characters
  if (nchar(filename) > 40) {
    stop("The file name must be less than 40 characters.")
  }

  # Check if filename contains invalid characters
  if (grepl(invalid_chars, filename)) {
    stop("The file name contains invalid characters. Invalid characters are: <>:\"/\\|?*")
  }

  TRUE
}


#' @author andrea Argentini
#' @title Validate report folder path
#' @description Validate and prepare the report folder. Creates the folder (and a "Result" subfolder) if it does not exist. Returns TRUE when the folder is available; otherwise the function will stop with an error.
#' @param report_folder Character scalar. Path to the folder to validate/create.
#' @return Logical TRUE if the folder is valid/created (invisibly); otherwise the function stops with an error.
#' @importFrom assertthat assert_that is.writeable
validate_folder <- function(report_folder) {
  # Define invalid characters for Windows file system
  invalid_chars <- "[<>:\"/\\|?*]"

  # # Check if the folder path contains invalid characters
  # if (grepl(invalid_chars, report_folder)) {
  #   stop("The folder path contains invalid characters. Invalid characters are: <>:\"/\\|?*")
  # }

  if (!dir.exists(file.path(report_folder))) {
    dir.create(file.path( report_folder),recursive = TRUE)
  }
  dir.create(file.path( report_folder, "Result"),recursive = TRUE)

  # Check if the folder path is writable
  #assertthat::assert_that(assertthat::is.writeable(report_folder), msg = "The folder path is not writable.")

  TRUE
}


#' @author andrea Argentini
#' @title Validate minimal parameters for the report
#' @description Validate a minimal set of parameters required to run the report pipeline. Checks types and existence of key files/values (e.g., input_file_tc, input_file_lip, fasta_file, design_file), validates formula syntax (uses flag_complex_formula), and ensures numeric/character parameters meet basic constraints. Stops with informative messages on failure; returns TRUE on success.
#' @param params List. A named list of parameters to validate. Expected keys include (but are not limited to):
#'   \describe{
#'     \item{input_file_tc}{character; path to TC input parquet or file (may be empty string to skip).}
#'     \item{input_file_lip}{character; path to LiP input file.}
#'     \item{fasta_file}{character; path to FASTA file.}
#'     \item{design_file}{character; path to design file.}
#'     \item{folder_prj}{character; project/folder path.}
#'     \item{description, title, subtitle, author}{character; metadata strings.}
#'     \item{formula}{character; model formula string (will be checked for disallowed complex terms).}
#'     \item{comparisons}{character; vector of comparison labels (must contain at least one).}
#'     \item{FC_thr}{numeric; fold-change threshold (must be positive).}
#'     \item{comparison_label}{character; labels for comparisons.}
#'     \item{poi}{character; optional proteins/IDs of interest.}
#'   }
#' @return Logical TRUE if all required checks pass; otherwise the function stops with an informative error message.
#' @importFrom assertthat assert_that is.string
validate_params_minimal <- function(params) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  is_empty <- function(x) {
    is.character(x) && length(x) == 1 && x == ''
  }

  check_path <- function(x) {
    if (length(x) == 1 && x == '') {
      TRUE
    } else {
      file.exists(x)
    }
  }

   check_select_group <- function(x) {

    # Must be a named list
    if (!is.list(x) || is.null(names(x))) {
      return(FALSE)
    }else{
      return (TRUE)
    }
   }
  requirements <- list(
    input_file_tc = list(
      type = "string",
      check = function(x) check_path(x) ,
      msg = "Input file TC does not exist."
    ),
     input_file_lip = list(
      type = "string",
      check = function(x) file.exists(x),
      msg = "Input file LiP does not exist or is not specified."
    ),
     fasta_file = list(
      type = "string",
      check = function(x) file.exists(x),
      msg = "Input file does not exist or is not specified."
    ),
    design_file = list(
      type = "string",
      check = function(x) file.exists(x),
      msg = "Design file does not exist or is not specified."
    ),
    folder_prj = list(type= "string"
    )
    ,
    description= list(
      type = "string"
    ),
    title= list(
      type = "string"
    ),
    subtitle = list(
      type = "string"
    ),
    author = list(
      type = "string"
    ),
   formula = list(
    type = "string",
      check = function(x) {
      res <- flag_complex_formula(as.formula(x))
      if (res$flag) {
        stop(paste(
          "Formula contains complex terms:",
          paste(res$problematic_terms, collapse = ", ")
        ))
      }
      TRUE
        }
  ),
     comparisons = list(
      type = "character",
      check = function(x) length(x) >= 1 && all(x != ''),
      msg = "Comparisons must contain at least one value."
    ),
     FC_thr = list(
      type = "numeric",
      check = function(x) x > 0,
      msg = "FC_thr must be a positive number."
    ),
     comparison_label = list(
      type = "character",
      check = function(x) length(x) >= 1 && all(x != ''),
      msg = "Comparison label must contain at least one value."
    ),
     poi = list(
      type = "character"
    )


  )

  for (p in names(requirements)) {
    val <- params[[p]]
    req <- requirements[[p]]

    # Type check
    if (req$type == "string") {
      assertthat::assert_that(assertthat::is.string(val), msg = paste0("'", p, "' must be a string."))
    } else if (req$type == "numeric") {
      assertthat::assert_that(is.numeric(val), msg = paste0("'", p, "' must be numeric."))
    } else if (req$type == "logical") {
      assertthat::assert_that(is.logical(val), msg = paste0("'", p, "' must be logical (TRUE/FALSE)."))
    } else if (req$type == "character") {
      assertthat::assert_that(is.character(val), msg = paste0("'", p, "' must be a character vector."))
    } else if (req$type == "list") {
      assertthat::assert_that(is.list(val), msg = paste0("'", p, "' must be a list."))
    }

    # Value check (if provided)
    if (!is.null(req$check)) {
      assertthat::assert_that(req$check(val), msg = req$msg %||% paste0("Invalid value for '", p, "'."))
    }
  }
  TRUE
}


#' @author Andrea Argentini
#' @title merge_default_parameters
#' @description Merge a user-supplied parameter list with package defaults. The function reads defaults from config/default_parameter.yaml in the package and fills any parameters missing from params_int with those default values, returning the merged parameter list.
#' @param params_int Named list of input parameters provided by the user (overrides defaults)
#' @return Named list: the merged parameters where missing entries from params_int are filled using the package defaults
#' @importFrom yaml read_yaml

merge_default_parameters <- function  ( params_int  ){

  yaml_path <- system.file("config", "default_parameter.yaml", package = "dialipar")
  default_p <- read_yaml(yaml_path )

  miss <- base::setdiff(names(default_p$params), names(params_int))
   for (a in miss) {
     params_int[[a]] <- default_p$params[[a]] }

  return (params_int)
}

#' @author andrea argentini
#' @title Render a DIA-LiPA report using a Quarto template
#' @description Render a full DIA-LiPA HTML report from processed input files and a Quarto template.
#' The function validates parameters and paths, prepares a temporary rendering workspace (copies the package template and RDS results),
#' executes a Quarto render with the provided parameters, and copies the generated HTML (and resource folder) to the requested output folder.
#' It also initializes logging to both console and a file inside the report folder.
#' @param params_report Named list of analysis and report parameters. See Details for required and commonly used entries.
#' @param template Character; the filename of the Quarto template inside the package template folder (e.g. "Template_.qmd")
#' @param report_folder Character; directory where the final rendered report and resources should be written (will be created if missing)
#' @param report_filename Character; output HTML filename for the rendered report (e.g. "my_dialip_report.html")
#' @details
#' Required and commonly used elements inside params_report:
#' \describe{
#'   \item{input_file_tc}{Path to TC parquet file (can be an empty string '' if not provided and LiP-only processing is desired).}
#'   \item{input_file_lip}{Path to LiP parquet file (required).}
#'   \item{fasta_file}{Path to FASTA file used for sequence annotation (required).}
#'   \item{design_file}{Path to the experiment design file (tsv/arrow readable) (required).}
#'   \item{formula}{Model formula as a string (e.g. "~ Condition + Batch"). Complex formula terms (interactions, functions) are rejected.}
#'   \item{comparisons}{Character vector of comparisons to test (e.g. c("ConditionB-ConditionA")).}
#'   \item{comparison_label}{Character vector of labels for the comparisons (used for naming outputs).}
#'   \item{FC_thr}{Numeric; fold-change threshold used in plotting/reporting (e.g. 1).}
#'   \item{adjpval_thr}{Numeric; adjusted p-value threshold for significance (e.g. 0.05).}
#'   \item{poi}{Optional character vector of protein accessions (Uniprot IDs) to highlight in plots.}
#' }
#' The function will merge params_report with package defaults (from config/default_parameter.yaml) and then validate the minimal required fields.
#' It creates a temporary working directory, copies the package Quarto template and supporting files, saves intermediate RDS objects used by the template,
#' runs quarto::quarto_render with execute_params set to the final parameter list, and then copies the resulting HTML and resource folder to report_folder.
#' @return Character scalar: full path to the rendered HTML report file (in report_folder). The function will stop with an error if rendering fails.
#' @examples
#' \dontrun{
#' params <- list(
#'   input_file_tc = '',
#'   input_file_lip = "data/my_lip.parquet",
#'   fasta_file = "data/uniprot.fasta",
#'   design_file = "data/design.tsv",
#'   formula = "~ Condition",
#'   comparisons = c("ConditionB-ConditionA"),
#'   comparison_label = c("B_vs_A"),
#'   FC_thr = 1,
#'   adjpval_thr = 0.05,
#'   poi = character(0)
#' )
#'
#' # Renders the report (writes my_report.html and associated resources into report_folder)
#' render_dialipa_report(params_report = params,
#'                     template = "Template_.qmd",
#'                     report_folder = file.path(tempdir(), "dialip_report"),
#'                     report_filename = "my_report.html")
#' }
#' @export
#' @importFrom quarto quarto_render
#' @importFrom fs file_move
#' @importFrom logger log_info log_threshold log_appender log_formatter INFO appender_console appender_file
#' @importFrom yaml as.yaml
#' @importFrom utils modifyList
#' @importFrom withr with_dir
#' @importFrom assertthat assert_that is.string

render_dialipa_report <- function(params_report, template, report_folder, report_filename ) {

  # Validate parameters
  validate_template( template)
  validate_folder(report_folder)
  validate_filename( filename = report_filename)


  params_report <- merge_default_parameters(params_report)

  validate_params_minimal(params_report)


  logger::log_threshold(logger::INFO)
  logger::log_appender(logger::appender_console)
  logger::log_formatter(logger::formatter_glue)
  logfile <- file.path(report_folder, "logfile_dialipa.log")
  file.create(logfile)
  logger::log_appender(logger::appender_file(logfile ), index = 2)

  log_info ('DIA-LiPA start  ...')

  ## Drafting the  flow
  if (params_report$input_file_tc == ''){
      inputproc  <- parse_input ( params_report$input_file_tc, params_report$input_file_lip ,  dual = FALSE, params_report$design_file)

  }else{
          inputproc  <- parse_input ( params_report$input_file_tc, params_report$input_file_lip ,  dual = TRUE, params_report$design_file)
  }
  if (inputproc$status == 1) stop(inputproc$error)
  fastaproc <- read_fasta_ann(params_report$fasta_file )
  
  if (fastaproc$status == 1) stop(fastaproc$error)

  if (inputproc$diann_flag == TRUE) {
      annproc <- annotate_diann( inputproc$design, inputproc$lip, inputproc$tc, fastaproc$result)

  }else{
      annproc <- annotate_spectronaut( inputproc$design, inputproc$lip, inputproc$tc, fastaproc$result)

  }
  log_info( paste('Dim annotate_diann ', dim(annproc$result), collapse = ' '))
  if (annproc$status == 1) stop(annproc$error)
  consproc <- consensus_normalisation(annproc$result)
  
  if (consproc$status == 1) stop(consproc$error)

  LiP_annotated <- consproc$normalized %>%  filter(Pipeline=="LiP")
  TC_annotated <-  consproc$normalized %>% filter( Pipeline=="TC")

  log_info( paste('Dim LiP_annotated ', dim(LiP_annotated), collapse = ' '))
  
  log_info( paste('Dim TC_annotated ', dim(TC_annotated), collapse = ' '))

  coverageproc_lip  <- calculate_coverages(LiP_annotated)
  coverageproc_tc <- calculate_coverages(TC_annotated)
  complete_report <- bind_rows(coverageproc_lip$result, coverageproc_tc$result)
  log_info( paste('Dim complete_report ', dim(complete_report), collapse = ' '))

  diann_col <- c('Run', 
      'Precursor.Id', 
      'pep_type',
      'total_repeats',
      'repeat_nr',
      'start',
      'end',
      'Modified.Sequence', 
      'Stripped.Sequence', 
      'Accession',
      'Protein.Group',
      'Protein.Names',
      'Genes',
      'pep_type',
      'Proteotypic')
  
  
  tcproc <- input_qf (coverageproc_tc$result , inputproc$design , columns_not_wide = diann_col, flag_tc = TRUE )

  tc_proc_scaling <- processing_tc_qfeat(tcproc$qf_pe, inputproc$design )
  ## TODO hard oded 
  LiP_annotated_corr <- 
  coverageproc_lip$result %>% 
      mutate( ID = paste(Protein.Group, Treatment),
            abundance_adjustment = tc_proc_scaling$adj_scaling_df$abundance_adjustment[ match(ID,   tc_proc_scaling$adj_scaling_df$ID)] %>% 
              ifelse(is.na(.), 0, .),
            adjPQ = normPQ-abundance_adjustment ) 

  res_lip <- input_qf ( LiP_annotated_corr , 
                     inputproc$design , 
                      columns_not_wide = diann_col, 
                      flag_tc = FALSE )
  res_de <-  msqrob_model(pe = res_lip$qf_p, params = params_report, layer = 'precursor' )

  res_DE_ <-  lapply(params_report$comparisons, dep_volcano_barcode,
                    data= res_de$q_feat  ,
                    params = params_report ,
                    df_anno =LiP_annotated_corr ,
                    layer= 'precursor' )
  names(res_DE_) <-  params_report$comparison_label

  ## -- template creation legacy code.

  template_source_folder <- system.file("quarto_template", package = "dialipar")
  if (template_source_folder == "") {
    stop("Template folder not found in the package.")
  }


  # Create a unique temporary working directory
  temp_work_dir <- file.path(tempdir(), paste0("quarto_temp_", Sys.getpid()))
  dir.create(temp_work_dir, recursive = TRUE, showWarnings = FALSE)
  log_info('Temp folder created : {temp_work_dir}')

  # Copy the entire template folder content to the temporary directory
  # This copies all files and subfolders (e.g., resource folders with JS/CSS files)
  success <- file.copy(from = template_source_folder,
                       to = temp_work_dir,
                       recursive = TRUE)
  if (!success) {
    stop("Failed to copy the template folder to the temporary directory.")
  }
  log_info('Copy template file ...done')
  saveRDS(  coverageproc_lip$result , file.path(temp_work_dir,basename(template_source_folder), 'lip_comp.RDS'  ))
  params_report$lip_rep <-   file.path(temp_work_dir,basename(template_source_folder),'lip_comp.RDS'  )

  saveRDS(complete_report, file.path(temp_work_dir,basename(template_source_folder), 'complete_.RDS'  ))
  params_report$complete_rep <-   file.path(temp_work_dir,basename(template_source_folder),'complete_.RDS'  )

  saveRDS(res_DE_, file.path(temp_work_dir,basename(template_source_folder), 'resDE_.RDS'  ))
  params_report$res_DE <-   file.path(temp_work_dir,basename(template_source_folder),'resDE_.RDS' )
  log_info('Copy Rds results ...done')
   # Construct the path to the copied template file in the temp directory.
  # Assumes that the template file is directly inside the copied folder.
  temp_template_path <- file.path(temp_work_dir, basename(template_source_folder), template)
  if (!file.exists(temp_template_path)) {
    stop("Template file not found in the temporary directory: ", temp_template_path)
  }

  path <- file.path(temp_work_dir, basename(template_source_folder))

  tryCatch({
    with_dir(path, {
      quarto_render(
        input = temp_template_path,
        output_format = "html",
        output_file = report_filename,
        execute_params = params_report,
        quarto_args = c("--output-dir", path)
      )
    })
  }, error = function(e) {
    print("Error in Quarto rendering:")
    print(e$message)
    print("Cleaning Temp folder")
    unlink(temp_work_dir, recursive = TRUE)
    stop(e)
  })

  resource_folder_name <- paste0(tools::file_path_sans_ext(report_filename), "_files")
  rendered_report_path <- file.path(path, report_filename)

  if (!dir.exists(report_folder)) {
    dir.create(report_folder, recursive = TRUE)
  }

  log_info('Copying rendered html report ...')
  # Copy the rendered HTML report to the target folder
  file.copy(from = rendered_report_path, to = file.path(report_folder, report_filename), overwrite = TRUE)
  # If a resource folder was generated, copy it as well
  temp_resource_path <- file.path(temp_work_dir, resource_folder_name)
  if (dir.exists(temp_resource_path)) {
    file.copy(from = temp_resource_path,
              to = file.path(report_folder, resource_folder_name),
              recursive = TRUE, overwrite = TRUE)
  }

  log_info('Cleaning temp folder ...')
  # Optionally, remove the temporary working directory to clean up
  unlink(temp_work_dir, recursive = TRUE)
  return (-1)
}
