
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

#' @author Andrea Argentini
#' @title annotate_diann
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

 annotate_diann <- function(annotation_df, report_file1, report_file2 = NULL,fasta_ann){

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
#' @importFrom stringr str_split
#' @importFrom logger log_info

   parse_input <- function( input_parquet_tc, input_parquet_lip ,  dual , input_design){
     
      is_diann <- function(df) {
          "Run" %in% colnames(df)
        }
     
      tryCatch( expr = {
               
                ## code here
                if (dual == TRUE){
                    log_info('Reading Tc and Lip from SEPARATE parquet files ...')
                    TC_report <- read_parquet(input_parquet_tc)
                    LiP_report <- read_parquet(input_parquet_lip)
                  }else{
                    log_info('Reading both LiP and TC from ONE parquet file ...')
                    LiP_report <- read_parquet( input_parquet_lip)
                    TC_report <- NULL
                }
                  diann_flag <- is_diann(LiP_report)
                  log_info('Reading experiment Design file  ...')
        					design  <- read_tsv_arrow(input_design)
									## check design file. 
                  # 
                  col_design_required <-  c('Run',	'Pipeline', 	'Treatment',	'Condition'	,'Replicate', 'CondRep')
                  checkdesign <- check_design_requirement(design , col_design_required)
                  if (checkdesign$status == 1){
                    return( list(error= checkdesign$error , status= 1,lip =NULL ))
                    
                  }else{
                    return( list(error= '', status= 0,lip = LiP_report ,
                          tc = TC_report,
                          design = design,
                          diann_flag = diann_flag))
                  }
                  
                ## good exit
      },error = function(err){
                    print(paste(" Input Parquet :  ",err))
                    return( list(error= err, status= 1,lip =NULL ))
                  } )
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
#' @title consensus_normalisation
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



consensus_normalisation <- function(report){

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
#' @title dfToWideMsqrob
#' @description Convert a long precursor-level data.frame to wide format suitable for Qfeature/msqrob: keep specified feature columns and pivot the given quantification column by Run into sample columns.
#' @param data Data frame or tibble in long format containing at minimum the columns specified in wide_colums, "Run", and the quantification column named in precursorquan
#' @param precursorquan Character; name of the quantification column to spread into wide columns (e.g., "Precursor.Quantity")
#' @param wide_colums Character vector of column names to retain as identifier/feature columns prior to pivoting
#' @return A data.frame (tibble) in wide format with one row per feature (identifiers in wide_colums) and one column per Run containing the quantification values
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data

dfToWideMsqrob <- function(data, precursorquan, wide_colums) {
  data %>%
    select(
       wide_colums,
      .data[[precursorquan]]
    ) %>%
    pivot_wider(
      names_from = Run,
      values_from = .data[[precursorquan]]
    )
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

    return (list(error= '', status= 0,q_feat = pe  ))


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


make_volcano_plot <- function(df, params, title, annotation_fields , poi_vis , to_save  ) {
  # Compose tooltip string based on requested annotation fields
  if (to_save== F ){
    df <- df %>%
      mutate(
        tooltip_text = apply(
          select(., all_of(annotation_fields)),
          1,
          function(row) {
            paste(paste0(annotation_fields, ": ", row), collapse = "<br>")
          }
        )
      )
    # plotly native code 

  hline_y <- -log10(
  0.05 * sum(distinct(df, precursor.Id, .keep_all = TRUE)$adjPval <= 0.05, na.rm = TRUE) /
  nrow(distinct(df, precursor.Id, .keep_all = TRUE))
  )
 
  shape_map <- c(
  "Semi-Tryptic" = "circle",
  "Tryptic" = "square",
  "Non-Tryptic" = "triangle-up"
  )
  
  # Get interest vector
  interest_vec <- as.character(df$interest)

  # Identify POI (points of interest) levels
  poi_levels <- setdiff(unique(interest_vec), c("Significant", "Not Significant"))
  poi_levels <- poi_levels[!is.na(poi_levels)]  # remove any NA

  # Number of POI levels
  n_poi <- length(poi_levels)

  # Generate POI colors dynamically using RColorBrewer or ggplot palettes
  # Here using pal_npg from ggplot2 / ggprism / ggsci
  poi_colors <- pal_npg()(n_poi)

  # Create a named vector mapping all levels of interest to colors
  level_colors <- c(
    setNames(poi_colors, poi_levels),  # map POIs to dynamic colors
    "Significant" = "red",
    "Not Significant" = "black"
  )
  
  df$plotly_symbol <- shape_map[ as.character(df$pep_type) ]  
  plot_ly(
    data = df,
    x = ~logFC,
    y = ~-log10(pval),
    text = ~tooltip_text,
    color = ~interest,
    colors = level_colors,
    type = "scatter",
    mode = "markers",
    marker = list(
       symbol = df$plotly_symbol, 
      size = 6,
      line = list(width = 0.5, color = "grey")
    ),
    hoverinfo = "text"
  ) %>%
    add_lines(
      x = c(- params$FC_thr, - params$FC_thr),
      y = c(0, max(-log10(df$pval), na.rm = TRUE)),
      line = list(color = "grey", dash = "dot"),
      name = paste0("Log2FC = ", params$FC_thr),
      inherit = FALSE,
      showlegend = FALSE
    )   %>% add_lines(
      x = c( params$FC_thr,  params$FC_thr),
      y = c(0, max(-log10(df$pval), na.rm = TRUE)),
      line = list(color = "grey", dash = "dot"),
      name = paste0("Log2FC = ", params$FC_thr),
      inherit = FALSE,
      showlegend = FALSE
    ) %>% 
    add_lines(
      x = range(df$logFC, na.rm = TRUE),
      y = c(hline_y, hline_y),
      line = list(color = "grey", dash = "dot"),
      name = "Adjusted p-value threshold",
      inherit = FALSE,
      showlegend = FALSE
    ) %>% layout(
       title = list(text = "Differential LiP precursors",  x = 0.5),
      xaxis = list(title = "Log2(Fold change)", zeroline = FALSE),
      yaxis = list(title = "-Log10(P value)", zeroline = FALSE),
      legend = list(title = list(text = "Legend"), orientation = "v")
    )
  }else{

    ggplot(df, aes(x = logFC, y = -log10(pval))) +
      geom_point(size = 1,
             
             aes(colour= interest,  shape = pep_type , alpha = interest),
             show.legend = T) +
    theme_bw() +
    geom_vline(xintercept = c(-params$FC_thr, params$FC_thr), col = "grey") +
    geom_hline(yintercept = -log10(0.05*sum(distinct(df, precursor.Id, .keep_all = T)$adjPval<=0.05, na.rm = T)/nrow(distinct(df, precursor.Id, .keep_all = T))),
                 col = "grey")  +
    scale_colour_manual(values = c(pal_npg()(length(poi_vis)), ifelse("Significant" %in% df$interest, "red", "black"), "black")) +
    scale_shape_manual(values = c(20, 18, 15) , guide = guide_legend(override.aes = list(size = 2))) +
    scale_alpha_manual(values = c(rep(1, length(poi_vis)), 0.2, 0.1), 
                     guide = guide_legend(override.aes = list(alpha = 1))) +
    scale_x_continuous(breaks = -100:100*2) +
    geom_text_repel( data =  df %>% filter(interest != "Not Significant"),
                  aes(label = Genes,
                      colour = interest),
                  size = 2,
                  segment.size = 0.25,
                  show.legend = F) +
    labs(title = title,
       x = expression(Log[2](Fold~change)),
       y = expression(-Log[10](P~value)),
        colour = "Protein",
       alpha = "Protein",
       shape = "Peptide type")
  }
 
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
#browser()
barcode_plotly_list <- lapply(barcode, function(el) {
   #el$gg is the ggplot, el$colour_mapping_significance is the mapping
  plotly_from_ggplot(el$gg, el$colour_mapping_significance, title_center = 0.0, top_margin = 80, tooltip = "text")
})
  
return ( list( toptable =DEall , volcano = volcano, volcano2file = p_toFile , barcode_gg = barcode , barcode_plty = barcode_plotly_list,  POI = POIs ) )

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
#' @importFrom stats setNames



plot_barcode <- function(POI, input, group_column, DE_result,  indicate_direction = F){
groups <- input %>% pull(group_column) %>% unique()
#browser()
signif_names <- c(
  "Not Significant",
  paste(groups[1], "Up"),
  paste(groups[2], "Up"),
  paste(groups[1], "Missing"),
  paste(groups[2], "Missing"),
  "Significant",
  "Missing"
)
signif_colours <- c(
  "grey40",
  "lightslateblue",
  "orange2",
  "yellow2",
  "lightskyblue",
  "orange",
  "blue"
)
colour_mapping_significance <- setNames(signif_colours, signif_names)
colour_mapping_type <- c("Non-proteotypic, internally repeating" = "green",
                    "Non-proteotypic" = "red",
                    "Internally repeating" = "purple",
                    "Proteotypic" = "grey80")

part_vis <-
  input %>%
  filter(Accession %in% POI) %>%
  group_by(Precursor.Id, .data[[group_column]] ) %>% #change to group_column
  mutate(median_abundance = median(normPQ, na.rm = T)) %>%
  distinct(Precursor.Id, .data[[group_column]], repeat_nr, .keep_all = T) %>%
  group_by(Precursor.Id) %>%
  mutate(coverage = round(coverage*100, 2),
         directionality = ifelse(length(unique(.data[[group_column]]))==2,
                                 ifelse(median_abundance[.data[[group_column]] ==groups[2]]>median_abundance[ .data[[group_column]] ==groups[1]],
                                        groups[2],
                                        groups[1]),
                                 ifelse(unique( .data[[group_column]])==groups[2], #check if correct
                                        groups[1],
                                        groups[2])),
         completeness = ifelse(length(unique(.data[[group_column]]))==2,
                               "Up",
                               "Missing")) %>%
  distinct(Precursor.Id, repeat_nr, .keep_all = T) %>%
  left_join(DE_result %>% select(precursor.Id, logFC, pval, adjPval),
            join_by(Precursor.Id== precursor.Id)) %>%
  mutate(significance =
           ifelse(indicate_direction,
             ifelse(any(adjPval <= 0.05, all(is.na(adjPval), completeness=="Missing"), na.rm = T),
                    paste(directionality, completeness),
                    "Not Significant"),
             ifelse(is.na(adjPval), # check if completeness=="Missing"
                    "Missing",
                    ifelse(adjPval<=0.05,
                           "Significant",
                           "Not Significant"))) %>%
           factor(levels =  c("Not Significant",
                              paste(rep(groups, 2), rep(c("Missing", "Up"), each = 2)),
                              "Missing", "Significant")),
           type = ifelse(Proteotypic==0 & total_repeats>1, "Non-proteotypic, internally repeating",
                       ifelse(total_repeats>1, "Internally repeating",
                              ifelse(Proteotypic==0, "Non-proteotypic",
                                     "Proteotypic"))) %>%
           factor(levels = c("Non-proteotypic, internally repeating",
                             "Non-proteotypic",
                             "Internally repeating",
                             "Proteotypic")) ,
                  section = ifelse(significance == "Not Significant",
                          "Not Significant",
                          ifelse(completeness=="Missing",
                                 "Missing",
                                 "Significant")) %>% 
                  factor(levels = c("Not Significant", "Significant", "Missing")),
                 peptype = case_match(pep_type,
                              "SemiTryptic" ~ "Semi-Tryptic",
                              "Tryptic" ~ "Tryptic",
                              "NonTryptic" ~ "Non-Tryptic")
                            ) %>%
  group_by(Stripped.Sequence) %>%
  mutate(tier = match(Precursor.Id, sort(unique(Precursor.Id))),
         max_tier = max(tier)) %>%
  ungroup() %>%
  arrange(significance, type)

  barplot_obj <- make_barplot(part_vis, POI_ = POI, colour_mapping_significance, colour_mapping_type)
 
  return (list(
        gg = barplot_obj,
        colour_mapping_significance = colour_mapping_significance,
        colour_mapping_type = colour_mapping_type
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
make_barplot <- function ( df, POI_ , colour_mapping_significance, colour_mapping_type){
 
   ggplot(df , aes(x=start, y = 1)) +
  geom_rect(aes(xmin = start,
                xmax = end,
                ymin = (1/max_tier)*(tier-1),
                ymax = (1/max_tier)*(tier),
                fill = significance,
                group = significance,
               colour = type, 
              text = paste0(
              "Peptide type: ", peptype, "<br>",
              "Significance: ", significance, "<br>"
               )),
              ,
            linewidth = 0.2) +
  scale_x_continuous(breaks = c(0:1000*10^(floor(log10(df$length[1]-1)))), 
                       limits = c(0, df$length[1]),
                       expand = c(0,0)
                       ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = paste("Significant changes by precursor: "),
       subtitle = paste(df$Genes[grep(POI_, df$Accession)[1]], "/", POI_, ": ", df$coverage[1], "% coverage", sep = ""),
       x = "Residue",
       y = "Precursors",
       fill = "Significance",
       colour = "Proteotypicity") +
  scale_fill_manual(values = colour_mapping_significance) +
  scale_colour_manual(values = colour_mapping_type, guide = "none") + # hide colour legend 
  #scale_colour_manual(values = colour_mapping_type,
  #                    guide = guide_legend(override.aes = list(fill = "transparent"))) +
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
}

### to be added 
plotly_from_ggplot <- function(ggp, colour_mapping_significance, title_center = 0.0, top_margin = 70, tooltip = "text") {
  # Requires: library(plotly)
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