library(testthat)

test_that("parsing_base", {

    params_start <- list()
    params_start$description <-  "DIA-LiPA from DIA-NN input (TC 1 LiP in the same parquet file)"
    params_start$design_file <-  test_path("SampleAnnotation.txt")
    params_start$input_file_lip <-  test_path("LiP_small.parquet")
    params_start$input_file_tc <-  ''
    params_start$fasta_file <-  test_path('SP_9606_PK.fasta' )
    params_start$folder_prj <- 'Check'
    params_start$title <-  "Dev Report "
    params_start$subtitle <-  "DIA-LiPA"
    params_start$author <-  "The GateKeeper" 
    params_start$formula <-  '~  -1 + Condition '
    params_start$comparisons <- c('ConditionRapa_LiP - ConditionDMSO_LiP')
    params_start$FC_thr <- 1
    params_start$adjPval_thr < 0.1
    params_start$comparison_label <- c('Rapa - Dmso')
    params_start$paired <- TRUE
    params_start$poi <- c('P62942','Q02790','Q00688')
   res <- validate_params_minimal(params_start )
   #print(params_start)
  # Perform the tests
  expect_true(res)

} )


test_that("validate_params_minimal_complexformula", {
  params_start <- list()

    params_start$description <-  "DIA-LiPA from DIA-NN input (TC 1 LiP in the same parquet file)"
    params_start$design_file <-  test_path("SampleAnnotation.txt")
    params_start$input_file_lip <-  test_path("LiP_small.parquet")
    params_start$input_file_tc <-  ''
    params_start$fasta_file <-  test_path('SP_9606_PK.fasta' )
    params_start$folder_prj <- 'Check'
    params_start$title <-  "Dev Report "
    params_start$subtitle <-  "DIA-LiPA"
    params_start$author <-  "The GateKeeper" 
    params_start$formula <-  '~  -1 + Condition:Time '
    params_start$comparisons <- c('ConditionRapa_LiP - ConditionDMSO_LiP')
    params_start$FC_thr <- 1
    params_start$adjPval_thr < 0.1
    params_start$comparison_label <- c('Rapa - Dmso')
    params_start$paired <- TRUE
    params_start$poi <- c('P62942','Q02790','Q00688')

  
  expect_error(
    validate_params_minimal(params_start),
    regexp = "Formula contains complex terms"   # checks the error message
  )
})




test_that("validate_params_minimal_missingcomparison", {
  params_start <- list()
  params_start$formula <- '~ Condition'  # something invalid
    params_start$comparisons <- ''
    params_start$FC_thr <- 1
    params_start$adjPval_thr < 0.1
    params_start$comparison_label <- c('Rapa - Dmso')
    params_start$poi <- c('P62942','Q02790','Q00688')
    params_start$design_file <-  test_path("SampleAnnotation.txt")
    params_start$input_file_lip <-  test_path("LiP_small.parquet")
     params_start$input_file_tc <- ''
    params_start$fasta_file <-  test_path('SP_9606_PK.fasta' )
    params_start$folder_prj <- 'Check '
    params_start$description <-  "DIA-LiPA from DIA-NN input (TC 1 LiP in the same parquet file)"
    params_start$title <-  "Dev Report "
    params_start$subtitle <-  "DIA-LiPA"
    params_start$author <-  "The GateKeeper" 
    params_start$paired <- TRUE

  
  expect_error(
    validate_params_minimal(params_start),
    regexp = "Comparisons must contain at least one value."   # checks the error message
  )
})


test_that("validate_params_minimal_missingLiPfile", {
  params_start <- list()
  params_start$formula <- '~ Condition' 
    params_start$comparisons <- ''
    params_start$FC_thr <- 1
    params_start$adjPval_thr < 0.1
    params_start$comparison_label <- c('Rapa - Dmso')
    params_start$poi <- c('P62942','Q02790','Q00688')
    params_start$design_file <-  test_path("SampleAnnotation.txt")
    params_start$input_file_lip <- ''
    params_start$input_file_tc <- ''
    params_start$fasta_file <-  test_path('SP_9606_PK.fasta' )
    params_start$folder_prj <- 'Check'
    params_start$description <-  "DIA-LiPA from DIA-NN input (TC 1 LiP in the same parquet file)"
    params_start$title <-  "Dev Report "
    params_start$subtitle <-  "DIA-LiPA"
    params_start$author <-  "The GateKeeper" 

  
  expect_error(
    validate_params_minimal(params_start),
    regexp = "Input file LiP does not exist or is not specified."   # checks the error message
  )
})