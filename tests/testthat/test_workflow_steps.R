library(testthat)

library(arrow)

test_that("dialipar_checkdesignI", {

   params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation.txt')
   #params$fasta_ <-  'SP_9606_PK.fasta' 

  out <-  parse_input( params$input_tc, params$input_lip ,  dual = TRUE , params$design)
  # Perform the tests
  #print(import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = lst_wide_columns  )$error)
  expect_equal(out$status,0) # Ensure the data is loaded  

} 
)

test_that("dialipar_checkdesignII", {

   params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation_wrong.txt')
   #params$fasta_ <-  'SP_9606_PK.fasta' 

  out <-  parse_input( params$input_tc, params$input_lip ,  dual = TRUE , params$design)
  # Perform the tests
  #print(import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = lst_wide_columns  )$error)
  expect_equal(out$status,1) # Ensure the data is loaded  

} 
)


test_that("dialipar_fastaCheck", {
			params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation.txt')
   params$fasta_file <-  test_path('9606_wrong.fasta') 

  	suppressWarnings({
  	fastaproc <- read_fasta_ann(params$fasta_file )
     } )
   ##print(fastaproc$result %>%  distinct(Accession)  )
   expect_true(all(is.na(fastaproc$result %>%  distinct(Accession))))
    
      #})
  	
} )


test_that("dialipar_workflowI", {
			params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation.txt')
   params$fasta_file <-  test_path('SP_9606_PK.fasta') 

			inputproc <-  parse_input( params$input_tc, params$input_lip ,  dual = TRUE , params$design)
  	suppressWarnings({
  		fastaproc <- read_fasta_ann(params$fasta_file )})
  	annproc <- annotate_spectronaut( inputproc$design, inputproc$lip, inputproc$tc, fastaproc$result)
			res_norm <- consensus_normalisation(annproc$result)
			LiP_annotated <- res_norm$normalized %>%  filter(Pipeline=="LiP")
  	tc_annotated <- res_norm$normalized %>%  filter(Pipeline=="TC")
   saveRDS(tc_annotated,test_path('tc-processed.Rds') )
   expect_equal(annproc$status,0) # spectronaut data is annotated
   expect_equal(res_norm$status,0) # normalization is ok 
			expect_equal(dim(res_norm$normalized)[1],12) # dim normalized data
   expect_equal(dim(res_norm$normalized)[2],30) # dim normalized data
  	expect_equal(dim(res_norm$normalized %>%  filter(Pipeline=="LiP"))[1],8) # dim norm data LiP
   expect_equal(dim(res_norm$normalized %>%  filter(Pipeline=="TC"))[1],4) # dim norm data TC
} )



test_that("dialipar_qfI", {
			params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation.txt')
   params$fasta_file <-  test_path('SP_9606_PK.fasta') 
   tc_input <-readRDS(test_path('TC_processed_full.Rds'))
  
			inputproc <-  parse_input( params$input_tc, params$input_lip ,  dual = TRUE , params$design)
			coverageproc_tc <- calculate_coverages(tc_input)
  
			col  <- c('Run', 'Precursor.Id', 'total_repeats','repeat_nr','Modified.Sequence', 
          'Stripped.Sequence', 'Accession', 'Protein.Group','Protein.Names','Genes','pep_type','Proteotypic')
   tcproc <- input_qf(coverageproc_tc$result  , inputproc$design , 
																							columns_not_wide = col, 
																							flag_tc = TRUE )
    expect_equal (dim(tcproc$qf_pe[['precursor']])[1], 182 )
				expect_equal( dim(colData(tcproc$qf_pe)) [1], 8 )
  		expect_equal( dim(colData(tcproc$qf_pe)) [2], 6)

			
} )