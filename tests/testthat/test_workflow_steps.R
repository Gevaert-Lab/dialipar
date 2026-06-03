library(testthat)

library(arrow)

test_that("dialipar_checkdesignI", {

   params <- list()
   params$input_tc <- test_path('TC_small.parquet')
   params$input_lip <- test_path('LiP_small.parquet')
   params$design <- test_path('SampleAnnotation.txt')
   #params$fasta_ <-  'SP_9606_PK.fasta' 

  out <-  parse_input( params$input_tc, params$input_lip  , params$design)
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

  out <-  parse_input( params$input_tc, params$input_lip  , params$design)
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
   ## in case of fasta not seperated by |, we still retrieve correct accesion.
   expect_false(all(is.na(fastaproc$result %>%  distinct(Accession))))
    
      #})
  	
} )


test_that("dialipar_workflowI", {
			params <- list()
         params$input_tc <- test_path('TC_small.parquet')
         params$input_lip <- test_path('LiP_small.parquet')
         params$design <- test_path('SampleAnnotation.txt')
         params$fasta_file <-  test_path('SP_9606_PK.fasta') 

			a_ <-  parse_input( params$input_tc, params$input_lip , params$design)
  	      suppressWarnings({ 
                  qf_base <- create_qfeat_(a_$design, a_$lip, a_$tc, a_$diann_flag)})
        #print(qf_base)
        expect_equal(qf_base$status,0)
        expect_equal(length(qf_base$result),18)
          suppressWarnings({  qf_norm <- normalization_scaling_factor(qf_base$result) })
        expect_equal(dim(assay(qf_norm$result[['proteins_tc']]))[1],39)
        expect_equal(dim(assay(qf_norm$result[['precursors_lip_norm']]))[1],49)	  
} )

test_that("dialipar_workflow_Unpaired", {
			params <- list()
         params$input_tc <- test_path('TC_small.parquet')
         params$input_lip <- test_path('LiP_small.parquet')
         params$design <- test_path('SampleAnnotation.txt')
         params$fasta_file <-  test_path('SP_9606_PK.fasta') 
         params$formula <-  '~  -1 + Treatment'
         params$comparisons <-  c('TreatmentRapa - TreatmentDMSO')
   
			a_ <-  parse_input( params$input_tc, params$input_lip , params$design)
  	      suppressWarnings({ 
                  qf_base <- create_qfeat_(a_$design, a_$lip, a_$tc, a_$diann_flag)})
        #print(qf_base)
        expect_equal(qf_base$status,0)
        expect_equal(length(qf_base$result),18)
        suppressWarnings({  qf_norm <- normalization_scaling_factor(qf_base$result) })
        expect_equal(dim(assay(qf_norm$result[['proteins_tc']]))[1],39)
        expect_equal(dim(assay(qf_norm$result[['precursors_lip_norm']]))[1],49)
         suppressWarnings({ a <-  msqrob_model(pe = qf_norm$result, params = params, layer = 'precursors_lip_norm' )
        b  <-  msqrob_model(pe = a$q_feat, params = params, layer = 'proteins_tc' )
         })
        qf_unpair <- calculate_lip_usage(b$q_feat, i_lip = "precursors_lip_norm", 
                                                   i_tc = "proteins_tc",
                                                   contrasts = base::colnames(b$contr_exp))
         expect_equal(dim(assay(qf_unpair$qf[['proteins_tc']]))[1],39)
        expect_equal(dim(assay(qf_unpair$qf[['precursors_lip_norm']]))[1],49)
         } )