#test on a data set with significant AUCs (AML)
#Could also implement a scheme to select features by their p-values
#get p-value from NES: 1 - pnorm(q = abs(NES), mean = 0, sd = 1)

MRs_per_fold <- function(Risk_group_classification_object = NULL, 
                         k, 
                         c,
                         pipeline_iter = 1,
                         p_value_threshold = 0.05,
                         maximize_diff_null_test = FALSE,
                         use_null = TRUE,
                         savefile = TRUE,
                         loadfile = TRUE){
  
  res = list()
  final_res = list()
  
  if (loadfile == TRUE) Risk_group_classification_object = readRDS('./intermediate_files/training_classification_results.rda')
  
  for (pipeline_iterations in pipeline_iter){
    
    if (is.null(Risk_group_classification_object) == FALSE) k = which(1:max(k) %in% as.numeric(gsub(pattern = 'k', replacement = '', x = sapply(names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC']]), function(i) strsplit(i, '_')[[1]][1])))) 
    
    
    for (clust_iter in c){
      
      list_features = list()
      list_mean_NES = list()
      mat_features = matrix(NA, nrow = 0, ncol = 4, dimnames = list(NULL, c('features', 'votes', 'mean_NES', 'p-val')))
      
      for (fold_iter in k){
        
        #Load the Risk_group_classification_object (specific to the pipeline iteration number, the fold number and the cluster number)
        #Extract the AUC values for the validation and the negative control
        
        if (is.null(Risk_group_classification_object)) Risk_group_classification_object = readRDS(paste0('./intermediate_files/validation_p', pipeline_iterations, '_k', fold_iter, '_c', clust_iter, '.rda'))
        
        list_AUC = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC']][[paste0('k', fold_iter, '_c', clust_iter)]]
        list_AUC_neg_ctrl = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC_neg_ctrl']][[paste0('k', fold_iter, '_c', clust_iter)]]
        if (use_null == TRUE) list_AUC_null = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC_null']][[paste0('k', fold_iter, '_c', clust_iter)]]
        
        vector_means_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, mean) 
        if (use_null == TRUE) vector_means_AUC_null = apply(list_AUC_null, 2, mean)
        
        #Check how many AUCs were computed (if only one, stochastic method will be false)
        if ( length(list_AUC[,1]) == 1 ) { #*** CHECK TO SEE IF THIS TRICK WORKS WHEN THERE IS ONLY ONE RESULT PER ROW (THE COMMA MAY INTRODUCE AN INDEXING ERROR)
          stochastic_method = FALSE
        }else{
          stochastic_method = TRUE  
        }
        
        
        #if true, there are two AUC distributions: one for the negative control and one for the validation. p-value is calculated by a t-test
        #(stochastic_method == FALSE) -- if false, there is one AUC distribution for the negative control and one single point for the negative control, p-value is calculated by computing the AUC's z-score in the negative control's distribution (minimum: 1/nb samples)
        
        #//////////////////////
        #Do not perform t-test to compute p-value even when stochastic method == TRUE
        #t.test can raise an error if there is no variability in the distributions (for example all 1's)
        
        if (stochastic_method == TRUE)  list_AUC = apply(list_AUC, 2, mean)
        # nb_features_tested = ncol(list_AUC)
        # #invert AUCs if systematically under 0.5
        # if (sum(apply(list_AUC, 2, mean) < 0.5)/ncol(list_AUC) > 0.7) { #if AUC is systematically below 0.5, invert it
        #   list_AUC = 1 - list_AUC
        #   list_AUC_neg_ctrl = 1 - list_AUC_neg_ctrl
        #   if (use_null) list_AUC_null = 1 - list_AUC_null
        # }
        # #compute p-values by t-test
        # vector_p_vals = sapply(1:ncol(list_AUC), function(i) t.test(list_AUC[,i], list_AUC_neg_ctrl[,i])$p.value )
        # if (use_null) vector_p_vals_null = sapply(1:ncol(list_AUC), function(i) t.test(list_AUC[,i], list_AUC_null[,i])$p.value )
        # 
        # vector_AUC = apply(list_AUC, 2, mean)
        # vector_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, mean)
        # 
        # if (use_null) vector_AUC_null = apply(list_AUC_null, 2, mean)
        
        stochastic_method = FALSE
        
        if (stochastic_method == FALSE) { #stochastic_method FALSE
          nb_features_tested = length(list_AUC)
          #invert AUC signs if systematically under 0.5
          if (sum(list_AUC < 0.5)/length(list_AUC) > 0.7) { #if AUC is systematically of the wrong sign, invert the signs
            list_AUC = 1 - list_AUC
            list_AUC_neg_ctrl = 1 - list_AUC_neg_ctrl
            if (use_null) list_AUC_null = 1 - list_AUC_null
          }  
          
          vector_AUC = list_AUC
          vector_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, mean)
          if (use_null) vector_AUC_null = apply(list_AUC_null, 2, mean)
          
          #compute p-value by modeling AUC_neg_ctrl as a normal distribution
          vector_sds_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, sd)
          vector_means_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, mean)
          vector_p_vals = sapply(1:length(list_AUC), function(i) pnorm(q = list_AUC[i], mean = vector_means_AUC_neg_ctrl[i], sd = vector_sds_AUC_neg_ctrl[i], lower.tail = F)) 
          
          if (use_null){
            vector_sds_AUC_null = apply(list_AUC_null, 2, sd)
            vector_means_AUC_null = apply(list_AUC_null, 2, mean)
            vector_p_vals_null = sapply(1:length(list_AUC), function(i) pnorm(q = list_AUC[i], mean = vector_means_AUC_null[i], sd = vector_sds_AUC_null[i], lower.tail = F)) 
          }
          
        }
        if (is.na(vector_p_vals[1])){
          vector_p_vals[1] = 0
        }
        
        
        #Condition for selecting the top features: minimum index value when continuity from 1 is broken in:
        #1- indices for significant p-values
        #2- indices for AUCs above negative control
        
        
        
        Top_features_selected = FALSE
        COUNTER = 2
        Top_features = NULL
        nb_features_to_test = length(vector_AUC)
        
        if (maximize_diff_null_test == FALSE){ #will select a number of proteins that are consistently above the negative control and the null model (if available)
          if (use_null == TRUE){
            while (Top_features_selected == FALSE){
              
              #if ( (vector_AUC[COUNTER] >= vector_means_AUC_neg_ctrl[COUNTER]) & (vector_p_vals[COUNTER] < p_value_threshold) & (vector_AUC[COUNTER] >= vector_means_AUC_null[COUNTER]) & (vector_p_vals_null[COUNTER] < p_value_threshold) & (COUNTER <= nb_features_to_test) ){
              if ( (vector_AUC[COUNTER] >= vector_means_AUC_neg_ctrl[COUNTER]) & (vector_AUC[COUNTER] >= vector_means_AUC_null[COUNTER]) & (vector_p_vals_null[COUNTER] < p_value_threshold) & (COUNTER <= nb_features_to_test) ){
                COUNTER = COUNTER + 1  
              } else{
                if (COUNTER > 2){
                  Top_features = 1:(COUNTER - 1)
                  Top_features = names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]])[Top_features]
                }
                Top_features_selected = TRUE
              }
            }
          }else{
            
            while (Top_features_selected == FALSE){
              
              #if ( (vector_AUC[COUNTER] >= vector_means_AUC_neg_ctrl[COUNTER]) && (vector_p_vals[COUNTER] < p_value_threshold) & (COUNTER <= nb_features_to_test) ){
              if ( (vector_AUC[COUNTER] >= vector_means_AUC_neg_ctrl[COUNTER]) & (COUNTER <= nb_features_to_test) ){
                COUNTER = COUNTER + 1  
              } else{
                if (COUNTER > 2){
                  Top_features = 1:(COUNTER - 1)
                  Top_features = names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]])[Top_features]
                }
                Top_features_selected = TRUE
              }
            }  
          }
          
        } else{ #will only maximize the p-value between the null model and the test AUC and select the protein with lowest p-value that is also above the negative control
          
          COUNTER = 2
          ordered_p_vals = order(vector_p_vals_null, decreasing = F)
          while (Top_features_selected == FALSE){
            index = ordered_p_vals[COUNTER]
            
            
            if ( (vector_AUC[index] >= vector_means_AUC_neg_ctrl[index]) & (vector_p_vals_null[index] < p_value_threshold)  ){
              Top_features = 1:index
              Top_features = names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]])[Top_features]
              Top_features_selected = TRUE
            } else{
              COUNTER = COUNTER + 1
              if (COUNTER >= length(ordered_p_vals)) {
                Top_features_selected = TRUE
                Top_features = NULL
              }
            }
          }
        }
        
        
        
        list_features[[fold_iter]] = Top_features
        list_mean_NES[[fold_iter]] = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]][Top_features]
        
        
        
      }#end fold_iter
      
      ##############3
      #Part 2: merging the features from each fold (within the cluster number and the pipeline iteration)
      
      #if no significant features were found, return an empty list
      if (is.null(Top_features)) {
        res[[paste0('iter_', pipeline_iterations)]][['significant_MRs']][[paste0('c', clust_iter)]] = list()
        if (savefile == TRUE) saveRDS(res, paste0('./intermediate_files/significant_MRs_p', pipeline_iterations, '_c', clust_iter, '.rda'))
      } else{
        
        
        mat_features = matrix(data = NA, nrow = length(unique(unlist(list_features))), ncol = 3, dimnames = list(unique(unlist(list_features)), c('features', 'votes', 'mean_NES')) )
        unique_features = unique(unlist(list_features))
        row.names(mat_features) = unique_features
        mat_features[,'features'] = unique_features
        mat_features[,'votes'] = sapply(unique_features, function(i) table(unlist(list_features))[i] )
        
        list_mean_NES = lapply(unique_features, function(i) mean(unlist(list_mean_NES)[which(unlist(list_features) == i)]) )
        list_mean_NES = unlist(list_mean_NES)
        
        mat_features[,'mean_NES'] = list_mean_NES
        mat_features = mat_features[,2:3]
        mat_features = apply(mat_features, c(1,2), function(i) as.numeric(i))
        mat_features = mat_features[order(-mat_features[,'votes'], -mat_features[,'mean_NES']),]
        
        final_ranked_list_features = row.names(mat_features)
        MR_range = 1:length(final_ranked_list_features)
        
        
        
        if (length(final_ranked_list_features) > 1){ 
          #final_feature_selection = names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]])[Top_features]  #find the names of the features
          final_feature_selection = final_ranked_list_features
          res[[paste0('iter_', pipeline_iterations)]][['significant_MRs']][[paste0('c', clust_iter)]] = final_ranked_list_features
          res[[paste0('iter_', pipeline_iterations)]][['Unsorted_top_Features']][[paste0('c', clust_iter)]] = mat_features
          if (savefile == TRUE) saveRDS(res, paste0('./intermediate_files/significant_MRs_p', pipeline_iterations, '_c', clust_iter, '.rda'))
        }
        else{
          res[[paste0('iter_', pipeline_iterations)]][['significant_MRs']][[paste0('c', clust_iter)]] = list()
          if (savefile == TRUE) saveRDS(res, paste0('./intermediate_files/significant_MRs_p', pipeline_iterations, '_c', clust_iter, '.rda'))
        }
        
        
      }#end if (more than 0 features found)
      
      
      
      final_res = c(final_res, res)
      res = list()
      
    }#end clust_iter
  }#end pipeline_iterations
  return(final_res)
}#end function





#validation of the features found in MRs_per_fold in the heldout data set
#for each fold, if there are significant features, use them to classify the heldout data
#cluster heldout data using the same clustering algorithm
#find the cluster corresponding to 
Internal_validation <- function(Risk_group_classification_object = NULL,
                                pipeline_iter = 1,
                                nb_folds = 5,
                                c,
                                classification_algorithm = c('random_forest', 'logistic_regression', 'lda','svm', 'xgboost', 'neural_nets'), 
                                nb_iterations = 100,
                                interactome,
                                max_features_to_test = NULL,
                                mRNA_control = FALSE,
                                compute_null = TRUE,
                                pAUC_range = c(0,1),
                                savefile = TRUE,
                                loadfile = TRUE){
  #normalize heldout data with the entire data set
  #compute viper
  #use the MRs to classify the heldout data
  #keep vector of predictions and vector of truth values for final ROC curve
  require(viper)
  require(randomForest)
  require(ROCR)
  require(DESeq2)
  require(MASS)
  require(e1071) 
  require(pROC)
  
  skip_validation = FALSE
  
  if (loadfile == TRUE){
    Risk_group_classification_object = readRDS('./intermediate_files/folds.rda')
    interactome = Risk_group_classification_object$Experimental_Settings$interactome
  }
  
  if (classification_algorithm == 'random_forest' || classification_algorithm == 'xgboost' || classification_algorithm == 'neural_nets'){
    stochastic_method = TRUE
  } else {
    stochastic_method = FALSE
  }
  
  for (clust_iter in c){
    for (pipeline_iterations in pipeline_iter){
      if (loadfile == TRUE) GEM = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$folds$GEM
      GEM = as.matrix(GEM)
      
      final_ranked_list_features = readRDS(paste0('./intermediate_files/significant_MRs_p', pipeline_iterations, '_c', clust_iter, '.rda'))
      final_ranked_list_features = final_ranked_list_features[[paste0('iter_', pipeline_iterations)]]$significant_MRs[[paste0('c', clust_iter)]]
      
      if (is.null(max_features_to_test) == FALSE) final_ranked_list_features = final_ranked_list_features[1:min(max_features_to_test, length(final_ranked_list_features))]
      
      if (length(final_ranked_list_features) < 2) skip_validation = TRUE
      MR_range = 1:length(final_ranked_list_features)
      
      #Assignment for training patients and test patients sets
      test_patients_group_1 = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$heldout_patients_group_1
      test_patients_group_2 = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$heldout_patients_group_2
      training_patients_group_1 = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$patient_group_1
      training_patients_group_2 = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$patient_group_2
      
      training_patients_group_1 = setdiff(training_patients_group_1, test_patients_group_1)
      training_patients_group_2 = setdiff(training_patients_group_2, test_patients_group_2)
      
      list_res_AUC = list()
      list_res_pred = list()
      list_res_truth = list()
      list_AUC_neg_ctrl = list()
      list_pred_neg_ctrl = list()
      list_truth_neg_ctrl = list()
      
      list_AUC_null = list()
      list_pred_null = list()
      list_truth_null = list()

      if (skip_validation == FALSE){
        # vps_classification = viperSignature(GEM[,c(training_patients_group_1, training_patients_group_2, test_patients_group_1, test_patients_group_2)],
        #                                     GEM[,c(training_patients_group_1, training_patients_group_2)],
        #                                     per = 1000, 
        #                                     method = 'zscore',
        #                                     verbose = F, bootstrap = F)
        # 
        # vps_classification$signature[which(is.na(vps_classification$signature))] <- 0 #to make sure there are no NA values in the signature (would happen if you have a gene with 0 counts in all samples)
        # vps_classification$nullmodel[which(is.na(vps_classification$nullmodel))] <- 0
        # 
        #remove features with no variation (usually they correspond to rows that sum to 0)
        #features with no variation will not efficiently classify
        #invariant features may end up in the negative control, but they will mess up its classification
        
        # vps_classification$signature = as.matrix(vps_classification$signature)
        # vps_classification$signature[which(is.na(vps_classification$signature))] = 0 #needs to be a matrix, if a data.frame "unspecified columns" error will come up
        # vps_classification$signature = vps_classification$signature[-which(apply(vps_classification$signature, 1, sd) == 0),] #remove genes that don't vary
        
        #mRNA_control == TRUE : for classifying using mRNA expression as predictors instead of viper activities
        
        vp_classification = fast_viper(full_matrix = GEM, 
                                       test = c(training_patients_group_1, training_patients_group_2, test_patients_group_1, test_patients_group_2), 
                                       ref = c(training_patients_group_1, training_patients_group_2),
                                       interactome = interactome, 
                                       per = 100, 
                                       return_GES = mRNA_control)
        
        
        if (mRNA_control == TRUE) {
          vp_classification[which(is.na(vp_classification))] <- 0 #consider removing rows with low variance (these features may lead to errors in classification)
          vp_classification = vp_classification[which(apply(vp_classification, 1, sd) != 0),] #remove rows with SD == 0 (potential problem: global SD may not be 0, but SD within subset fold may be 0)
          }
        
        final_ranked_list_features = intersect(final_ranked_list_features, row.names(vp_classification))
      }
      
      
      
      #PREPARING LISTS OF PARTIAL RESULTS
      #partial results (fold_iter and clust_iter - specific) are stored in matrix_res, list_pred_validation and list_truth_validation
      if (stochastic_method == FALSE){
        matrix_res = matrix(data = 0, nrow = 1, ncol = length(final_ranked_list_features))
        list_pred_validation = matrix(data = NA, nrow = 1*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
        list_truth_validation = matrix(data = NA, nrow = 1*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      }else{
        matrix_res = matrix(data = 0, nrow = nb_iterations, ncol = length(final_ranked_list_features))   
        list_pred_validation = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
        list_truth_validation = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      }
      
      
      matrix_AUC_neg_ctrl = matrix(data = 0, nrow = nb_iterations, ncol = length(final_ranked_list_features))
      matrix_pred_neg_ctrl = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      matrix_truth_neg_ctrl = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      
      matrix_AUC_null = matrix(data = 0, nrow = nb_iterations, ncol = length(final_ranked_list_features))
      matrix_pred_null = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      matrix_truth_null = matrix(data = NA, nrow = nb_iterations*length(c(test_patients_group_1, test_patients_group_2)), ncol = length(final_ranked_list_features))
      
      ################################
      
      if (skip_validation == FALSE){
        cat('\n')
        cat(paste0('Classification, pipeline iteration = ', pipeline_iterations, ' c = ', clust_iter, '\n'))
        pb = txtProgressBar(min = 1, max = length(final_ranked_list_features), initial = 1, style = 3)
        
        for (iter_features in 2:length(final_ranked_list_features)){
          setTxtProgressBar(pb, iter_features)
          
          if (compute_null == TRUE){
            mat_group_1_training_null = matrix(NA, nrow = nb_iterations, ncol = length(training_patients_group_1))
            mat_group_2_training_null = matrix(NA, nrow = nb_iterations, ncol = length(training_patients_group_2))
            mat_group_1_validation_null = matrix(NA, nrow = nb_iterations, ncol = length(test_patients_group_1))
            mat_group_2_validation_null = matrix(NA, nrow = nb_iterations, ncol = length(test_patients_group_2))
            
            mat_group_1_training_null = sapply(1:nb_iterations, function(i) sample(x = c(training_patients_group_1, training_patients_group_2), size = length(training_patients_group_1), replace = F))
            mat_group_2_training_null = sapply(1:nb_iterations, function(i) sample(x = c(training_patients_group_1, training_patients_group_2), size = length(training_patients_group_2), replace = F))
            mat_group_1_validation_null = sapply(1:nb_iterations, function(i) sample(x = c(test_patients_group_1, test_patients_group_2), size = length(test_patients_group_1), replace = F))
            mat_group_2_validation_null = sapply(1:nb_iterations, function(i) sample(x = c(test_patients_group_1, test_patients_group_2), size = length(test_patients_group_2), replace = F))
          }
          
          
          for (iter_replicates_validation in 1:nb_iterations){
            
            matrix_row_indices = seq( 1 + (iter_replicates_validation*length(c(test_patients_group_1, test_patients_group_2)) - length(c(test_patients_group_1, test_patients_group_2))) , iter_replicates_validation*length(c(test_patients_group_1, test_patients_group_2)) )
            #only run the validation once if the classification method is NOT stochastic, otherwise run nb_iterations times
            if ( (iter_replicates_validation == 1 & stochastic_method == FALSE) || (stochastic_method == TRUE) 
            ){
             
              classification_res = classification(classification_method = classification_algorithm,
                                                  vp_matrix = vp_classification,
                                                  group_1_training = training_patients_group_1, 
                                                  group_2_training = training_patients_group_2,
                                                  group_1_validation = test_patients_group_1, 
                                                  group_2_validation = test_patients_group_2,
                                                  features = final_ranked_list_features[1:iter_features],
                                                  pAUC_range = pAUC_range)
              
              
              if (is.null(classification_res$AUC) == TRUE) {
                if (stochastic_method == TRUE){
                  matrix_res[iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                  list_pred_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = NA  
                  list_truth_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                }else{
                  matrix_res[1, rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                  list_pred_validation[1:length(c(test_patients_group_1, test_patients_group_2)) , rank(MR_range)[which(MR_range == iter_features)]  ] = NA  
                  list_truth_validation[1:length(c(test_patients_group_1, test_patients_group_2)) , rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                }
                
              }
              else{
                if (stochastic_method == TRUE){
                  matrix_res[iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$AUC
                  list_pred_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$predictions
                  list_truth_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$truth  
                }
                else{
                  matrix_res[1, rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$AUC
                  list_pred_validation[1:length(c(test_patients_group_1, test_patients_group_2)) , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$predictions
                  list_truth_validation[1:length(c(test_patients_group_1, test_patients_group_2)) , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$truth  
                }
                
              }
              
            }
            
            neg_control_features = sample(x = row.names(vp_classification), size = nrow(vp_classification), replace = F) #at each iteration, a random regulator is selected
            
           
            classification_res_neg_ctrl = classification(classification_method = classification_algorithm,
                                                         vp_matrix = vp_classification,
                                                         group_1_training = training_patients_group_1, 
                                                         group_2_training = training_patients_group_2,
                                                         group_1_validation = test_patients_group_1, 
                                                         group_2_validation = test_patients_group_2,
                                                         features = neg_control_features[1:iter_features],
                                                         pAUC_range = pAUC_range) 
            
            
            
            
            if (compute_null == TRUE){
              group_1_training_classification_null = mat_group_1_training_null[, iter_replicates_validation ]
              group_2_training_classification_null = mat_group_2_training_null[, iter_replicates_validation ]
              group_1_validation_classification_null = mat_group_1_validation_null[, iter_replicates_validation ]
              group_2_validation_classification_null = mat_group_2_validation_null[, iter_replicates_validation ]

              
            classification_res_null = classification(classification_method = classification_algorithm,
                                                     vp_matrix = vp_classification,
                                                     group_1_training = group_1_training_classification_null, 
                                                     group_2_training = group_2_training_classification_null,
                                                     group_1_validation = group_1_validation_classification_null, 
                                                     group_2_validation = group_2_validation_classification_null,
                                                     features = final_ranked_list_features[1:iter_features],
                                                     pAUC_range = pAUC_range)   
            }
         
            
            matrix_row_indices = seq( 1 + (iter_replicates_validation*length(c(test_patients_group_1, test_patients_group_2)) - length(c(test_patients_group_1, test_patients_group_2))) , iter_replicates_validation*length(c(test_patients_group_1, test_patients_group_2)) )
            
         
            if (is.null(classification_res$AUC) == TRUE) {
              matrix_AUC_neg_ctrl[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = NA
              if (compute_null) matrix_AUC_null[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = NA
            }
            if (is.null(classification_res$AUC) == FALSE) {
              matrix_AUC_neg_ctrl[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$AUC
              if (compute_null) matrix_AUC_null[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_null$AUC
            }
            
            
            matrix_pred_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$predictions
            matrix_truth_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$truth
            
            if (compute_null == TRUE){
              matrix_pred_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_null$predictions
              matrix_truth_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_null$truth
            }
            
            
          }
        }
      }
      
      
      
      if (length(final_ranked_list_features) == 0 || skip_validation == TRUE) {
        list_res_AUC = NA
        list_pred = NA
        list_truth = NA
        
        list_AUC_neg_ctrl = NA
        list_pred_neg_ctrl = NA
        list_truth_neg_ctrl = NA
        
        if (compute_null){
          list_AUC_null = NA
          list_pred_null = NA
          list_truth_null = NA
        }
        
      }
      else{
        colnames(matrix_res) = final_ranked_list_features
        colnames(list_pred_validation) = final_ranked_list_features
        colnames(list_truth_validation) = final_ranked_list_features
        
        list_res_AUC = matrix_res
        list_res_pred = list_pred_validation
        list_res_truth = list_truth_validation
        
        list_AUC_neg_ctrl = matrix_AUC_neg_ctrl
        list_pred_neg_ctrl = matrix_pred_neg_ctrl
        list_truth_neg_ctrl = matrix_truth_neg_ctrl
        
        if (compute_null){
          list_AUC_null = matrix_AUC_null
          list_pred_null = matrix_pred_null
          list_truth_null = matrix_truth_null
        }
        else{
          list_AUC_null = list()
          list_pred_null = list()
          list_truth_null = list()
        }
        
        # colnames(matrix_AUC_neg_ctrl) = paste0('MRs_',MR_range)
        # colnames(matrix_pred_neg_ctrl) = paste0('MRs_',MR_range)
        # colnames(matrix_truth_neg_ctrl) = paste0('MRs_',MR_range)
      }
      
        res = list()
        res[[paste0('iter_', pipeline_iterations)]]$final_validation[[paste0('c', clust_iter)]] = list(list_res_AUC = list_res_AUC,
                                                                                                        list_res_pred = list_res_pred,
                                                                                                        list_res_truth = list_res_truth, 
                                                                                                        list_AUC_neg_ctrl = list_AUC_neg_ctrl,
                                                                                                        list_pred_neg_ctrl = list_pred_neg_ctrl,
                                                                                                        list_truth_neg_ctrl = list_truth_neg_ctrl,
                                                                                                        list_AUC_null = list_AUC_null,
                                                                                                        list_pred_null = list_pred_null,
                                                                                                        list_truth_null = list_truth_null)
        
        if (savefile == TRUE) saveRDS(res, paste0('./intermediate_files/final_result_p', pipeline_iterations, '_c', clust_iter, '.rda'))

      
    }#pipeline iteration
  }#clust_iter

  return(res)
}



Final_consolidation <- function(path = './intermediate_files', results_file_name = 'results.rda'){
  res = readRDS(paste0(path, '/training_classification_results.rda'))
  
  nb_pipeline_iterations = res$Experimental_Settings$nb_iterations_pipeline
  nb_clusters = res$Experimental_Settings$nb_clusters
  
  for (iter_pipeline in 1:nb_pipeline_iterations){
    for (iter_clust in 1:nb_clusters){
      tmp_res = readRDS(paste0(path, '/final_result_p', iter_pipeline, '_c', iter_clust, '.rda'))
      res[[paste0('iter_', iter_pipeline)]]$test_results[[paste0('Cluster_', iter_clust)]] = tmp_res[[paste0('iter_', iter_pipeline)]]$final_validation[[paste0('c', iter_clust)]]
    }
  }
  saveRDS(res, paste0(path, '/', results_file_name))
  return(res)
}





