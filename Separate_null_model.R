#Compute null model separately
#Can be computed in parallel for each value of k and c, to speed up calculation time

#changes to make:
#name permutation of samples method: null model
#name selection of non-significant TRs: negative control
#implement mRNA_control

separate_null_model <- function(Risk_group_classification_object, 
                                k, 
                                c, 
                                nb_iterations_null = 100, 
                                MR_range,
                                interactome, 
                                top_and_bottom = TRUE, 
                                equilibrate_top_bottom = FALSE,
                                mRNA_control = FALSE,
                                equilibrate_classes = FALSE, #different from folds- equilibrate classes. Here will drop off extra patients from the majority class
                                filename = NULL){
  
  set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
  
  if (length(MR_range) == 1) MR_range = 1:MR_range #input for MR_range must be a list of MRs listed by their ranked position
  
  #classification by significant MRs on permutated samples (null model)
  list_AUC_null = list()
  list_pred_null = list()
  list_truth_null = list()
  #classification by non-significant MRs on actual samples (negative controls)
  list_AUC_neg_ctrl = list()
  list_pred_neg_ctrl = list()
  list_truth_neg_ctrl = list()

  for (fold_iter in k){ #extract the training and validation patients in each fold
    GEM_training = as.matrix(Risk_group_classification_object$list_training_vst[[fold_iter]])
    GEM_validation = as.matrix(Risk_group_classification_object$list_validation_vst[[fold_iter]])
    
    for (clust_iter in c){#extract the training and validation patients in each cluster in each fold
      group_1_training = intersect(names(which(Risk_group_classification_object$list_c_training[[fold_iter]] == clust_iter)), Risk_group_classification_object$patient_group_1)
      group_2_training = intersect(names(which(Risk_group_classification_object$list_c_training[[fold_iter]] == clust_iter)), Risk_group_classification_object$patient_group_2)
      
      group_1_validation = intersect(names(Risk_group_classification_object$list_c_validation[[fold_iter]][which(Risk_group_classification_object$list_c_validation[[fold_iter]] == clust_iter)]) , Risk_group_classification_object$patient_group_1)
      group_2_validation = intersect(names(Risk_group_classification_object$list_c_validation[[fold_iter]][which(Risk_group_classification_object$list_c_validation[[fold_iter]] == clust_iter)]) , Risk_group_classification_object$patient_group_2)
      
      #Checking the conditions to make sure we have enough samples to run viper and classification
      not_enough_samples = FALSE
      if (length(group_1_training) <= 1) not_enough_samples = TRUE #OR
      if (length(group_2_training) <= 1) not_enough_samples = TRUE #OR
      if (length(group_1_validation) == 0){ #AND
        if (length(group_2_validation) == 0) not_enough_samples = TRUE
      } 
      if (sum(length(group_1_training), length(group_2_training)) <= 2) not_enough_samples = TRUE
      
      if (not_enough_samples == TRUE){
        list_AUC_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
        list_pred_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
        list_truth_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
        list_AUC_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
        list_pred_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
        list_truth_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
      }
      
      if (not_enough_samples == FALSE){
        if (length(group_1_training) <= 2 | length(group_2_training) <= 2){
          #Cannot compute a signature with 2 or fewer samples, this workaround will yield the correct result
          #This workaround may be slightly slower, for this reason I don't use it as the default method for computing the viper matrix
          vps = viperSignature(GEM_training[, c(group_1_training, group_2_training)], 
                               GEM_training[, c(group_1_training, group_2_training)], 
                               per = 1000, method = 'zscore', verbose = FALSE)
          
          if (mRNA_control == TRUE) vp = vps$signature
          if (mRNA_control == FALSE) vp = viper(vps, method = 'none', regulon = interactome, verbose = F)
          
          vp = vp[ , group_1_training]
        }
        #normal method for computing viper signature
        #run vipersignature on the actual training samples (not with scrambled labels)
        #null model will check classification accuracy using a genuine MR list on a list with scrambled labels
        if (length(group_1_training) > 2 & length(group_2_training) > 2){
          vps = viperSignature(GEM_training[,group_1_training],
                               GEM_training[,group_2_training], 
                               per = 1000, method = 'zscore', verbose = FALSE)
          
          if (mRNA_control == TRUE) vp = vps$signature
          if (mRNA_control == FALSE) vp = viper(eset = vps$signature, method = 'none', regulon = interactome, verbose = F)
          
        }
        

        if (is.null(dim(vp)) == TRUE){ #if there is only one sample in group 1, the result will be a dimensionless list, not a matrix
          if (top_and_bottom == TRUE & equilibrate_top_bottom == FALSE) stouffer = abs(vp)
          if (top_and_bottom == FALSE) stouffer = vp
          if (equilibrate_top_bottom == TRUE) {
            stouffer_top = sort(vp, decreasing = T)
            stouffer_bottom = sort(vp, decreasing = F)
            stouffer = stouffer_top
            stouffer[seq(2, length(stouffer_bottom), by = 2)] = stouffer_bottom[1:floor(length(stouffer_bottom)/2)]
            names(stouffer)[seq(2, length(stouffer_bottom), by = 2)] = names(stouffer_bottom)[1:floor(length(stouffer_bottom)/2)]
          }
          #will pick "null genes" (the ones that are not statistically significant - with lowest NES)
          stouffer_neg_ctrl = sort(abs(vp), decreasing = F)
        }
        else{ #normal case scenario where vp is a matrix with non-null dimensions
          if (top_and_bottom == TRUE & equilibrate_top_bottom == FALSE) stouffer = sort(apply( vp, 1, function(i) sum(abs(i)) ), decreasing = T)
          if (top_and_bottom == FALSE) stouffer = sort(apply( vp, 1, function(i) sum(i) ), decreasing = T)
          if (equilibrate_top_bottom == TRUE) {
            stouffer_top = sort(apply( vp, 1, function(i) sum(i) ), decreasing = T)
            stouffer_bottom = sort(apply( vp, 1, function(i) sum(i) ), decreasing = F)
            stouffer = stouffer_top
            stouffer[seq(2, length(stouffer_bottom), by = 2)] = stouffer_bottom[1:floor(length(stouffer_bottom)/2)]
            names(stouffer)[seq(2, length(stouffer_bottom), by = 2)] = names(stouffer_bottom)[1:floor(length(stouffer_bottom)/2)]
          }
          #negative control: will pick "null genes" (the ones that are not statistically significant - with lowest NES)
          stouffer_neg_ctrl = sort(abs(apply( vp, 1, function(i) sum(i))), decreasing = F)
        }
      
        #full viper matrix (signature with all patients) to use for classification
        #patient labels are not important, as long as we keep separate the training samples from the validation samples
        #labels will later be shuffled to create the null model
        vps_validation_null = viperSignature(GEM_validation[,c(group_1_training, group_2_training, group_1_validation, group_2_validation)],
                                             GEM_validation[,c(group_1_training, group_2_training)],
                                             per = 1000,
                                             method = 'zscore',
                                             verbose = F)
        
        if (mRNA_control == TRUE) vp_validation_null = vps_validation_null$signature
        if (mRNA_control == FALSE) vp_validation_null = viper( eset = vps_validation_null,
                                                               regulon = interactome,
                                                               method = 'none',
                                                               verbose = F)  
        
        #remove negative controls for which the expression/inferred activity does not vary
        #without this step, some of the selected negative controls would have 0 variance and the random forest classifier would stall - the program will not crash, only stall forever
        #mostly/only a problem with mRNA_control == TRUE
        if (mRNA_control == TRUE) stouffer_neg_ctrl = stouffer_neg_ctrl[  which(apply(vp_validation_null[names(stouffer_neg_ctrl) , ], 1, function(i) length(unique(i))) > (0.7*ncol(vp_validation_null))  )  ]
        

      #initialize the matrix that will have the final AUC results (real MRs, permutated samples)
      matrix_AUC_null = matrix(data = 0, nrow = nb_iterations_null, ncol = length(MR_range))
      matrix_pred_null = matrix(data = NA, nrow = nb_iterations_null*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
      matrix_truth_null = matrix(data = NA, nrow = nb_iterations_null*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
      
      #initialize matrices to store the results for classification with non significant MRs (negative controls)
      matrix_AUC_neg_ctrl = matrix(data = 0, nrow = nb_iterations_null, ncol = length(MR_range))
      matrix_pred_neg_ctrl = matrix(data = NA, nrow = nb_iterations_null*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
      matrix_truth_neg_ctrl = matrix(data = NA, nrow = nb_iterations_null*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
      
      cat('\n')
      cat(paste0('Null model classification, k = ', fold_iter, ' c = ', clust_iter, '\n'))
      pb = txtProgressBar(min = 1, max = nb_iterations_null, initial = 1, style = 3)
      #iterations for the null model, keeping c and k constant
      
      first_iteration = TRUE  
      for (null_iter in 1:nb_iterations_null){
        #time calculation
        if (first_iteration == TRUE) {
          time_init = as.integer(Sys.time()) 
          second_iteration = TRUE
        }
        if (first_iteration == FALSE & second_iteration == TRUE){
          time_final = as.integer(Sys.time()) - time_init 
          cat(paste0('\nTime for one iteration (seconds): ', time_final))
          cat(paste0('\nTime to complete null model (seconds): ', time_final * (nb_iterations_null) ))
          cat(paste0('\nTime to complete null model (minutes): ', round((1/60) * time_final * (nb_iterations_null),2) ))
          cat(paste0('\nTime to complete null model (hours): ', round((1/3600)* time_final * (nb_iterations_null),2) ))
          cat(paste0('\nTime to complete null model (days): ', round((1/86400)* time_final * (nb_iterations_null),2) ))
          cat(paste0('\nCurrent time: ', Sys.time() ))
          cat(paste0('\nEstimated time to complete null model: ', Sys.time() + (time_final * (nb_iterations_null - 1) ), '\n'  ))
          second_iteration = FALSE
        }
        first_iteration = FALSE
        
        setTxtProgressBar(pb, null_iter)
          
        #randomize the training samples, keeping the group sizes the same
        group_1_training_null = sample(x = c(group_1_training, group_2_training), size = length(group_1_training), replace = F)
        group_2_training_null = setdiff(c(group_1_training, group_2_training), group_1_training_null)
        
        
        
        #in case of LOOCV, or any case where there is only 1 validation sample: assign the validation sample randomly in group 1 or 2
        if (length(which(Risk_group_classification_object$list_c_validation[[fold_iter]] == clust_iter)) == 1){ 
          if (sample(x = 1:2, size = 1) == 1){ #if random number between 1 and 2 is == 1, validation sample is assigned to group 1
            group_1_validation_null = names(Risk_group_classification_object$list_c_validation[[fold_iter]][Risk_group_classification_object$list_c_validation[[fold_iter]] == clust_iter])
            group_2_validation_null = NULL
          }
          else{ #if random number between 1 and 2 is == 2, validation sample is assigned to group 2
            group_1_validation_null = NULL
            group_2_validation_null = names(Risk_group_classification_object$list_c_validation[[fold_iter]][Risk_group_classification_object$list_c_validation[[fold_iter]] == clust_iter])
          }
        }
        #if there is more than 1 sample in the validation set (not LOOCV)
        else{ 
          group_1_validation_null = sample(x = c(group_1_validation, group_2_validation), size = length(group_1_validation), replace = F)
          group_2_validation_null = setdiff(c(group_1_validation, group_2_validation) , group_1_validation_null)
        }
        
            
        ########################################################################################################################################
        ########################################################################################################################################  
            
          #COMPUTE AUC FOR CLASSIFICATION BETWEEN NULL GROUPS 1 AND 2 WITH NULL LIST MRs
          
          features = names(stouffer)
          features_neg_ctrl = names(stouffer_neg_ctrl)
          
          if (MR_range[1] == 1){
            features_range = MR_range[2:length(MR_range)]
          }
          else{
            features_range = MR_range
          }
          

          #final nested loop (4 loops!) - loop for k, loop for c, loop for iterations null, loop for iteration features
          
          for (iter_features in features_range){
            
            if (equilibrate_classes == TRUE){
              if (length(group_1_training) > length(group_2_training)) group_1_training = sample(group_1_training, length(group_2_training), F)
              if (length(group_2_training) > length(group_1_training)) group_2_training = sample(group_2_training, length(group_1_training), F)
              if (length(group_1_validation) > length(group_2_validation)) group_1_validation = sample(group_1_validation, length(group_2_validation), F)
              if (length(group_2_validation) > length(group_1_validation)) group_2_validation = sample(group_2_validation, length(group_1_validation), F)
              
              if (length(group_1_training_null) > length(group_2_training_null)) group_1_training = sample(group_1_training_null, length(group_2_training_null), F)
              if (length(group_2_training_null) > length(group_1_training_null)) group_2_training = sample(group_2_training_null, length(group_1_training_null), F)
              if (length(group_1_validation_null) > length(group_2_validation_null)) group_1_validation = sample(group_1_validation_null, length(group_2_validation_null), F)
              if (length(group_2_validation_null) > length(group_1_validation_null)) group_2_validation = sample(group_2_validation_null, length(group_1_validation_null), F)
            }
            
            
            #classification of permutated samples with significant MRs
            rf_res = Random_forest_classification(vp_matrix = vp_validation_null,
                                                  group_1_training = group_1_training_null, group_2_training = group_2_training_null,
                                                  group_1_validation = group_1_validation_null, group_2_validation = group_2_validation_null,
                                                  features = features[1:iter_features])
            
            
            #classification of non-permutated samples with non-significant MRs
            rf_res_neg_ctrl = Random_forest_classification(vp_matrix = vp_validation_null,
                                                  group_1_training = group_1_training, group_2_training = group_2_training,
                                                  group_1_validation = group_1_validation, group_2_validation = group_2_validation,
                                                  features = sample(x = features_neg_ctrl, size = iter_features, replace = F))
            
            #features = features_neg_ctrl[1:iter_features])
            
            matrix_row_indices = seq( 1 + (null_iter*length(c(group_1_validation_null, group_2_validation_null)) - length(c(group_1_validation_null, group_2_validation_null))) , null_iter*length(c(group_1_validation_null, group_2_validation_null)) )
            
            #rank(MR_range): will sort the list of indices corresponding to the MR range in the input
            if (is.null(rf_res$AUC) == TRUE) {
              matrix_AUC_null[ null_iter, rank(MR_range)[which(MR_range == iter_features)] ] = NA
              matrix_AUC_neg_ctrl[ null_iter, rank(MR_range)[which(MR_range == iter_features)] ] = NA
            }
            else{
              matrix_AUC_null[ null_iter, rank(MR_range)[which(MR_range == iter_features)] ] = rf_res$AUC
              matrix_AUC_neg_ctrl[ null_iter, rank(MR_range)[which(MR_range == iter_features)] ] = rf_res_neg_ctrl$AUC
            }
            
            
            matrix_pred_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = rf_res$predictions
            matrix_truth_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = rf_res$truth
            
            matrix_pred_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = rf_res_neg_ctrl$predictions
            matrix_truth_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = rf_res_neg_ctrl$truth
            
          }#end iter_features
          
          colnames(matrix_AUC_null) = paste0('MRs_',MR_range)
          colnames(matrix_pred_null) = paste0('MRs_',MR_range)
          colnames(matrix_truth_null) = paste0('MRs_',MR_range)
          
          list_AUC_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_AUC_null
          list_pred_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_pred_null
          list_truth_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_truth_null
          
          colnames(matrix_AUC_neg_ctrl) = paste0('MRs_',MR_range)
          colnames(matrix_pred_neg_ctrl) = paste0('MRs_',MR_range)
          colnames(matrix_truth_neg_ctrl) = paste0('MRs_',MR_range)
          
          list_AUC_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_AUC_neg_ctrl
          list_pred_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_pred_neg_ctrl
          list_truth_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_truth_neg_ctrl
          
        }

      }#end null_iter
    }#end for (clust_iter in c)
  }#end for (fold_iter in k)
  
  res = list(list_AUC_null = list_AUC_null,
             list_pred_null = list_pred_null, 
             list_truth_null = list_truth_null,
             list_AUC_neg_ctrl = list_AUC_neg_ctrl,
             list_pred_neg_ctrl = list_pred_neg_ctrl,
             list_truth_neg_ctrl = list_truth_neg_ctrl)
  
  if (is.null(filename) == FALSE) saveRDS(object = res, file = filename)
  
  return(res)
}



for (k_iter in 1:3){
  for (c_iter in 1:2){
    test_null <- separate_null_model(readRDS(paste0('./test/vst_test_k3c2_k', k_iter, '.rda')), k = k_iter, c = c_iter, nb_iterations_null = 25, MR_range = 1:10, interactome = pruneRegulon(network_vst), top_and_bottom = T, equilibrate_top_bottom = T, mRNA_control = F, filename = paste0('./test/null_test_k3c2_k', k_iter,'c', c_iter, '_MR_', 1, '_', 10, '.rda' ))
    test_null <- separate_null_model(readRDS(paste0('./test/vst_test_k3c2_k', k_iter, '.rda')), k = k_iter, c = c_iter, nb_iterations_null = 25, MR_range = 11:20, interactome = pruneRegulon(network_vst), top_and_bottom = T, equilibrate_top_bottom = T, mRNA_control = F, filename = paste0('./test/null_test_k3c2_k', k_iter,'c', c_iter, '_MR_', 11, '_', 20, '.rda' ))
  }
}

test_null <- separate_null_model(test_vst_vp, k = 2, c = 1, nb_iterations_null = 25, MR_range = 1:100, interactome = pruneRegulon(network_vst), top_and_bottom = T, equilibrate_top_bottom = T, mRNA_control = F)

test_null <- separate_null_model(readRDS(paste0('./test/vst_test_k3c2_k', 2, '.rda')), k = 2, c = 1, nb_iterations_null = 25, MR_range = 1:10, interactome = pruneRegulon(network_vst), top_and_bottom = T, equilibrate_top_bottom = T, mRNA_control = F, filename = paste0('./test/null_test_k3c2_k', k_iter,'c', c_iter, '_MR_', 1, '_', 10, '.rda' ))

plot(apply(test_null$list_AUC_null$k2_c1, 2, zscore), type = 'l', ylim = c(0,1))
lines(x = 1:100, y = apply(test_null$list_AUC_neg_ctrl$k2_c1, 2, 'zscore'))



null_NC_AAML0531_R_NR_T_LOexp1_k10c1_k4_MR_1_50 <- separate_null_model(Risk_group_classification_object = readRDS('./incomplete_data/vst_NC_AAML0531_R_NR_T_LOexp1_k10c1_k4.rda'), 
                                                                       k = 4, 
                                                                       c = 1, 
                                                                       nb_iterations_null = 100, 
                                                                       MR_range = 1:50, 
                                                                       interactome = pruneRegulon(network_AAML0531_TF_coTF), 
                                                                       top_and_bottom = F, 
                                                                       equilibrate_top_bottom = F,
                                                                       mRNA_control = F, 
                                                                       filename = './incomplete_data/null_NC_AAML0531_R_NR_T_LOexp1_k10c1_k4_MR_1_50.rda')


test_null <- separate_null_model(Risk_group_classification_object = test_vst, k = 2, c = 2, MR_range = 1:50, interactome = pruneRegulon(network_TF_coTF), top_and_bottom = F, equilibrate_top_bottom = F, mRNA_control = F, nb_iterations_null = 100)

test_val <- separate_validation(Risk_group_classification_object = test_vst, k = 2, c = 2, MR_range = 1:50, interactome = pruneRegulon(network_TF_coTF), top_and_bottom = F, equilibrate_top_bottom = F, mRNA_control = F, nb_iterations = 100)
vp_test <- viper(test_null$signature, pruneRegulon(network_TF_coTF), method = 'none', verbose = T)


