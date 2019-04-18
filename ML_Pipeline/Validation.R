#FUNCTION TO COMPUTE VALIDATION AUC WITH INCREASING NUMBERS OF MRs

#Risk_group_classification_object: either a fully consolidated vst and clustered object or a partial object (as long as the object has data for the correct k, c and pipeline_iteration, the function will work) -- this object is the ouput of the function "Initial_normalization_clust()"
#k: single or multiple value(s) of k for the computation. can enter a single number (i.e "3" in a 5-fold CV) or multiple numbers (i.e "1:5" in a 5-fold CV)
#c: a single or multiple value(s) of c for the computation. like with k, can enter a single value like "2" or a range like "1:3" or a series like "c(1,3,5:7)"
#pipeline_iteration: iteration in the pipeline (for testing), can be a single value or a series of values (e.g.: 1 or 1:10 or 4:8 or c(2,3,9))
#nb_iterations: number of times to repeat random forest to get a distribution of AUC values
#MR_range: range of ranked MRs to test. You can break up the computation into smaller computations (1:25, 26:50, 51:75, 76:100) vs (1:100)
#interactome: output of the aracne2regulon function. should wrap pruneRegulon() around the interactome
#top_and_bottom: if TRUE, will rank list of MRs this way: top MRs are the ones with largest absolute NES values, if FALSE, top MRs are the ones with largest NES values (large negative NESes are not at the top of the list)
#equilibrate_top_bottom: if TRUE, will rank the list of MRs this way: first MR from top (largest positive NES), second MR from bottom (largest negative NES), second MR from top (second-largest positive NES), second MR from bottom (second-largest negative NES), etc.
#mRNA_control: if TRUE, will compute the pipeline on an mRNA signature not a viper matrix
#filename: will save the result in a file. filename should follow a specific naming scheme: example: prefix_kc_k_iter_c_iter_MRrange_1_100.rda / HiLoRisk_k3c2_k2c1_MRrange_1_100.rda

#need to change this function to take in a range of different classification methods, random forests being one of them
#change all references of 'rf_' to "ClassResult_", same output from other classification methods, if possible
#Implement Parallelization (divide MRs into bins of equal sizes, randomly distribute across cores)

#consider implementing an alternative, optional feature selection scheme: pick most significant regulators, then iterate through all regulators to find the best, then the second best, etc

#subdivide validation into two functions:
#first part: will prepare the data for the validation (make sure there is enough samples, organize the training and validation data, compute the signatures and the viper activity)
#second part is the validation with the results acquisition --- same function can be called for running validation, negative control, null model and final inter-iteration validation
validation <- function(Risk_group_classification_object = NULL, 
                       k, 
                       c,
                       pipeline_iter = 1,
                       classification_algorithm = c('random_forest', 'logistic_regression', 'lda','svm', 'xgboost', 'neural_nets'),
                       nb_iterations = 100, #number of times classification AUC will be computed for each combination of features
                       MR_range = 1:50, #range of ranked MRs for which to compute AUC
                       interactome = NULL,
                       GEM = NULL,
                       top_and_bottom = TRUE, 
                       equilibrate_top_bottom = FALSE,
                       mRNA_control = FALSE, 
                       random_negative_control = TRUE,
                       compute_null = TRUE,
                       equilibrate_classes = FALSE, #different from folds- equilibrate classes. Here will drop off extra patients from the majority class
                       pAUC_range = c(0,1),
                       save_path = './intermediate_files',
                       savefile = TRUE,
                       loadfile = TRUE){
  
  require(viper)
  require(randomForest)
  require(ROCR)
  require(DESeq2)
  require(MASS)
  require(e1071)
  require(pROC)
  
  
  set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
  
  if (is.null(Risk_group_classification_object) || loadfile == TRUE){
    load_risk_group_classification_file = TRUE #will load an rda object with the risk group classification object
  } else {
    load_risk_group_classification_file = FALSE 
  }
  
  if (classification_algorithm == 'random_forest' | classification_algorithm == 'xgboost'){
    stochastic_method = TRUE
  }
  else {
    stochastic_method = FALSE
  }
  
  val_res = list()
  
  if (length(MR_range) == 1) MR_range = 1:MR_range
  
  for (pipeline_iterations in pipeline_iter){
    message('Computing iteration ', pipeline_iterations)
    
    list_MRs = list() #list of features for classification
    
    list_AUC = list() #list of AUC values for each iteration (the mean and standard deviation will be taken to analyze the results)
    list_pred = list() #list of probabilities that each patient in validation cohort belongs to group_1
    list_truth = list() #list of truth values (1 corresponds to group_1, 0 corresponds to group_2)
    
    list_AUC_neg_ctrl = list()
    list_pred_neg_ctrl = list()
    list_truth_neg_ctrl = list()
    
    if (compute_null == TRUE){
      list_AUC_null = list()
      list_pred_null = list()
      list_truth_null = list()
    }
    
    if (load_risk_group_classification_file == TRUE){
      folds_object = readRDS('./intermediate_files/folds.rda')
      interactome = folds_object$Experimental_Settings$interactome
      GEM = folds_object[[paste0('iter_',pipeline_iterations)]]$folds$GEM  
    }
    
    
    #iteration through each value of k
    for (fold_iter in k){
      message('Computing fold  ', fold_iter)
      
      #two possible inputs: in the function call or by loading a data file
      if (load_risk_group_classification_file == TRUE) Risk_group_classification_object = readRDS(paste0(save_path, '/norm_clust_p', pipeline_iterations, '_k', fold_iter, '.rda'))

      
      #Gene expression matrices (vst-normalized) are fold-specific
      GEM_training = as.matrix(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_training_vst[[fold_iter]])
      GEM_validation = as.matrix(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_validation_vst[[fold_iter]])
      
      #iteration through each value of c
      for (clust_iter in c){
        #training and validation groups are cluster-specific
        group_1_training = intersect( names(which(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_c_training[[fold_iter]] == clust_iter)) , folds_object[[paste0('iter_', pipeline_iterations)]]$patient_group_1 )
        group_2_training = intersect( names(which(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_c_training[[fold_iter]] == clust_iter)) , folds_object[[paste0('iter_', pipeline_iterations)]]$patient_group_2 )
        group_1_validation = intersect(names(which(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_c_validation[[fold_iter]] == clust_iter)), folds_object[[paste0('iter_', pipeline_iterations)]]$patient_group_1)
        group_2_validation = intersect(names(which(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]]$list_c_validation[[fold_iter]]== clust_iter)), folds_object[[paste0('iter_', pipeline_iterations)]]$patient_group_2)
        
        if (compute_null == TRUE){
          mat_group_1_training_null = matrix(NA, nrow = nb_iterations, ncol = length(group_1_training))
          mat_group_2_training_null = matrix(NA, nrow = nb_iterations, ncol = length(group_2_training))
          mat_group_1_validation_null = matrix(NA, nrow = nb_iterations, ncol = length(group_1_validation))
          mat_group_2_validation_null = matrix(NA, nrow = nb_iterations, ncol = length(group_2_validation))
          
          mat_group_1_training_null = sapply(1:nb_iterations, function(i) sample(x = c(group_1_training, group_2_training), size = length(group_1_training), replace = F))
          mat_group_2_training_null = sapply(1:nb_iterations, function(i) sample(x = c(group_1_training, group_2_training), size = length(group_2_training), replace = F))
          mat_group_1_validation_null = sapply(1:nb_iterations, function(i) sample(x = c(group_1_validation, group_2_validation), size = length(group_1_validation), replace = F))
          mat_group_2_validation_null = sapply(1:nb_iterations, function(i) sample(x = c(group_1_validation, group_2_validation), size = length(group_2_validation), replace = F))
        }
        
##################################
        #partial results (fold_iter and clust_iter - specific) are stored in matrix_res, list_pred_validation and list_truth_validation
        
        if (equilibrate_classes == TRUE) stochastic_method = TRUE
        
        if (stochastic_method == FALSE){
          matrix_res = matrix(data = 0, nrow = 1, ncol = length(MR_range))
          list_pred_validation = matrix(data = NA, nrow = 1*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
          list_truth_validation = matrix(data = NA, nrow = 1*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        }else{
          matrix_res = matrix(data = 0, nrow = nb_iterations, ncol = length(MR_range))   
          list_pred_validation = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
          list_truth_validation = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        }
        
        
        matrix_AUC_neg_ctrl = matrix(data = 0, nrow = nb_iterations, ncol = length(MR_range))
        matrix_pred_neg_ctrl = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        matrix_truth_neg_ctrl = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        
        matrix_AUC_null = matrix(data = 0, nrow = nb_iterations, ncol = length(MR_range))
        matrix_pred_null = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        matrix_truth_null = matrix(data = NA, nrow = nb_iterations*length(c(group_1_validation, group_2_validation)), ncol = length(MR_range))
        
        
        #LIST TOP MRs
        #If there aren't enough samples, validation will be skipped
        not_enough_samples = FALSE
        if (length(group_1_training) <= 1) not_enough_samples = TRUE #OR
        if (length(group_2_training) <= 1) not_enough_samples = TRUE #OR
        if (length(group_1_validation) == 0){ #AND
          if (length(group_2_validation) == 0) not_enough_samples = TRUE
        } 
        if (sum(length(group_1_training), length(group_2_training)) <= 2) not_enough_samples = TRUE
        if (not_enough_samples == FALSE){
          if (length(group_1_training) <= 2 | length(group_2_training) <= 2){
            #Cannot compute a signature with 2 or fewer samples, this workaround will yield the correct result
            #This workaround may be slightly slower, for this reason I don't use it as the default method for computing the viper matrix
            
            vp = fast_viper(full_matrix = GEM_training[, c(group_1_training, group_2_training)], 
                test = c(group_1_training, group_2_training), 
                ref = c(group_1_training, group_2_training), 
                interactome = interactome,
                return_GES = mRNA_control)
            
            if (mRNA_control == TRUE) {
              vp[which(is.na(vp))] <- 0
              vp = vp[which(apply(vp, 1, sd) != 0),]
            }
            
            vp = vp[ , group_1_training]
          }
          else{ #enough samples
            
            vp = fast_viper(full_matrix = GEM_training, 
                            test = group_1_training, 
                            ref = group_2_training, 
                            interactome = interactome, 
                            return_GES = mRNA_control)
            
            if (mRNA_control == TRUE){
              vp[which(is.na(vp))] <- 0
              vp = vp[which(apply(vp, 1, sd) != 0),]
            } 
            
          }
          
          
          
          if (top_and_bottom == TRUE & equilibrate_top_bottom == FALSE) stouffer = sort(abs(apply( vp, 1, sum )), decreasing = T)
          if (top_and_bottom == FALSE) stouffer = sort(apply( vp, 1, sum ), decreasing = T)
          
          #create two lists: one sorted with decreasing = TRUE and on with decreasing = FALSE
          #merge the two lists staggering them by 1 index
          if (equilibrate_top_bottom == TRUE) {
            stouffer_top = sort(apply( vp, 1, sum ), decreasing = T)
            stouffer_bottom = sort(apply( vp, 1, sum ), decreasing = F)
            stouffer = stouffer_top
            stouffer[seq(2, length(stouffer_bottom), by = 2)] = stouffer_bottom[1:floor(length(stouffer_bottom)/2)]
            names(stouffer)[seq(2, length(stouffer_bottom), by = 2)] = names(stouffer_bottom)[1:floor(length(stouffer_bottom)/2)]
          }
          
          non_sig_reg <- names(sort(abs(apply(vp, 1, sum)), decreasing = F)) #non-significant "regulators" (protein activies with |NES|'s close to 0)
          
          list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]] = stouffer
          
          #VALIDATION ########################################
          
          vp_c = fast_viper(full_matrix = cbind(GEM_training, GEM_validation), 
                            test = c(group_1_training, group_2_training, group_1_validation, group_2_validation), 
                            ref = c(group_1_training, group_2_training), 
                            interactome = interactome, 
                            return_GES = mRNA_control)
          
          
          if (mRNA_control == TRUE) {
            vp_c[which(is.na(vp_c))] <- 0
            vp_c = vp_c[which(apply(vp_c, 1, sd) != 0),]
            }

          vp_training_c = vp_c[ , c(group_1_training, group_2_training)]
          vp_validation_c = vp_c[ , c(group_1_validation, group_2_validation)]
          
          features = names(list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]])
          
          if (MR_range[1] == 1){
            features_range = MR_range[2:length(MR_range)]
          }
          else{
            features_range = MR_range
          }
          
          cat('\n')
          cat(paste0('Classification, k = ', fold_iter, ' c = ', clust_iter, '\n'))
          pb = txtProgressBar(min = features_range[1], max = features_range[length(features_range)], initial = features_range[1], style = 3)
          
          
          #case: if the features to be tested are not in the classification viper matrix: skip classification
          if (sum(features_range > nrow(vp_c)) > 0){ #if some feature indices are above the number of features in the viper matrix
            features_range = min(features_range) : nrow(vp_c) 
            if (length(features_range) == 1) stop('ERROR: MR range out of range of the viper matrix. Set MR_range option so that it does not go beyond ', nrow(vp_c))
            if (features_range[1] > features_range[2]) stop('ERROR: MR range out of range of the viper matrix. Set MR_range option so that it does not go beyond ', nrow(vp_c))
          }
          
          for (iter_features in features_range){

            setTxtProgressBar(pb, iter_features)
            
            for (iter_replicates_validation in 1:nb_iterations){
            
                group_1_training_classification = group_1_training
                group_2_training_classification = group_2_training
                group_1_validation_classification = group_1_validation
                group_2_validation_classification = group_2_validation
                
                #IMPLEMENT EQUILIBRATE CLASSES -- NOTE THAT list_pred_res and list_pred_truth WILL HAVE DIFFERENT LENGTHS AT EACH ITERATION
                #May not need to equilibrate classes, as it is already implemented at the patient subsetting step
                # if(equilibrate_classes == TRUE){
                #   equilibration_type = sample(c('oversample_minority', 'undersample_majority'), 1)
                #   if (equilibration_type == 'oversample_minority')
                # }
          
                
                if (compute_null == TRUE){
                  group_1_training_classification_null = mat_group_1_training_null[, iter_replicates_validation ]
                  group_2_training_classification_null = mat_group_2_training_null[, iter_replicates_validation ]
                  group_1_validation_classification_null = mat_group_1_validation_null[, iter_replicates_validation ]
                  group_2_validation_classification_null = mat_group_2_validation_null[, iter_replicates_validation ]
                }
                
                
              matrix_row_indices = seq( 1 + (iter_replicates_validation*length(c(group_1_validation_classification, group_2_validation_classification)) - length(c(group_1_validation_classification, group_2_validation_classification))) , iter_replicates_validation*length(c(group_1_validation_classification, group_2_validation_classification)) )
              #only run the validation once if the classification method is NOT stochastic, otherwise run nb_iterations times
              if ( (iter_replicates_validation == 1 & stochastic_method == FALSE) ||
                   (stochastic_method == TRUE) 
              ){
                classification_res = classification(classification_method = classification_algorithm,
                                        vp_matrix = vp_c,
                                        group_1_training = group_1_training_classification, 
                                        group_2_training = group_2_training_classification,
                                        group_1_validation = group_1_validation_classification, 
                                        group_2_validation = group_2_validation_classification,
                                        features = features[1:iter_features],
                                        pAUC_range = pAUC_range)
                
        
                
                if (is.null(classification_res$AUC) == TRUE) {
                  if (stochastic_method == TRUE){
                    matrix_res[iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                    list_pred_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = NA  
                    list_truth_validation[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                  }else{
                    matrix_res[1, rank(MR_range)[which(MR_range == iter_features)]  ] = NA
                    list_pred_validation[1:length(c(group_1_validation_classification, group_2_validation_classification)) , rank(MR_range)[which(MR_range == iter_features)]  ] = NA  
                    list_truth_validation[1:length(c(group_1_validation_classification, group_2_validation_classification)) , rank(MR_range)[which(MR_range == iter_features)]  ] = NA
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
                    list_pred_validation[1:length(c(group_1_validation_classification, group_2_validation_classification)) , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$predictions
                    list_truth_validation[1:length(c(group_1_validation_classification, group_2_validation_classification)) , rank(MR_range)[which(MR_range == iter_features)]  ] = classification_res$truth  
                  }
                  
                }
            
                
              }

              
              if (random_negative_control == TRUE) neg_control_features = sample(x = non_sig_reg, size = length(non_sig_reg), replace = F) #at each iteration, a random regulator is selected
              if (random_negative_control == FALSE) neg_control_features = non_sig_reg
              
              classification_res_neg_ctrl = classification(classification_method = classification_algorithm,
                                                           vp_matrix = vp_c,
                                                           group_1_training = group_1_training_classification, 
                                                           group_2_training = group_2_training_classification,
                                                           group_1_validation = group_1_validation_classification, 
                                                           group_2_validation = group_2_validation_classification,
                                                           features = neg_control_features[1:iter_features],
                                                           pAUC_range = pAUC_range) #sample(features, iter_features, F)) 
             
              
              
              if (compute_null == TRUE){
                classification_null = classification(classification_method = classification_algorithm,
                                                     vp_matrix = vp_c,
                                                     group_1_training = group_1_training_classification_null, 
                                                     group_2_training = group_2_training_classification_null,
                                                     group_1_validation = group_1_validation_classification_null, 
                                                     group_2_validation = group_2_validation_classification_null,
                                                     features = features[1:iter_features],
                                                     pAUC_range = pAUC_range) 
              }
          
              
              
              if (is.null(classification_res$AUC) == TRUE) {
                matrix_AUC_neg_ctrl[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = NA
                if (compute_null) matrix_AUC_null[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = NA
                }
              if (is.null(classification_res$AUC) == FALSE) {
                matrix_AUC_neg_ctrl[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$AUC
                if (compute_null) matrix_AUC_null[ iter_replicates_validation, rank(MR_range)[which(MR_range == iter_features)] ] = classification_null$AUC
                }
            

              matrix_pred_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$predictions
              matrix_truth_neg_ctrl[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_res_neg_ctrl$truth
              
              if (compute_null == TRUE){
                matrix_pred_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_null$predictions
                matrix_truth_null[matrix_row_indices , rank(MR_range)[which(MR_range == iter_features)] ] = classification_null$truth
              }
              
              
            }#END iter_replicates_validation (replicates maintaining # of features (top MRs) constant)
          }#END iter_features
        }#END classification / validation
        
        if (not_enough_samples == TRUE){
          list_AUC[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          list_pred[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          list_truth[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]] = NA
          
          list_AUC_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          list_pred_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          list_truth_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          
          if (compute_null){
            list_AUC_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
            list_pred_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
            list_truth_null[[paste0('k', fold_iter, '_c', clust_iter)]] = NA
          }
          
        }
        else{
          colnames(matrix_res) = names(list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]][MR_range])
          colnames(list_pred_validation) = names(list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]][MR_range])
          colnames(list_truth_validation) = names(list_MRs[[paste0('k',fold_iter,'_c',clust_iter)]][MR_range])
          
          list_AUC[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_res
          list_pred[[paste0('k', fold_iter, '_c', clust_iter)]] = list_pred_validation
          list_truth[[paste0('k', fold_iter, '_c', clust_iter)]] = list_truth_validation
          
          
          colnames(matrix_AUC_neg_ctrl) = paste0('MRs_',MR_range)
          colnames(matrix_pred_neg_ctrl) = paste0('MRs_',MR_range)
          colnames(matrix_truth_neg_ctrl) = paste0('MRs_',MR_range)
          
          list_AUC_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_AUC_neg_ctrl
          list_pred_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_pred_neg_ctrl
          list_truth_neg_ctrl[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_truth_neg_ctrl
          
          
          if (compute_null){
            list_AUC_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_AUC_null
            list_pred_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_pred_null
            list_truth_null[[paste0('k', fold_iter, '_c', clust_iter)]] = matrix_truth_null
          } 
          else{
            list_AUC_null = list()
            list_pred_null = list()
            list_truth_null = list()
          }
          
        }
          
          if (savefile == TRUE){
            res = list()
            res[[paste0('iter_', pipeline_iterations)]] = list(list_MRs = list_MRs,
                                                          list_AUC = list_AUC,
                                                          list_pred = list_pred,
                                                          list_truth = list_truth, 
                                                          list_AUC_neg_ctrl = list_AUC_neg_ctrl,
                                                          list_pred_neg_ctrl = list_pred_neg_ctrl,
                                                          list_truth_neg_ctrl = list_truth_neg_ctrl,
                                                          list_AUC_null = list_AUC_null,
                                                          list_pred_null = list_pred_null,
                                                          list_truth_null = list_truth_null)
            
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$classification_algorithm = classification_algorithm
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$mRNA_control = mRNA_control
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$random_negative_control = random_negative_control
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$equilibrate_classes = equilibrate_classes
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$top_and_bottom = top_and_bottom
            res[[paste0('iter_', pipeline_iterations)]]$Classification_Settings$equilibrate_top_bottom = equilibrate_top_bottom
            
            saveRDS(res, paste0('./intermediate_files/validation_p', pipeline_iterations, '_k', fold_iter, '_c', clust_iter, '.rda'))
            
          }
        
###########################
      }#end clust_iter
    }#end fold_iter (k)
    
    res = list(list_MRs = list_MRs,
               list_AUC = list_AUC,
               list_pred = list_pred,
               list_truth = list_truth, 
               list_AUC_neg_ctrl = list_AUC_neg_ctrl,
               list_pred_neg_ctrl = list_pred_neg_ctrl,
               list_truth_neg_ctrl = list_truth_neg_ctrl,
               list_AUC_null = list_AUC_null,
               list_pred_null = list_pred_null,
               list_truth_null = list_truth_null)
    
    
    val_res[[paste0('iter_', pipeline_iterations)]] = res
  }#end pipeline_iteration
  
  
  
  return(val_res)
  
}
 



#implement a function that merges all the normalization/clustering files and the validation files
#try to implement a function that detects the number of pipeline iterations, number of folds and number of clusters automatically
#because the number of pipeline iterations might not be known (need to divide 1/fraction of heldout samples)
#the number of clusters might not be known either because (e.g. if you set clusters = 2, but MClust only finds 1 cluster, the setting nb_clusters will be set incorrectly to 2)

merge_norm_clust_validation_files <- function(path = './intermediate_files', remove_intermediate_files = FALSE){
  res = readRDS(paste0(path, '/folds.rda'))
  
  nb_pipeline_iterations = res$Experimental_Settings$nb_iterations_pipeline
  nb_folds = res$Experimental_Settings$nb_folds
  nb_clusters = res$Experimental_Settings$nb_clusters
  
  for (iter_pipeline in 1:nb_pipeline_iterations){
    for (iter_folds in 1:nb_folds){
      
      tmp_norm = tryCatch({
        readRDS(paste0(path, '/norm_clust_p', iter_pipeline, '_k', iter_folds, '.rda'))
      }, error = function(e) {
        message(paste0('Warning: ', path, '/norm_clust_p', iter_pipeline, '_k', iter_folds, '.rda does not exist'))
        return(NULL)
      })
        
      res[[paste0('iter_', iter_pipeline)]]$list_training_vst = tmp_norm[[paste0('iter_', iter_pipeline)]]$list_training_vst
      res[[paste0('iter_', iter_pipeline)]]$list_validation_vst = tmp_norm[[paste0('iter_', iter_pipeline)]]$list_validation_vst
      res[[paste0('iter_', iter_pipeline)]]$list_c_training = tmp_norm[[paste0('iter_', iter_pipeline)]]$list_c_training
      
      for (iter_clust in 1:nb_clusters){
        tmp_validation = tryCatch({
          readRDS(paste0(path, '/validation_p', iter_pipeline, '_k', iter_folds, '_c', iter_clust, '.rda'))
        }, error = function(e) {
          message(paste0('Warning: ', path, '/validation_p', iter_pipeline, '_k', iter_folds, '_c', iter_clust, '.rda does not exist'))
          return(NULL)
        })
        
        res[[paste0('iter_', iter_pipeline)]]$list_MRs[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_MRs[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_AUC[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_AUC[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_pred[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_pred[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_truth[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_truth[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_AUC_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_AUC_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_pred_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_pred_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_truth_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_truth_neg_ctrl[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_AUC_null[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_AUC_null[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_pred_null[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_pred_null[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$list_truth_null[[paste0('k', iter_folds, '_c', iter_clust)]] = tmp_validation[[paste0('iter_', iter_pipeline)]]$list_truth_null[[paste0('k', iter_folds, '_c', iter_clust)]]
        res[[paste0('iter_', iter_pipeline)]]$Classification_Settings = tmp_validation[[paste0('iter_', iter_pipeline)]]$Classification_Settings
      }
    }
  }
  saveRDS(res, paste0(path, '/training_classification_results.rda'))
  return(res)
}


