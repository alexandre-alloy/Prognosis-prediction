#DEPRECATED -- DELETE FILE

#function for placing patients into folds and clusters


#function to compute vst on each individual fold

Separate_vst <- function(GEM, interactome,
                         patient_group_1 = NULL, patient_group_2 = NULL,
                         folds_list, 
                         folds_to_compute = NULL,
                         nb_clusters = NULL, 
                         clustering_algo = c('consensus clustering', 'Mclust', 'none'), 
                         clustering_subset = NULL, 
                         filename = NULL){
  
  GEM = as.matrix(GEM)
  list_training_vst = list() 
  list_validation_vst = list()
  list_c_training = list()
  list_c_validation = list()
  
  #1- COMPUTE VST MATRICES
  for (fold_iter in folds_to_compute){
    training_samples = unlist(folds_list[["list_training"]][fold_iter])
    validation_samples = unlist(folds_list[["list_validation"]][fold_iter])
    
    if (run_vst == TRUE){
      list_training_vst[[fold_iter]] = DESeq2::varianceStabilizingTransformation(GEM[ , training_samples])
      list_validation_vst[[fold_iter]] = DESeq2::varianceStabilizingTransformation(GEM[ , c(training_samples, validation_samples)])  
    }
    if (run_vst == FALSE){
      list_training_vst[[fold_iter]] = GEM[ , training_samples]
      list_validation_vst[[fold_iter]] = GEM[ , c(training_samples, validation_samples)]
    }
    
    
  #2- Compute VIPER matrices
    GEM_training = as.matrix(list_training_vst[[fold_iter]])
    GEM_validation = as.matrix(list_validation_vst[[fold_iter]])
    
    vps_training_validation = viperSignature(cbind(GEM_training, GEM_validation), 
                                             GEM_training,
                                             per = 0,
                                             method = 'mean',
                                             verbose = F)
    
    vp_training_validation = viper(vps_training_validation,
                                   regulon = interactome, 
                                   method = 'none', 
                                   verbose = F)
  
    vps_training = viperSignature(GEM_training, GEM_training, per = 0, method = 'zscore', verbose = F)
    vp_training = viper(eset = vps_training, regulon = interactome, method = 'none', verbose = F)
    vp_validation = vp_training_validation[ , unlist(folds_list[["list_validation"]][fold_iter])]
    
    #3-Clustering
    if (nb_clusters == 1) clustering_algo = 'none'
    if (clustering_algo == 'none'){  #if there is only 1 cluster: no clustering and no validation set assignments to cluster
      c = list()
      c = rep(1, ncol(vp_training) )
      names(c) = colnames(vp_training)
      
      c_validation = list()
      c_validation = rep(1, length(unlist(folds_list[["list_validation"]][fold_iter])) )
      names(c_validation) = colnames(vp_validation)
    }
    
    else{ #if there is a clustering algorithm (number of clusters > 1)
      c = clustering(vp_matrix = vp_training, algorithm = clustering_algo, nb_clusters = nb_clusters, clustering_subset = clustering_subset)  
      
      #5- PREDICT CLUSTER # FOR VALIDATION PATIENTS
      c_validation = assign_validation_patients_to_clusters(vp_training = vp_training, 
                                                            vp_validation = vp_validation, 
                                                            cluster_assignments = c)
    }
    
    list_c_training[[fold_iter]] = c   #Assign clustering results to a list indexed by fold numbers
    
    if (length(c_validation) == 1){ #when running LOOCV, ensure patient names are returned, not the cluster number
      names(c_validation) = unlist(folds[["list_validation"]][fold_iter])
    }
    
    list_c_validation[[fold_iter]] = c_validation  
    
    #Identify the top MRs for classifying the two groups
    actual_nb_clusters = length(unique(c)) #with Mclust sometimes you get fewer clusters than the set parameter
    
    
    
    
    }
  
  
  res = list(folds = folds_list,
             list_training_vst = list_training_vst,
             list_validation_vst = list_validation_vst,
             list_c_training = list_c_training,
             list_c_validation = list_c_validation,
             patient_group_1 = patient_group_1,
             patient_group_2 = patient_group_2)

  if (is.null(filename) == FALSE) saveRDS(object = res, file = filename)
  
  return(res)
      
}


test_vst <- Separate_vst(GEM = GEM_counts, interactome = pruneRegulon(network_vst), patient_group_1 = a$patient_group_1, patient_group_2 = a$patient_group_2, 
                         folds_list = a$folds, folds_to_compute = 3, nb_clusters = 2, clustering_algo = 'Mclust', filename = './test_vst.rda'
                          )
  
  
for (i in c(2, 5, 10, 15, 20, 25)){
  assign(x = paste0('folds_', i), 
         value = Separate_folds_function(patient_group_1 = High_risk, patient_group_2 = Low_risk, CV_type = 'kfold', nb_folds = i)
           )
}



