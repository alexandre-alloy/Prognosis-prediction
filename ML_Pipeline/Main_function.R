library(viper)
library(ConsensusClusterPlus)
library(randomForest)
library(cluster)
library(ROCR)
library(mclust)
library(DESeq2)
library(MASS) #for lda
library(e1071) #for svm
library(FNN) #for get.knn

#Initial processing functions:

#1- Initial_processing_folds:
#Divides patients/samples into training/validation sets in each fold
#If heldout data > 0, will hold out a fraction of the data set for testing the model

#2- Initial_normalizations_clust:
#Runs VST on each fold (if run_vst == FALSE : keep the matrix as in the input (e.g. cpm))
#Runs VIPER on each fold (group_1 vs group_2 to determine top MRs)
#Runs clustering on each fold

#3- Consolidate_initial_processing
#consolidates all the partial vst results into one object / file

#Returns:
# - Fold composition (training patients and validation patients)
# - VST matrices for each fold (training and training + validation) 
# - List of top MRs as computed by VIPER
# - Clustering results in each fold

#Inputs:
#(Initial_processing_folds)
#patient_group_1: list of patients belonging to group 1 (the group with the phenotype e.g. "high risk", "metastatic", "low survival", etc)
#patient_group_2: list of patients belonging to group 2 (the group used as a reference e.g. "low risk", "non metastatic", "healthy", etc)
#CV_type: cross-validation type, k-fold, LOOCV or monte carlo
#nb_folds: number of folds for k-fold CV or number of iterations for monte carlo. ignored if CV_type is "LOOCV"
#validation_fold_size: size of the validation folds relative to the entire data set - value between 0+ and 1- (e.g. 0.3) only used in monte carlo CV, ignored otherwise (in LOOCV validation test is 1 patient; in k-folds it is nb_patients * (1/k))
#filename: if not NULL, will save the results in the designated path

#(Initial_processing_vst_clust)
#GEM: Gene expression matrix (RNA seq matrix of RAW COUNTS) not normalized
#interactome: network resulting from the aracne2regulon function
#patient_group_1: list of patients belonging to group 1 (the group with the phenotype e.g. "high risk", "metastatic", "low survival", etc)
#patient_group_2: list of patients belonging to group 2 (the group used as a reference e.g. "low risk", "non metastatic", "healthy", etc)
#compute_folds: computing vst is very slow. computing vst 5, 10, 20 or n times can be prohibitively slow. with compute_folds you can compute vst on individual folds and then merge the results at the end -- ONLY COMPUTE 1 FOLD AT A TIME AND SAVE EACH FOLD IN A SEPARATE FILE NAMED "experimentname_1.rda" , ... "experimentname_k.rda"
#nb_clusters: number of desired clusters. You should test several different values (1,2,3,4,5), depending on the size of your data set, and pick the model with lowest error
#clustering_algo: consensus clustering, Mclust or none. If nb_clusters == 1 then "none" is selected by default
#clustering subset: cluster on a subset of the features of the viper matrix. You can, for example, subset on the inferred activity of the TFs and coTFs only in a viper matrix containing TFs, coTFs and Sig
#classification_algorithm: random forest or binomial regression (only random forest is currently implemented)
#nb_clusters: number of clusters to divide the data set into (default = 1 - no clustering) Use when you suspect that your data will subdivide into clusters
#clustering_subset: for Mclust, will cluster based on a subset of TRs. for example, you can subset the viper matrix to the TFs and/or coTFs and ignore sig .. example: clustering_subset = c(TF, coTF) 
#mRNA_control: if TRUE, will run the pipeline on the mRNA expression matrix (for comparison of classification accuracy using only mRNA expression vs using VIPER activities)
#filename: if not NULL, will save the results in the designated path

#(Consolidate_initial_processing)
#directory: directory where all the partial files are
#filename_prefix: all filenames should have the same name (prefix) followed by an '_', followed by a number corresponding to a value of k_iter in the interval 1:k. do not include "_" in the prefix
#filename: if not NULL, will save the results in the designated path

#ORDER OF THE FUNCTIONS:
#1- Initial_processing_folds (calls Patient_groupings)
#2- Initial_normalization_clust (calls clustering)
#3- Validation (calls Random_forest_classification) 



#need to fix patient subsampling in equlibrate_classes and equilibrate_classes_heldout_data == TRUE

#implement option to oversample minority class or subsample majority class

Initial_processing_folds <- function(patient_group_1, patient_group_2,
                                    CV_type = c('kfold', 'LOOCV', 'monte_carlo'), 
                                    nb_folds = 5, validation_fold_size = NULL, nb_clusters = 1,
                                    holdout_data = 0.1, #fraction 0 to 1, % of samples to hold out for testing 
                                    compute_null = TRUE,
                                    equilibrate_classes = FALSE,
                                    downsample_majority_class = TRUE,
                                    oversample_minority_class = FALSE,
                                    SMOTE = FALSE, TOMEK = FALSE, 
                                    GEM = NULL,
                                    interactome = NULL,
                                    save_path = './intermediate_files',
                                    savefile = TRUE){

 
  
  
  if (holdout_data > 0 & 1%%holdout_data != 0 ) stop('ERROR: heldout fraction must divide 1 without a remainder (e.g. 0.05, 0.1, 0.2, etc)')
  if (holdout_data >= 1 | holdout_data < 0) stop('ERROR: heldout data must be a number between 0 and 1, incl.')
  
  CV_testing = FALSE #if TRUE: will reiterate the training/validation phase until all heldout data has been tested on
  
  if (SMOTE == TRUE){
    if (is.null(GEM)) stop('ERROR: if SMOTE is TRUE, you must pass a gene expression matrix (GEM)')
    if (sum(!c(patient_group_1, patient_group_2) %in% colnames(GEM))) stop('ERROR: the gene expression matrix (GEM) does not correspond to the patient groups 1 and 2')
    
    res_SMOTE = SMOTE(input_matrix = GEM, k = 5, sample_group_1 = patient_group_1, sample_group_2 = patient_group_2)
    GEM = res_SMOTE
    
    if (patient_group_1[1] %in% colnames(res_SMOTE)) condition_patient_group_1 = TRUE
    if (patient_group_2[1] %in% colnames(res_SMOTE)) condition_patient_group_1 = FALSE
    
    if (condition_patient_group_1){
      patient_group_1 = patient_group_1
      patient_group_2 = grep(pattern = 'SMOTE', x = colnames(test_SMOTE), value = T)
    } else {
      patient_group_2 = patient_group_2
      patient_group_1 = grep(pattern = 'SMOTE', x = colnames(test_SMOTE), value = T)
    }
    
    if (TOMEK == TRUE){
      res_TOMEK = TOMEK(input_matrix = GEM, sample_group_1 = patient_group_1, sample_group_2 = patient_group_2)
      GEM = res_TOMEK
      patient_group_1 = colnames(GEM)[which(colnames(GEM) %in% patient_group_1)]
      patient_group_2 = colnames(GEM)[which(colnames(GEM) %in% patient_group_2)]
    }
    
  }
  
  patients = c(patient_group_1, patient_group_2) #initial list with all the patients, patients will be removed from that list in the while loop
  heldout_patients = NULL  #list of heldout patients within each iteration, it is reset at the end of the iteration
  cumulative_heldout_patients = NULL #cumulative list of heldout patients in previous iterations, it is not reset at the end of the iteration (to prevent the same patients to be heldout more than once)
  
  #check the value in holdout_data (must be 1, 0 or a value in between)
  if (holdout_data > 0) {
    CV_testing = TRUE
    list_heldout_folds = list()
    }
  
  if (holdout_data != 0) nb_iterations_pipeline = 1/holdout_data
  if (holdout_data == 0) nb_iterations_pipeline = 1
  ######
  
  list_folds = list() #the folds with the heldout patients will be stored in this list, the items will be stored in a slot named "iter_1, iter_2, etc"
  
  list_folds$Experimental_Settings$CV_type = CV_type
  list_folds$Experimental_Settings$nb_folds = nb_folds
  list_folds$Experimental_Settings$nb_clusters = nb_clusters
  list_folds$Experimental_Settings$holdout_data = holdout_data
  list_folds$Experimental_Settings$nb_iterations_pipeline = nb_iterations_pipeline
  list_folds$Experimental_Settings$equilibrate_classes = equilibrate_classes
  list_folds$Experimental_Settings$interactome = interactome
    
  
  for (pipeline_iterations in 1:nb_iterations_pipeline){
    
    
    heldout_patients_group_1 = sample(x = setdiff(patient_group_1, cumulative_heldout_patients), size = floor((holdout_data) * length(patient_group_1)), replace = F)
    heldout_patients_group_2 = sample(x = setdiff(patient_group_2, cumulative_heldout_patients), size = floor((holdout_data) * length(patient_group_2)), replace = F)  
    cumulative_heldout_patients = c(cumulative_heldout_patients, heldout_patients_group_1, heldout_patients_group_2)
    patients = c(patient_group_1, patient_group_2)
  
    
    
    if ( ( (CV_type == 'kfold') | (CV_type == 'LOOCV') | (CV_type == 'monte_carlo') ) == FALSE){
      stop('ERROR: CV_type incorrectly set. Must be either kfold, LOOCV or monte_carlo')
    }
    
    if (CV_type == 'kfold') overlapping = FALSE
    if (CV_type == 'monte_carlo') overlapping = TRUE
    if (CV_type == 'LOOCV'){
      nb_folds = length(patients)
      overlapping = FALSE
    }
    
    #calls function Patient_groupings to create the samplings in the k-folds, returns a list of patients in each fold
    folds = Patient_groupings(samples = patients,
                              samples_group_1 = patient_group_1,
                              samples_group_2 = patient_group_2,
                              nb_iterations = nb_folds, 
                              validation_fold_size = validation_fold_size, 
                              overlapping = overlapping,
                              equilibrate_classes = equilibrate_classes, 
                              subsample_majority_class = downsample_majority_class, 
                              oversample_minority_class = oversample_minority_class, 
                              SMOTE = SMOTE, 
                              GEM = GEM)
    
    res = list(folds = folds,
               patient_group_1 = patient_group_1,
               patient_group_2 = patient_group_2)
    
    if (exists('heldout_patients_group_1')){
      res$heldout_patients_group_1 = heldout_patients_group_1
      res$heldout_patients_group_2 = heldout_patients_group_2
    }
    
    list_folds[[paste0('iter_', pipeline_iterations)]] = res
    
  } #end the nb_iterations_pipeline loop
  
  if (savefile == TRUE) {
    if (dir.exists(save_path) == FALSE) dir.create(save_path)
    saveRDS(object = list_folds, file = paste0(save_path, '/folds.rda'))
    }
  
  
  return(list_folds)
}



Initial_normalization_clust <- function(Initial_processing_folds_object, #object name resulting from Initial_processing_folds
                                        GEM = NULL, #RNA seq counts matrix (if normalized, set run_vst to FALSE, if not normalized, set run_vst to TRUE)
                                        interactome = NULL, #aracne2regulon output
                                        compute_folds = 3, #compute individual folds 
                                        nb_clusters = 1, #preset number of clusters before run
                                        cluster_GEM = TRUE, #if TRUE will cluster on the vst-GEM; if FALSE will cluster on the VIPER matrix
                                        run_vst = TRUE,
                                        clustering_algo = c('consensus clustering', 'Mclust', 'UMAP', 'none'), #clustering algorithm #implement UMAP
                                        clustering_subset = NULL, #if not null, cluster viper matrix on subset of features listed here
                                        mRNA_control = FALSE, #run pipeline once with mRNA_control set to TRUE to compare classification accuracy for mRNA expression vs viper activity
                                        pipeline_iter = 1, #list of pipeline iterations to compute (e.g.: c(1,2,3), 1:10, 3, etc), default: 1 (if the pipeline is only run once)
                                        save_path = './intermediate_files', #save result in file name
                                        savefile = TRUE,
                                        loadfile = TRUE
                               ){
  
  require(viper)
  require(ConsensusClusterPlus)
  require(randomForest)
  require(cluster)
  require(ROCR)
  require(mclust)
  require(DESeq2)
  require(MASS)
  require(e1071) 
  require(FNN) 
  
  if (loadfile == TRUE) {
    Initial_processing_folds_object = readRDS(file = './intermediate_files/folds.rda')
    interactome = Initial_processing_folds_object$Experimental_Settings$interactome
    }
  
  
  normalization_res = list()
  
  
  for (pipeline_iterations in pipeline_iter){
    message(paste0('Computing pipeline iteration ', pipeline_iterations, '\n'))
    if (is.null(Initial_processing_folds_object[[paste0('iter_',pipeline_iterations)]]$folds$GEM) == FALSE || loadfile == TRUE) GEM = Initial_processing_folds_object[[paste0('iter_',pipeline_iterations)]]$folds$GEM
    
    GEM = as.matrix(GEM)
    patient_group_1 = Initial_processing_folds_object[[ paste0('iter_',pipeline_iterations) ]]$patient_group_1
    patient_group_2 = Initial_processing_folds_object[[ paste0('iter_',pipeline_iterations) ]]$patient_group_2
    
    #lists containing the cluster assignments
    list_c_training = list()
    list_c_validation = list()
    
    #lists containing the vst-normalized matrices *** change variable names to end in "normalized" instead of "vst"
    list_training_vst = list()
    list_validation_vst = list()
    
    #1- ASSIGN PATIENTS TO FOLDS
    folds = Initial_processing_folds_object[[paste0('iter_',pipeline_iterations)]]$folds
    
    
    res = list() #the list that will contain all the results (patient subsets and normalized matrices)
    #results that are common to all folds (but not all pipeline iterations)
    #will be copied only once in the norm_clust object, when fold_iter == 1
    res_all_folds = list()
    res_all_folds$folds = folds
    res_all_folds$patient_group_1 = patient_group_1
    res_all_folds$patient_group_2 = patient_group_2
    res_all_folds$heldout_patients_group_1 = Initial_processing_folds_object[[paste0('iter_', pipeline_iterations) ]]$heldout_patients_group_1
    res_all_folds$heldout_patients_group_2 = Initial_processing_folds_object[[paste0('iter_', pipeline_iterations) ]]$heldout_patients_group_2
    
    for (fold_iter in compute_folds){
      message(paste0('Computing fold ', fold_iter, '\n'))
      #2- COMPUTE VST MATRICES
      training_samples = unlist(folds[["list_training"]][fold_iter])
      validation_samples = unlist(folds[["list_validation"]][fold_iter])
     
          
      if (run_vst == TRUE){
        vst_res = DESeq2::varianceStabilizingTransformation(GEM[ , c(training_samples, validation_samples)]) #may violate the independance principle between training and validation sets
        #alternative is to compute vst (training + validation) and vst(training) separately
        list_training_vst[[fold_iter]] = vst_res[ , training_samples]
        list_validation_vst[[fold_iter]] = vst_res[ , validation_samples]
      }
      if (run_vst == FALSE){ #if you are supplying a matrix that is already normalized
        list_training_vst[[fold_iter]] = GEM[ , training_samples]
        list_validation_vst[[fold_iter]] = GEM[ , validation_samples]
      }
      
      #3- COMPUTE GLOBAL VIPER MATRICES
      #for each fold, compute the viper matrix for the training samples against the average of the training samples
      GEM_training = as.matrix(list_training_vst[[fold_iter]])
      GEM_validation = as.matrix(list_validation_vst[[fold_iter]])
      #vp_training_validation will be used to compute the AUC for classification with a number of top MRs
      vps_training_validation = viperSignature(cbind(GEM_training, GEM_validation), 
                                               GEM_training, #reference is GEM_training, so that any new patient will have its viper activity computed against the same reference
                                               per = 1000,
                                               method = 'zscore',
                                               verbose = T, bootstrap = FALSE)
      
      if (mRNA_control == FALSE){
        vp_training_validation = viper(vps_training_validation,
                                       regulon = interactome, 
                                       method = 'none', 
                                       verbose = T)  
      }
      else{
        vp_training_validation = vps_training_validation$signature #mRNA control: viper signature matrix
      }
      
      
      #vp training is independent of validation
      vps_training = viperSignature(GEM_training, GEM_training, per = 1000, method = 'zscore', verbose = T, bootstrap = F)
      if (mRNA_control == FALSE) vp_training = viper(eset = vps_training, regulon = interactome, method = 'none', verbose = T)
      if (mRNA_control == TRUE) vp_training = vps_training$signature
      
      #vp validation is not totally independent of training (e.g. new patients will be computed against the average of the training data)
      vp_validation = vp_training_validation[ , unlist(folds[["list_validation"]][fold_iter])]
      
      #4- CLUSTERING GLOBAL TRAINING VIPER MATRIX
      if (nb_clusters == 1) clustering_algo = 'none'
      if (clustering_algo == 'none'){  #if there is only 1 cluster: no clustering and no validation set assignments to cluster
        c = list()
        c = rep(1, ncol(vp_training) )
        names(c) = colnames(vp_training)
        
        c_validation = list()
        c_validation = rep(1, length(unlist(folds[["list_validation"]][fold_iter])) )
        names(c_validation) = colnames(vp_validation)
      }
      
      else{ #if there is a clustering algorithm (number of clusters > 1)
        
        if (cluster_GEM == TRUE){
          c = clustering(input_matrix = list_training_vst[[fold_iter]], algorithm = clustering_algo, nb_clusters = nb_clusters, clustering_subset = clustering_subset)
        }
        else{
          c = clustering(input_matrix = vp_training, algorithm = clustering_algo, nb_clusters = nb_clusters, clustering_subset = clustering_subset)  
        }
        
        
        names(c) = training_samples
        if (length(unique(c)) == 1){
          c_validation = list()
          c_validation = rep(1, length(unlist(folds[["list_validation"]][fold_iter])) )
          names(c_validation) = colnames(vp_validation)
        }
        else{ #if there are more than 1 cluster
          #5- PREDICT CLUSTER # FOR VALIDATION PATIENTS  *******change 'assign_validation_patients_to_clusters' function to use any other kind of classification algorithm (lda, logit etc)
          if (cluster_GEM == TRUE){
            
            c_validation = assign_validation_patients_to_clusters(vp_training = list_training_vst[[fold_iter]], 
                                                                  vp_validation = list_validation_vst[[fold_iter]] , 
                                                                  cluster_assignments = c)  
          }
          else{
            c_validation = assign_validation_patients_to_clusters(vp_training = vp_training, 
                                                                  vp_validation = vp_validation, 
                                                                  cluster_assignments = c)    
          }
          
        }
        
      }
      
      
      list_c_training[[fold_iter]] = c   #Assign clustering results to a list indexed by fold numbers
      
      if (length(c_validation) == 1){ #when running LOOCV, ensure patient names are returned, not the cluster number
        names(c_validation) = unlist(folds[["list_validation"]][fold_iter])
      }
      
      list_c_validation[[fold_iter]] = c_validation  
      
      #Identify the top MRs for classifying the two groups
      actual_nb_clusters = length(unique(c)) #with Mclust sometimes you get fewer clusters than the set parameter
      
      
      
      if (fold_iter == 1){
        res[[paste0('iter_', pipeline_iterations)]] = res_all_folds
        res[[paste0('iter_', pipeline_iterations)]]$list_training_vst = list_training_vst
        res[[paste0('iter_', pipeline_iterations)]]$list_validation_vst = list_validation_vst
        res[[paste0('iter_', pipeline_iterations)]]$list_c_training = list_c_training
        res[[paste0('iter_', pipeline_iterations)]]$list_c_validation = list_c_validation
      }
      else{
        res[[paste0('iter_', pipeline_iterations)]]$list_training_vst = list_training_vst
        res[[paste0('iter_', pipeline_iterations)]]$list_validation_vst = list_validation_vst
        res[[paste0('iter_', pipeline_iterations)]]$list_c_training = list_c_training
        res[[paste0('iter_', pipeline_iterations)]]$list_c_validation = list_c_validation
      }
      
      #when running on the cluster, submit 1 job per fold and pipeline iteration. A file will be saved for all the intermediate steps.
      if (savefile == TRUE) saveRDS(object = res, file = paste0(save_path, '/norm_clust_p', pipeline_iterations, '_k', fold_iter, '.rda'))
      
    }#end fold_iter (to make the folds, compute vst and generate clusters)
    
    #If not running on the cluster, you can run all the folds and pipeline iterations at once and the function will return a merged norm_clust object
    
    res = res_all_folds
    res$list_training_vst = list_training_vst
    res$list_validation_vst = list_validation_vst
    res$list_c_training = list_c_training
    res$list_c_validation = list_c_validation
    
    normalization_res[[paste0('iter_', pipeline_iterations)]] = res

  }#end loop for each iteration

  
  return(normalization_res)
}


