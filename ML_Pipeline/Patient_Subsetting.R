#k-fold (no overlap) or Monte-carlo (with overlap) subsamplings 
#validation_fold_size only used for Monte Carlo subsampling (with overlapping samples) so that an arbitrary validation fold size can be set (in k-folds the fold size is the only parameter)

#oversampling of minority class: naive resampling.. consider implementing SMOTE or ADASYN
#SMOTE and ADASYN might yield interesting results for samples for which only few samples exist (e.g. death without remission)

Patient_groupings <- function(samples, 
                              samples_group_1, 
                              samples_group_2, 
                              nb_iterations, 
                              validation_fold_size, 
                              overlapping = TRUE, 
                              equilibrate_classes = TRUE, 
                              subsample_majority_class = TRUE, 
                              oversample_minority_class = FALSE,
                              SMOTE = FALSE, GEM = NULL){
  
  if ( ((validation_fold_size >= 1) || (validation_fold_size <= 0)) && overlapping == TRUE) stop('ERROR: Validation_fold_size must be a fraction of all samples between 0 and 1')
  
  if (subsample_majority_class == TRUE & oversample_minority_class == TRUE & equilibrate_classes == TRUE) stop('ERROR: subsample_majority_class == TRUE AND oversample_minority_class == TRUE, one of them must be false and the other true if equilibrate_classes == TRUE')
  if (subsample_majority_class == FALSE & oversample_minority_class == FALSE & equilibrate_classes == TRUE) stop('ERROR: subsample_majority_class == FALSE AND oversample_minority_class == FALSE, one of them must be false and the other true if equilibrate_classes == TRUE')
  
  if (SMOTE == TRUE) equilibrate_classes = TRUE
  if ((SMOTE == TRUE) && (is.null(GEM) == TRUE)) stop('ERROR: With SMOTE == TRUE, you need to supply a GEM')
  
  set.seed(as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31))
  
  list_training = list()
  list_validation = list()
  
  list_validation_group_1 = list()
  list_training_group_1 = list()
  list_validation_group_2 = list()
  list_training_group_2 = list()
  
  if (equilibrate_classes == FALSE){
    validation_fold_size = validation_fold_size * length(samples)
    for (i in 1:nb_iterations){
      if (overlapping == TRUE){
        list_validation[[i]] = sample(x = samples, size = validation_fold_size, replace = T)
        list_training[[i]] = sample(x = setdiff(samples, list_validation[[i]]), size = length(samples) - validation_fold_size, replace = T) 
      }
      else{ #overlapping == FALSE
        list_validation[[i]] = sample(x = setdiff(samples, unlist(list_validation)), size = floor(1/nb_iterations*length(samples)), replace = F)
        list_training[[i]] = setdiff(samples, list_validation[[i]])
      }
    }#end nb_iterations
  }
  
  
  if ((equilibrate_classes == TRUE) && (SMOTE == FALSE)){
    
    if (subsample_majority_class == TRUE) nb_samples_to_equilibrate_to = min(length(samples_group_1), length(samples_group_2))
    if (oversample_minority_class == TRUE) nb_samples_to_equilibrate_to = max(length(samples_group_1), length(samples_group_2))
    
    for (i in 1:nb_iterations){
      if (overlapping == TRUE){
        list_validation_group_1[[i]] = sample(x = samples_group_1, size = nb_samples_to_equilibrate_to*validation_fold_size, replace = F)
        list_training_group_1[[i]] = sample(setdiff(samples_group_1, list_validation_group_1[[i]]), size = nb_samples_to_equilibrate_to - length(list_validation_group_1[[i]]), replace = F )
        list_validation_group_2[[i]] = sample(x = samples_group_2, size = nb_samples_to_equilibrate_to*validation_fold_size, replace = F)
        list_training_group_2[[i]] = sample(setdiff(samples_group_2, list_validation_group_2[[i]]), size = nb_samples_to_equilibrate_to - length(list_validation_group_2[[i]]), replace = F )
      }
      
      
      if (overlapping == FALSE){
        if (subsample_majority_class == TRUE){
          if (length(samples_group_1) > length(samples_group_2)){
            subsampled_samples_group_1 = sample(x = samples_group_1, size = nb_samples_to_equilibrate_to, replace = F)
            subsampled_samples_group_2 = samples_group_2
          }
          if (length(samples_group_1) < length(samples_group_2)){
            subsampled_samples_group_2 = sample(x = samples_group_2, size = nb_samples_to_equilibrate_to, replace = F)
            subsampled_samples_group_1 = samples_group_1
          }
          if (length(samples_group_1) == length(samples_group_2)){
            subsampled_samples_group_1 = samples_group_1
            subsampled_samples_group_2 = samples_group_2
          }
        }
        
        if (oversample_minority_class == TRUE){
          if (length(samples_group_1) < length(samples_group_2)){
            subsampled_samples_group_1 = sample(x = samples_group_1, size = nb_samples_to_equilibrate_to, replace = T)
            subsampled_samples_group_2 = samples_group_2
          }
          if (length(samples_group_1) > length(samples_group_2)){
            subsampled_samples_group_2 = sample(x = samples_group_2, size = nb_samples_to_equilibrate_to, replace = T)
            subsampled_samples_group_1 = samples_group_1
          }
          if (length(samples_group_1) == length(samples_group_2)){
            subsampled_samples_group_1 = samples_group_1
            subsampled_samples_group_2 = samples_group_2
          }
        }  
        
        #if overlap is false and equilibrate classes is true, the validation samples in the minority class will be depleted after a few iterations
        #in order to avoid this problem, oversample/subsample the minority/majority classes in the validation groups and the training groups
        #1- find the proportion of unique samples (not repeated)
        pct_unique_samples = min(length(samples_group_1), length(samples_group_2))/length(subsampled_samples_group_1)
        
        sampling_size = floor((1/nb_iterations)*nb_samples_to_equilibrate_to)
        unique_sampling_size_validation = floor(pct_unique_samples * sampling_size)
        repeated_sampling_size_validation = sampling_size - unique_sampling_size_validation
        
        #2-assign unique validation samples in the proportion of the sub/over sampling
        #3-fill up the validation group to the target size by repeating/resampling
        
        if (((length(samples_group_1) < length(samples_group_2)) & oversample_minority_class == TRUE) || ((length(samples_group_1) > length(samples_group_2)) & subsample_majority_class == TRUE)){
          list_validation_group_1[[i]] = sample(x = setdiff(samples_group_1, unlist(list_validation_group_1)), size = unique_sampling_size_validation, replace = F)
          list_validation_group_1[[i]] = c(list_validation_group_1[[i]], sample(x = list_validation_group_1[[i]], size = repeated_sampling_size_validation, replace = T))
          list_validation_group_2[[i]] = sample(x = setdiff(subsampled_samples_group_2, unlist(list_validation_group_2)), size = sampling_size, replace = F)
        }
        if (((length(samples_group_1) > length(samples_group_2)) & oversample_minority_class == TRUE) || ((length(samples_group_1) < length(samples_group_2)) & subsample_majority_class == TRUE)){
          list_validation_group_2[[i]] = sample(x = setdiff(samples_group_2, unlist(list_validation_group_2)), size = unique_sampling_size_validation, replace = F)
          list_validation_group_2[[i]] = c(list_validation_group_2[[i]], sample(x = list_validation_group_2[[i]], size = repeated_sampling_size_validation, replace = T)) 
          list_validation_group_1[[i]] = sample(x = setdiff(subsampled_samples_group_1, unlist(list_validation_group_1)), size = sampling_size, replace = F)
        }
        
        
        
        #4- randomly sample with replacement the training samples (samples that are not in the validation group)
        if (((length(samples_group_1) < length(samples_group_2)) & oversample_minority_class == TRUE) || ((length(samples_group_1) > length(samples_group_2)) & subsample_majority_class == TRUE)  ){
          list_training_group_1[[i]] = sample(subsampled_samples_group_1[ !subsampled_samples_group_1 %in% list_validation_group_1[[i]] ], size = length(subsampled_samples_group_1) - length(list_validation_group_1[[i]]), T)
          list_training_group_2[[i]] = subsampled_samples_group_2[ !subsampled_samples_group_2 %in% list_validation_group_2[[i]] ]
        }
        if (((length(samples_group_1) > length(samples_group_2)) & oversample_minority_class == TRUE) || ((length(samples_group_1) < length(samples_group_2)) & subsample_majority_class == TRUE)){
          list_training_group_2[[i]] = sample(subsampled_samples_group_2[ !subsampled_samples_group_2 %in% list_validation_group_2[[i]] ], size = length(subsampled_samples_group_2) - length(list_validation_group_2[[i]]), T)
          list_training_group_1[[i]] = subsampled_samples_group_1[ !subsampled_samples_group_1 %in% list_validation_group_1[[i]] ]
        }
      }
      list_validation[[i]] = c(list_validation_group_1[[i]], list_validation_group_2[[i]])
      list_training[[i]] = c(list_training_group_1[[i]], list_training_group_2[[i]])
    }
  }
  
  if (equilibrate_classes == TRUE & SMOTE == TRUE & overlapping == FALSE){
    for (i in 1:nb_iterations){
      sampling_size = floor((1/nb_iterations)*length(samples_group_1))
      list_validation_group_1[[i]] = sample(x = setdiff(samples_group_1, unlist(list_validation_group_1)), size = sampling_size, replace = F)
      list_validation_group_2[[i]] = sample(x = setdiff(samples_group_2, unlist(list_validation_group_2)), size = sampling_size, replace = F)
      list_training_group_1[[i]] = setdiff(samples_group_1, list_validation_group_1[[i]])
      list_training_group_2[[i]] = setdiff(samples_group_2, list_validation_group_2[[i]])  
      
      list_validation[[i]] = c(list_validation_group_1[[i]], list_validation_group_2[[i]])
      list_training[[i]] = c(list_training_group_1[[i]], list_training_group_2[[i]])
    }
    
  }
  
  res = list(list_training = list_training,
             list_validation = list_validation,
             GEM = GEM)
  return(res)
  
}




multiple_intersect <- function(...){
  args = list(...)
  res = unlist(args[1])
  for (i in 2:length(args)){
    #print(res); print(args[i])
    res = intersect(res, unlist(args[i]))
  }
  return(unlist(res))
}


