#directory: where all the rda files are for the experiment
#exp_name: the name of the experiment, which is used as prefix before all k and c values
#important to prefix all validation files with "val_" and all null model files with "null_"
#k_range: the entire range of k's (e.g. if doing 5-fold CV, input "1:5")
#c_range: the entire range of c's (e.g. if testing for 3 clusters, input "1:3"; if testing without clustering, input "1")
#MR_range_list_1: a list first elements in each MR_range tested. (e.g. if tested 1:10, 11:20 and 21:30, input "c(1,11,21)")
#MR_range_list_2: a list second elements in each MR_range tested. (e.g. if tested 1:10, 11:20 and 21:30, input "c(10,20,30)")

final_consolidation <- function(directory, 
                                exp_name, 
                                k_range, c_range, 
                                MR_range_list_1, MR_range_list_2, 
                                include_null = FALSE,
                                nb_iterations_val = 100, nb_iterations_null = 100, 
                                filename = NULL){
  k = max(k_range)
  c = max(c_range)
  skip_unlist_step = FALSE #TRUE if file doesn't exist
  res = list()
  
  print(paste0(directory, '/folds_', exp_name, '_k', k, '.rda'))
  tmp_vst_folds_clust_obj = readRDS(paste0(directory, '/folds_', exp_name, '_k', k, '.rda'))
  res$folds = tmp_vst_folds_clust_obj$folds
  res$patient_group_1 = tmp_vst_folds_clust_obj$patient_group_1
  res$patient_group_2 = tmp_vst_folds_clust_obj$patient_group_2
  res$heldout_patients_group_1 = tmp_vst_folds_clust_obj$heldout_patients_group_1
  res$heldout_patients_group_2 = tmp_vst_folds_clust_obj$heldout_patients_group_2
  
  for (k_iter in k_range){
    print(paste0(directory, '/vst_', exp_name, '_k', k, 'c', c, '_k', k_iter, '.rda'))
    
    tmp_vst_folds_clust_obj = readRDS(paste0(directory, '/vst_', exp_name, '_k', k, 'c', c, '_k', k_iter, '.rda'))
    res$list_c_training[[k_iter]] = tmp_vst_folds_clust_obj$list_c_training[[k_iter]]
    res$list_c_validation[[k_iter]] = tmp_vst_folds_clust_obj$list_c_validation[[k_iter]]
    
    for (c_iter in c_range){
      res$list_MRs[[paste0('k', k_iter, '_c', c_iter)]]
      tmp_AUC = list()
      tmp_pred = list()
      tmp_truth = list()
      tmp_AUC_null = list()
      tmp_pred_null = list()
      tmp_truth_null = list()
      tmp_AUC_neg_ctrl = list()
      tmp_pred_neg_ctrl = list()
      tmp_truth_neg_ctrl = list()
      for (MR in 1:length(MR_range_list_1)){
        
        MR_1 = MR_range_list_1[MR]
        MR_2 = MR_range_list_2[MR]
        
        tryCatch({
          tmp_val = readRDS(paste0(directory, '/', 'val_', exp_name, '_k', k, 'c', c, '_k', k_iter, 'c', c_iter, '_MR_', MR_1, '_', MR_2, '.rda'))
          
          
          tmp_AUC = cbind(tmp_AUC, tmp_val$list_AUC[[paste0('k', k_iter, '_c', c_iter)]])
          tmp_pred = cbind(tmp_pred, tmp_val$list_pred[[paste0('k', k_iter, '_c', c_iter)]])
          tmp_truth = cbind(tmp_truth, tmp_val$list_truth[[paste0('k', k_iter, '_c', c_iter)]])
          
          if (include_null == TRUE) {
            tmp_null = readRDS(paste0(directory, '/', 'null_', exp_name, '_k', k, 'c', c, '_k', k_iter, '_MR_', MR_1, '_', MR_2, '.rda'))
            print(paste0(directory, '/', 'null_', exp_name, '_k', k, 'c', c, '_k', k_iter, '_MR_', MR_1, '_', MR_2, '.rda'))
            
            tmp_AUC_null = cbind(tmp_AUC_null, tmp_null$list_AUC_null[[paste0('k', k_iter, '_c', c_iter)]])
            tmp_pred_null = cbind(tmp_pred_null, tmp_null$list_pred_null[[paste0('k', k_iter, '_c', c_iter)]])
            tmp_truth_null = cbind(tmp_truth_null, tmp_null$list_truth_null[[paste0('k', k_iter, '_c', c_iter)]])
            
            tmp_AUC_neg_ctrl = cbind(tmp_AUC_neg_ctrl, tmp_null$list_AUC_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
            tmp_pred_neg_ctrl = cbind(tmp_pred_neg_ctrl, tmp_null$list_pred_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
            tmp_truth_neg_ctrl = cbind(tmp_truth_neg_ctrl, tmp_null$list_truth_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
          }
          
          
          tmp_AUC_neg_ctrl = cbind(tmp_AUC_neg_ctrl, tmp_val$list_AUC_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
          tmp_pred_neg_ctrl = cbind(tmp_pred_neg_ctrl, tmp_val$list_pred_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
          tmp_truth_neg_ctrl = cbind(tmp_truth_neg_ctrl, tmp_val$list_truth_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]])
        }, error = function(e) { #if file doesn't exist
          # tmp_AUC = cbind(tmp_AUC, NA)
          # tmp_pred = cbind(tmp_pred, NA)
          # tmp_truth = cbind(tmp_truth, NA)
          # 
          # tmp_AUC_neg_ctrl = cbind(tmp_AUC_neg_ctrl, NA)
          # tmp_pred_neg_ctrl = cbind(tmp_pred_neg_ctrl, NA)
          # tmp_truth_neg_ctrl = cbind(tmp_truth_neg_ctrl, NA)
          skip_unlist_step <<- TRUE
        })
        

      }#MR_range
      if (skip_unlist_step == FALSE){
        res$list_MRs[[paste0('k', k_iter, '_c', c_iter)]] = tmp_val$list_MRs[[1]]
        res$list_AUC[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_AUC, c(1,2), unlist) #result is a list of lists. unlist to get lists of matrices
        res$list_pred[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_pred, c(1,2), unlist)
        res$list_truth[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_truth, c(1,2), unlist)
        
        if (include_null == TRUE){
          res$list_AUC_null[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_AUC_null, c(1,2), unlist)
          res$list_pred_null[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_pred_null, c(1,2), unlist)
          res$list_truth_null[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_truth_null, c(1,2), unlist)
        }
        
        
        res$list_AUC_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_AUC_neg_ctrl, c(1,2), unlist)
        res$list_pred_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_pred_neg_ctrl, c(1,2), unlist)
        res$list_truth_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = apply(tmp_truth_neg_ctrl, c(1,2), unlist)  
      }
      if (skip_unlist_step == TRUE){
        res$list_MRs[[paste0('k', k_iter, '_c', c_iter)]] = NA
        res$list_AUC[[paste0('k', k_iter, '_c', c_iter)]] = NA #result is a list of lists. unlist to get lists of matrices
        res$list_pred[[paste0('k', k_iter, '_c', c_iter)]] = NA
        res$list_truth[[paste0('k', k_iter, '_c', c_iter)]] = NA
        
        if (include_null == TRUE){
          res$list_AUC_null[[paste0('k', k_iter, '_c', c_iter)]] = NA
          res$list_pred_null[[paste0('k', k_iter, '_c', c_iter)]] = NA
          res$list_truth_null[[paste0('k', k_iter, '_c', c_iter)]] = NA  
        }
        
        
        res$list_AUC_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = NA
        res$list_pred_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = NA
        res$list_truth_neg_ctrl[[paste0('k', k_iter, '_c', c_iter)]] = NA
      }
      
      skip_unlist_step = FALSE
      
    }#c_iter
  }#k_iter
  
  if (is.null(filename) == FALSE) saveRDS(res, filename)
  return(res)
  
}#end function

test_consolidated <- final_consolidation(directory = './results/Consolidate/MC_R_NR_3CV_c4_diag_noTris21_noLowCounts_TF_coTF_T',
                                         exp_name = 'MC_R_NR_3CV_c4_diag_noTris21_noLowCounts_TF_coTF_T_2', 
                                         k_range = 1:3, 
                                         c_range = 1:4,
                                         MR_range_list_1 = c(1), 
                                         MR_range_list_2 = c(50), 
                                         nb_iterations_val = 100,
                                         nb_iterations_null = 100,
                                         filename = './results/Consolidate/result.rda')




