#test on a data set with significant AUCs (AML)
#Could also implement a scheme to select features by their p-values
#get p-value from NES: 1 - pnorm(q = abs(NES), mean = 0, sd = 1)

MRs_per_fold <- function(Risk_group_classification_object = NULL, 
                         k, 
                         c,
                         pipeline_iter = 1,
                         p_value_threshold = 0.05,
                         save_path = './intermediate_files'){
  
  res = list()
  final_res = list()
  
  if (is.null(Risk_group_classification_object)) Risk_group_classification_object = readRDS(paste0(save_path,'/classification_results.rda'))
  
  for (pipeline_iterations in pipeline_iter){
    for (fold_iter in k){
      for (clust_iter in c){
        #count the number of AUCs in the validation list (if 1, then p-value is calculated with a normally distributed null; if > 1, p-value is calculated by t-test)
        #find the largest sequence of features that yield significant AUCs
        #if p-values are normally distributed, the minimum p-value is equal to 1/nb observations
        #return p-values and return significant MRs
        
        list_AUC = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC']][[paste0('k', fold_iter, '_c', clust_iter)]]
        list_AUC_neg_ctrl = Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_AUC_neg_ctrl']][[paste0('k', fold_iter, '_c', clust_iter)]]
        
        vector_means_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, mean) 
        
        #Check how many AUCs were computed (if only one, stochastic method will be false)
        if ( length(list_AUC[,1]) == 1 ) {
          stochastic_method <- FALSE
        }else{
          stochastic_method <- TRUE  
          }
        
        
        #if true, there are two AUC distributions: one for the negative control and one for the validation. p-value is calculated by a t-test
        #(stochastic_method == FALSE) -- if false, there is one AUC distribution for the negative control and one single point for the negative control, p-value is calculated by computing the AUC's z-score in the negative control's distribution (minimum: 1/nb samples)
        if (stochastic_method == TRUE) {
          nb_features_tested = ncol(list_AUC)
          #invert signs of AUC if systematically under 0.5
          if (sum(apply(list_AUC, 2, mean) < 0.5)/ncol(list_AUC) > 0.7) { #if AUC is systematically of the wrong sign, invert the signs
            list_AUC = -list_AUC
            list_AUC_neg_ctrl = - list_AUC_neg_ctrl
          }
          #compute p-values by t-test
          vector_p_vals = sapply(1:ncol(list_AUC), function(i) t.test(list_AUC[,i], list_AUC_neg_ctrl[,i])$p.value )
          vector_p_vals = which(vector_p_vals <= p_value_threshold)
          
          vector_AUC_above_neg_ctrl = which(apply(list_AUC, 2, mean) > vector_means_AUC_neg_ctrl)

          
        }else { #stochastic_method FALSE
          nb_features_tested = length(list_AUC)
          #invert AUC signs if systematically under 0.5
          if (sum(list_AUC < 0.5)/length(list_AUC) > 0.7) { #if AUC is systematically of the wrong sign, invert the signs
            list_AUC = -list_AUC
            list_AUC_neg_ctrl = - list_AUC_neg_ctrl
          }  
          #compute p-value by modeling AUC_neg_ctrl as a normal distribution
          
          vector_sds_AUC_neg_ctrl = apply(list_AUC_neg_ctrl, 2, sd)
          vector_p_vals = sapply(1:length(list_AUC), function(i) pnorm(q = list_AUC[i], mean = vector_means_AUC_neg_ctrl[i], sd = vector_sds_AUC_neg_ctrl[i])) 
          vector_p_vals[which(vector_p_vals < 1/length(list_AUC_neg_ctrl[,1]))] = 1/length(list_AUC_neg_ctrl) #the smallest p-value must be 1/nb-negative controls in distribution
          vector_p_vals = which(vector_p_vals <= p_value_threshold)
          
          vector_AUC_above_neg_ctrl = which(list_AUC > vector_means_AUC_neg_ctrl)
        }
        
        #Condition for selecting the top features: minimum index value when continuity from 1 is broken in:
        #1- indices for significant p-values
        #2- indices for AUCs above negative control
        print(vector_p_vals)
        print(vector_AUC_above_neg_ctrl)
        Top_features = min( setdiff(0:nb_features_tested, vector_p_vals),
                            setdiff(0:nb_features_tested, vector_AUC_above_neg_ctrl) )
        
        print(Top_features)
        if (Top_features != 0){
          Top_features = 1:Top_features ; print(Top_features)  
        }
        
        
        if (length(Top_features) > 1){ 
          final_feature_selection = names(Risk_group_classification_object[[paste0('iter_', pipeline_iterations)]][['list_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]])[Top_features]  #find the names of the features
          res[[paste0('iter_', pipeline_iterations)]][['significant_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]] = final_feature_selection
          saveRDS(res, paste0(save_path, '/significant_MRs_p', pipeline_iterations, '_k', fold_iter, '_c', clust_iter, '.rda'))
        }
        else{
          res[[paste0('iter_', pipeline_iterations)]][['significant_MRs']][[paste0('k', fold_iter, '_c', clust_iter)]] = list()
          saveRDS(res, paste0(save_path, '/significant_MRs_p', pipeline_iterations, '_k', fold_iter, '_c', clust_iter, '.rda'))
        }
        
        final_res = c(final_res, res)
        res = list()
      }#end clust_iter
    }#end fold_iter
  }#end pipeline_iterations
  return(final_res)
}#end function



#validation of the features found in MRs_per_fold in the heldout data set
#for each fold, if there are significant features, use them to classify the heldout data
#cluster heldout data using the same clustering algorithm
#find the cluster corresponding to 
Internal_validation <- function(){
  
}





#identify the clusters that are closer (this step will be done sooner, in the clustering script, so that the cluster numbering scheme is consistent across folds)
#find the features that come up in the cluster x number of times (above a certain threshold)
consolidate_MRs_per_cluster <- function(){
  
}








