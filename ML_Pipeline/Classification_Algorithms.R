#Random forest
#Binomial regression

#implement parallelization:
#1- find number of cores; 
#2- partition features across the cores randomly (so that one core doesn't end up with the fastest computations (e.g. MRs = 1:10) ); 
#3- compute classification in parallel
#4- merge the result

classification <- function(classification_method = c('random_forest', 'logistic_regression', 'lda','svm', 'xgboost', 'neural_nets'), 
                           vp_matrix, 
                           group_1_training,
                           group_2_training, 
                           group_1_validation, 
                           group_2_validation, 
                           features,
                           return_rf = FALSE,
                           bagging = FALSE,
                           pAUC_range = c(0,1) ){
  
  dashes_in_colnames = length(grep('-', colnames(vp_matrix))) > 0
  vp_matrix = data.frame(vp_matrix)
  
  
  if (dashes_in_colnames == TRUE) colnames(vp_matrix) <- gsub('\\.', '-', colnames(vp_matrix))
  
  if (!classification_method %in% c('random_forest', 'logistic_regression', 'lda','svm', 'xgboost', 'neural_nets')) stop('Error: classification method must be one of the following: "random_forest", "logistic_regression", "lda","svm", "xgboost", "neural_nets"')
  
  if (classification_method == 'random_forest') return(Random_forest_classification(vp_matrix = vp_matrix,
                                                                             group_1_training = group_1_training,
                                                                             group_2_training = group_2_training,
                                                                             group_1_validation = group_1_validation,
                                                                             group_2_validation = group_2_validation,
                                                                             features = features,
                                                                             return_rf = return_rf,
                                                                             pAUC_range = pAUC_range
                                                                             ))
  
  
  if (classification_method == 'logistic_regression') return(Logistic_regression_classification(vp_matrix = vp_matrix,
                                                                                         group_1_training = group_1_training,
                                                                                         group_2_training = group_2_training,
                                                                                         group_1_validation = group_1_validation,
                                                                                         group_2_validation = group_2_validation,
                                                                                         features = features,
                                                                                         pAUC_range = pAUC_range))
  
  if (classification_method == 'lda') return(LDA(vp_matrix = vp_matrix,
                                          group_1_training = group_1_training,
                                          group_2_training = group_2_training,
                                          group_1_validation = group_1_validation,
                                          group_2_validation = group_2_validation,
                                          features = features,
                                          pAUC_range = pAUC_range))
  
  if (classification_method == 'svm') return(SVM(vp_matrix = vp_matrix,
                                          group_1_training = group_1_training,
                                          group_2_training = group_2_training,
                                          group_1_validation = group_1_validation,
                                          group_2_validation = group_2_validation,
                                          features = features,
                                          pAUC_range = pAUC_range))
  
}



Random_forest_classification <- function(vp_matrix, 
                                         group_1_training, 
                                         group_2_training, 
                                         group_1_validation, 
                                         group_2_validation, 
                                         features,
                                         return_rf = FALSE,
                                         pAUC_range = c(0,1)
                                         ){
  labels_training = c(rep(0, length(group_1_training)), rep(1, length(group_2_training)))
  labels_validation = c(rep(0, length(group_1_validation)), rep(1, length(group_2_validation)))

  rf = randomForest(x = t(vp_matrix[features , c(group_1_training, group_2_training)]), y = as.factor(labels_training),
                    keep.forest = T, do.trace = F , ntree = 1000)
  
  predictions = predict(object = rf, newdata = t(vp_matrix[features , c(group_1_validation, group_2_validation) ]), type = 'prob')
  predictions = predictions[c(group_1_validation, group_2_validation) , 2]
  
  if (length(group_1_validation) >= 1 & length(group_2_validation) >= 1){ #Only compute AUC if we have at least one sample for each group (otherwise prediction function throws an error saying that we have a number of classes not equal to 2)
    AUC = Calculate_AUC(predictions = predictions, truth = labels_validation, pAUC_range = pAUC_range )
  }
  else{
    AUC = NULL #in LOOCV, cannot compute AUC for a single data point, instead keep the predictions and truth values and compute the AUC later with all the patients pooled together
  }
  
  if (return_rf == TRUE) return(rf)
  
  rm(rf)
  
  return(list(AUC = AUC,
              truth = as.integer(labels_validation),
              predictions = as.double(predictions) ))
}




Logistic_regression_classification <- function(vp_matrix, 
                                               group_1_training, 
                                               group_2_training, 
                                               group_1_validation, 
                                               group_2_validation, 
                                               features,
                                               pAUC_range = c(0,1)){
  
  
  labels_training = c(rep(0, length(group_1_training)), rep(1, length(group_2_training))) 
  labels_validation = c(rep(0, length(group_1_validation)), rep(1, length(group_2_validation))) 

  glm_res = glm(formula = c(labels_training, labels_validation) ~ t(vp_matrix[ features , c(group_1_training, group_2_training, group_1_validation, group_2_validation) ]) , 
                family = binomial(link='logit'), 
                subset = which(colnames(vp_matrix) %in% c(group_1_training, group_2_training)) )

  
  df_colnames = colnames(vp_matrix) #keep the colnames in memory because when converting to a dataframe, dashes are converted to dots
  predictions_res = predict.glm(object = glm_res, newdata = data.frame(t(vp_matrix), row.names = df_colnames), type = 'response', probability = TRUE )[c(group_1_validation, group_2_validation)]

  if (length(group_1_validation) >= 1 & length(group_2_validation) >= 1){ #Only compute AUC if we have at least one sample for each group (otherwise prediction function throws an error saying that we have a number of classes not equal to 2)
    AUC = Calculate_AUC(predictions = predictions_res, truth = labels_validation , pAUC_range = pAUC_range)
  }
  else{
    AUC = NULL #in LOOCV, cannot compute AUC for a single data point, instead keep the predictions and truth values and compute the AUC later with all the patients pooled together
  }

  return(list(AUC = AUC,
              truth = as.integer(labels_validation),
              predictions = as.double(predictions_res)))
   
}




LDA <- function(vp_matrix,
                group_1_training,
                group_2_training,
                group_1_validation,
                group_2_validation,
                features,
                pAUC_range = c(0,1)){
  
  labels_training = c(rep(0, length(group_1_training)), rep(1, length(group_2_training))) 
  labels_validation = c(rep(0, length(group_1_validation)), rep(1, length(group_2_validation))) 
  
  vp_matrix = vp_matrix[features , c(group_1_training, group_2_training, group_1_validation, group_2_validation)]
  
  lda_res <- lda(x = t(vp_matrix), 
                  grouping = c(labels_training, labels_validation),
                  subset = 1:length(labels_training))
  
  predictions = predict(object = lda_res, newdata = as.data.frame(t(vp_matrix)), type = 'response' )$posterior[c(group_1_validation, group_2_validation),2]
  
  
  if (length(group_1_validation) >= 1 & length(group_2_validation) >= 1){ #Only compute AUC if we have at least one sample for each group (otherwise prediction function throws an error saying that we have a number of classes not equal to 2)
    AUC = Calculate_AUC(predictions = predictions, truth = labels_validation, pAUC_range = pAUC_range )
  }
  else{
    AUC = NULL #in LOOCV, cannot compute AUC for a single data point, instead keep the predictions and truth values and compute the AUC later with all the patients pooled together
  }
  
  return(list(AUC = AUC,
              truth = as.integer(labels_validation),
              predictions = as.double(predictions)))
  
  
}


SVM <- function(vp_matrix,
                group_1_training,
                group_2_training,
                group_1_validation,
                group_2_validation,
                features,
                pAUC_range = c(0,1)){
  
  labels_training = c(rep(0, length(group_1_training)), rep(1, length(group_2_training))) 
  labels_validation = c(rep(0, length(group_1_validation)), rep(1, length(group_2_validation))) 
  
  vp_matrix = vp_matrix[features , c(group_1_training, group_2_training, group_1_validation, group_2_validation)]
  
  svm_res <- svm(x = t(vp_matrix), 
                 y = c(labels_training, labels_validation), 
                 type = 'C-classification', 
                 kernel = 'radial', probability = T,
                 subset = 1:length(labels_training), 
                 scale = T )
  
  predictions = attr(predict(object = svm_res, newdata = t(vp_matrix), type = 'response', probability = TRUE ), 'probabilities')[c(group_1_validation, group_2_validation),2] 

  if (length(group_1_validation) >= 1 & length(group_2_validation) >= 1){ #Only compute AUC if we have at least one sample for each group (otherwise prediction function throws an error saying that we have a number of classes not equal to 2)
    AUC = Calculate_AUC(predictions = predictions, truth = labels_validation, pAUC_range = pAUC_range )
  }
  else{
    AUC = NULL #in LOOCV, cannot compute AUC for a single data point, instead keep the predictions and truth values and compute the AUC later with all the patients pooled together
  }
  
  #predictions = as.double(predictions[labels_validation])
  return(list(AUC = AUC,
              truth = as.integer(labels_validation),
              predictions = as.double(predictions)))

}



#Old function, cannot compute partial AUCs
# Calculate_AUC <- function(predictions, truth, return_roc = FALSE, add = FALSE, col = 'black', lwd = 1){
#   predictions = as.vector(predictions)
#   pred = prediction(predictions , truth)
#   perf_AUC = performance(pred,"auc") 
#   AUC = perf_AUC@y.values[[1]]
#   perf_ROC = performance(pred,"tpr","fpr")
#   if (return_roc == TRUE){
#     if (add == FALSE) ROCR::plot(perf_ROC, col = col, lwd = lwd)
#     if (add == TRUE) ROCR::plot(perf_ROC, add = TRUE, col = col, lwd = lwd)
#   }
#   return(AUC)
# }

#New function, can compute partial AUCs
Calculate_AUC <- function(predictions, truth, pAUC_range = pAUC_range){
  require(pROC)
  AUC = auc(response = truth, 
      predictor = predictions, 
      partial.auc = c(pAUC_range[1], pAUC_range[2]), 
      partial.auc.focus = "specificity", 
      allow.invalid.partial.auc.correct = TRUE, 
      partial.auc.correct = TRUE)
  #write.table(pAUC_range, file = './pAUC_range.txt')
  return(AUC)
}




