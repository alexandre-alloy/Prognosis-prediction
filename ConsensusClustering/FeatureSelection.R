feature_ranking <- function(list_VP, 
                            norm_GEM, 
                            list_networks,
                            patient_labels, 
                            cluster_assignments, 
                            cluster_to_test,
                            dist_clust_ranks, 
                            test_group, ref_group, 
                            prop_training = 0.7, prop_validation = 0.3, prop_test = 0.3,
                            test_samples = NULL,
                            nb_features = 50, 
                            nb_iterations = 10,
                            druggable_targets = NULL){
  
  require(MASS)
  require(parallel)
  require(randomForest)
  require(viper)
  require(DESeq2)
  
  
  list_AUC <- list()
  list_AUC_neg_ctrl <- list()
  list_AUC_neg_ctrl_SD <- list()
  list_AUC_null <- list()
  list_AUC_null_SD <- list()
  
  all_samples <- c(test_group, ref_group)
  nb_samples <- length(all_samples)
  
  if (is.null(test_samples)){
    test_set <- sample(all_samples, round(nb_samples * prop_test), F) #left out data set  
  } else {
    test_set <- test_samples
  }
  
  training_set <- setdiff(all_samples, test_set) #training data, which will be further split into training and validation sets
  
  dist_ranks <- dist_clust_ranks[,cluster_to_test]
  dist_ranks <- names(sort(dist_ranks, decreasing = T))
  nb_other_clusters <- length(dist_ranks) - 1
  param_nb_features <- nb_features
  
  for (i in 1:nb_iterations){
    nb_features <- param_nb_features
    print(i)
    
    write(i, './iter.txt')
    
    random_ref_nb_clusters <- sample(3:nb_other_clusters, 1, F) #at each iteration, a number of reference clusters is selected (sorted by distance from test cluster)
    # print(1:random_ref_nb_clusters)
    # print(dist_ranks[1:random_ref_nb_clusters])
    random_ref_nb_clusters <- dist_ranks[1:random_ref_nb_clusters]
    
    train <- sample(training_set, round(prop_training*length(training_set)), F)
    validation <- setdiff(training_set, train)
    
    ref <- names(cluster_assignments[which(cluster_assignments %in% random_ref_nb_clusters)])
    ref_train <- intersect(ref, train)
    ref_val <- intersect(ref, validation)
    
    test_group_training <- intersect(test_group, train)
    test_group_validation <- intersect(test_group, validation)
    
    # print(paste0('ref: ', ref_train))
    # print(paste0('ref_val: ', ref_val))
    # print(paste0('test_samples_training: ', test_group_training))
    # print(paste0('test_samples_validation: ', test_group_validation))
    
    if (length(ref_train) > length(test_group_training)){
      ref_train <- sample(ref_train, length(test_group_training), F)
    }
    
    if (length(ref_train) < length(test_group_training)){
      train <- sample(test_group_training, length(ref_train), F)
    }
    
    #Feature selection - vipersignature
    nb_networks <- length(list_networks)
    vps <- viperSignature(norm_GEM[,test_group_training], norm_GEM[,ref_train], per = 100, method = 'zscore', verbose = F)
    
    #nb_networks > 1 does not work -- if only one network is present, length will give the number of regulators
    # if (nb_networks > 1){ #multiple networks
    #   res_vp <- mclapply(1:nb_networks, function(i){
    #     res <- viper(eset = vps$signature, regulon = pruneRegulon(list_networks[[i]], 100), method = 'none', verbose = F, cores = 1)
    #     res <- rowSums(res)/sqrt(ncol(res)) #stouffer's integration
    #     return(res)
    #   }, mc.cores = nb_networks)
    #   
    #   features <- names(which(table(unlist(sapply(res_vp, function(i) names(i)))) == length(list_networks)))  
    #   
    #   mat_feats <- matrix(NA, nrow = length(features), ncol = length(list_networks), dimnames = list(features, 1:length(list_networks)))
    #   
    #   
    #   for (iter_net in 1:length(list_networks)){
    #     mat_feats[features, iter_net] <- res_vp[[iter_net]][features]
    #   }
    #   mat_feats <- apply(mat_feats, 2, function(i) rank(-i))
    #   mat_feats <- cbind(mat_feats, rowSums(mat_feats))
    #   mat_feats <- mat_feats[order(mat_feats[,ncol(mat_feats)]),]
    #   ordered_features <- row.names(mat_feats)
      
      
    #} else { #one network
      res_vp <- viper(eset = vps$signature, regulon = pruneRegulon(list_networks, 100), method = 'none', verbose = F, cores = 1)
      res_vp <- rowSums(res_vp)/sqrt(ncol(res_vp)) 
      ordered_features <- names(sort(res_vp, decreasing = T))
    #}
    
    
    
    
    mat_AUC <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    mat_AUC_neg_ctrl <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    mat_AUC_null <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    
    
    
    for (vp_iter in 1:length(list_VP)) {
      
      
      
      vp <- list_VP[[vp_iter]]
      
      ordered_features <- intersect(ordered_features, row.names(vp))
      
      if (!is.null(druggable_targets)){
        ordered_features <- intersect(druggable_targets, ordered_features)
      }
      
      if (nb_features > length(ordered_features)){
        nb_features <- length(ordered_features)
      }
      
      if (nb_features > 1){
      for (iter_feat in 2:nb_features){
        #write(iter_feat, './iter_feat.txt')
        #print(ordered_features[1:iter_feat])
        
        

        class_res <- classification(classification_method = 'random_forest', 
                                      vp_matrix = vp, 
                                      group_1_training = test_group_training, 
                                      group_2_training = ref_train, 
                                      group_1_validation = test_group_validation,
                                      group_2_validation = ref_val, 
                                      features = ordered_features[1:iter_feat])  
        
        
        class_res_neg_ctrl <- lapply(1:10, function(iter_neg_ctrl){
          random_feat <- sample(x = ordered_features, size = iter_feat, replace = F)
          print(random_feat)
          res <- classification(classification_method = 'random_forest', 
                         vp_matrix = vp, 
                         group_1_training = test_group_training, 
                         group_2_training = ref_train, 
                         group_1_validation = test_group_validation,
                         group_2_validation = ref_val, 
                         features = random_feat)
          return(res$AUC)})
        
        
        class_res_null <- lapply(1:10, function(iter_null){
          
          null_train_g1 <- sample(c(test_group_training, ref_train), size = length(test_group_training), replace = F)
          null_train_g2 <- setdiff(c(test_group_training, ref_train), null_train_g1)
          
          null_validation_g1 <- sample(c(test_group_validation, ref_val), size = length(test_group_validation), replace = F)
          null_validation_g2 <- setdiff(c(test_group_validation, ref_val), null_validation_g1)
          
          
          res <- classification(classification_method = 'random_forest', 
                         vp_matrix = vp, 
                         group_1_training = null_train_g1, 
                         group_2_training = null_train_g2, 
                         group_1_validation = null_validation_g1,
                         group_2_validation = null_validation_g2, 
                         features = ordered_features[1:iter_feat])
          return(res$AUC)
        })
        
        mat_AUC[iter_feat,vp_iter] <- class_res$AUC
        mat_AUC_neg_ctrl[iter_feat,vp_iter] <- mean(as.numeric(class_res_neg_ctrl), na.rm = T)
        mat_AUC_null[iter_feat,vp_iter] <- mean(as.numeric(class_res_null), na.rm = T)
        
      }
      } else {
      mat_AUC <- NA
      mat_AUC_neg_ctrl <- NA
      mat_AUC_null <- NA
    }
    }
    
    
    list_AUC[[i]] <- mat_AUC
    list_AUC_neg_ctrl[[i]] <- mat_AUC_neg_ctrl
    list_AUC_null[[i]] <- mat_AUC_null
    
    
  }
  
  return(list(list_AUC = list_AUC,
              list_AUC_neg_ctrl = list_AUC_neg_ctrl,
              list_AUC_null = list_AUC_null,
              test_set = test_set))
  
}



error_interval <- function(x, y, sd, 
                           col_interval = 'grey', col_line = 'black', 
                           nb_sd = 1, SE = FALSE, 
                           new_plot = TRUE, 
                           xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
                           title = '', xlab = 'nb_Features', ylab = 'AUC'){
  if (new_plot){
    plot(x = x, y = y, xlim = xlim, ylim = ylim, main = title, xlab = xlab, ylab = ylab, type = 'l', lwd = 3, col = col_line)
  } else {
    lines(x = x, y = y, type = 'l', lwd = 3, col = col_line)
  }
  polygon(x = c(x, rev(x)), y = c((y + (nb_sd * sd)), (rev(y) - rev(nb_sd * sd))  ), border = NA, col = col_interval)
  lines(x = x, y = y, type = 'l', lwd = 3, col = col_line)
}


merge_MRs_folds <- function(AUC, neg_ctrl, null, nb_features, pval_threshold = 0.05, list_VP, patients, test_set){
  
  patients <- setdiff(patients, test_set)
  
  pval <- sapply(2:nb_features, function(i){
    1-pnorm(mean(AUC[i,]), mean = mean(neg_ctrl[i,]), sd = sd(neg_ctrl[i,]))
  })
  stopping_point <- min(which(pval > pval_threshold))
  #print(stopping_point)
  if (is.infinite(stopping_point)){
    stopping_point <- nb_features
  }
  MRs_to_include <- row.names(AUC)[1:stopping_point]
  
  nb_VP <- length(list_VP)
  integrated_NES <- lapply(MRs_to_include, function(MR){
    lapply(list_VP, function(vp) {
      
      if (MR %in% row.names(vp)){
        mean(vp[MR,patients])    
      }
      
    })})
  
  #return(integrated_NES)
  
  
  list_MRs <- unlist(lapply(1:length(MRs_to_include), function(i) mean(unlist(integrated_NES[[i]]))))
  names(list_MRs) <- MRs_to_include
  
  return(list_MRs)
  
}



get_list_MRs <- function(res_classification, iterations = 1:10, list_VP, list_patients, nb_features = 40){
  
  
    a <- lapply(iterations, function(j) {
      
      res <- merge_MRs_folds(AUC = res_classification$list_AUC[[j]], 
                             neg_ctrl = res_classification$list_AUC_neg_ctrl[[j]], 
                             null = res_classification$list_AUC_null[[j]],
                             nb_features = nb_features, 
                             pval_threshold = 0.05, 
                             list_VP = list_VP, 
                             patients = list_patients, 
                             test_set =  res_classification$test_set)
    })
    
    unique(names(unlist(a)))
    mat_MRs <- matrix(data = 0, nrow = length(unique(names(unlist(a)))), ncol = 2, dimnames = list(unique(names(unlist(a))), c('Votes', 'NES')))
    for (i in a){
      mat_MRs[names(i),1] <- mat_MRs[names(i),1] + 1
      mat_MRs[names(i),2] <- mat_MRs[names(i),2] + i
    }

    mat_MRs <- data.frame(mat_MRs)
    mat_MRs <- mat_MRs[with(mat_MRs, order(-Votes, -NES)),]
  
  return(mat_MRs)
  
}




final_test <- function(res_classification, iteration = 3, cluster = 1, nb_iterations = 10, list_MRs, clustering, list_VP, nb_features = 20){
  
  
  list_AUC <- list()
  list_AUC_neg_ctrl <- list()
  list_AUC_null <- list()
  
  
  patients <- names(clustering)
  group1 <- names(which(clustering == cluster))
  ref <- names(which(clustering != cluster))
  
  test_set <- res_classification$test_set
  training_set <- setdiff(patients, test_set)
  
  
  for (i in 1:nb_iterations){
    
    print(i)
    write(paste0('Cluster: ', cluster, ' ; iteration: ', i), './iter.txt', append = T)
  
  group1_training <- intersect(group1, training_set)
  ref_training <- intersect(ref, training_set)
  group1_test <- intersect(group1, test_set)
  ref_test <- intersect(ref, test_set)
  
  training_group <- c(group1_training, ref_training)
  test_group <- c(group1_test, ref_test)
  
  ordered_features <- row.names(list_MRs)
  
  nb_features <- min(length(ordered_features), nb_features)
  
  # unordered_features <- table(unlist(lapply(list_VP, function(i) row.names(i))))
  # unordered_features <- which(unordered_features == length(list_VP))
  # unordered_features <- names(unordered_features)
  
  
  mat_AUC <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
  mat_AUC_neg_ctrl <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
  mat_AUC_null <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
  
  
  if (length(ref_training) > length(group1_training)){
    ref_training <- sample(ref_training, length(group1_training), F)
  }
  
  if (length(ref_training) < length(group1_training)){
    train <- sample(group1_training, length(ref_training), F)
  }
  
  
  
  for (vp_iter in 1:length(list_VP)) {
    
    vp <- list_VP[[vp_iter]]
    ordered_feats_intersect <- intersect(ordered_features, row.names(vp))
    #unordered_feats_intersect <- intersect(unordered_features, row.names(vp))
    unordered_feats_intersect <- row.names(vp)
    
    
    for (iter_feat in 2:nb_features){
  
      
      class_res <- classification(classification_method = 'lda', 
                                  vp_matrix = vp, 
                                  group_1_training = group1_training, 
                                  group_2_training = ref_training, 
                                  group_1_validation = group1_test,
                                  group_2_validation = ref_test, 
                                  features = ordered_feats_intersect[1:iter_feat]
      )
      
      class_res_neg_ctrl <- lapply(1:10, function(iter_neg_ctrl){
        res <- classification(classification_method = 'lda', 
                              vp_matrix = vp, 
                              group_1_training = group1_training, 
                              group_2_training = ref_training, 
                              group_1_validation = group1_test,
                              group_2_validation = ref_test, 
                              features = sample(x = unordered_feats_intersect, size = iter_feat, replace = F))
        return(res$AUC)
      })
      
      class_res_null <- lapply(1:10, function(iter_null){
        
        null_train_g1 <- sample(c(group1_training, ref_training), size = length(group1_training), replace = F)
        null_train_g2 <- setdiff(c(group1_training, ref_training), null_train_g1)
        
        null_test_g1 <- sample(c(group1_test, ref_test), size = length(group1_test), replace = F)
        null_test_g2 <- setdiff(c(group1_test, ref_test), null_test_g1)
        
        res <- classification(classification_method = 'lda', 
                              vp_matrix = vp, 
                              group_1_training = null_train_g1, 
                              group_2_training = null_train_g2, 
                              group_1_validation = null_test_g1,
                              group_2_validation = null_test_g2, 
                              features = ordered_feats_intersect[1:iter_feat])
        return(res$AUC)
      })
      
      mat_AUC[iter_feat,vp_iter] <- class_res$AUC
      mat_AUC_neg_ctrl[iter_feat,vp_iter] <- mean(as.numeric(class_res_neg_ctrl), na.rm = T)
      mat_AUC_null[iter_feat,vp_iter] <- mean(as.numeric(class_res_null), na.rm = T)
      
          
    }
  }
  
  list_AUC[[i]] <- mat_AUC
  list_AUC_neg_ctrl[[i]] <- mat_AUC_neg_ctrl
  list_AUC_null[[i]] <- mat_AUC_null
  }
  
  return(list(list_AUC = list_AUC,
              list_AUC_neg_ctrl = list_AUC_neg_ctrl,
              list_AUC_null = list_AUC_null,
              test_set = test_set))
  
}

#13 minutes (1 iteration)
a <- final_test(res_classification = res_classification, iteration = 3, cluster = 1, nb_iterations = 1, list_MRs = MRs_c1, clustering = final_clust$clustering[patient_passed_threshold], list_VP = list_VP)
final_time <- Sys.time()



rowMeans(a$list_AUC[[1]])
rowMeans(a$list_AUC_neg_ctrl[[1]])
rowMeans(a$list_AUC_null[[1]])

list_MRs <- list(MRs_c1, MRs_c2, MRs_c3, MRs_c4, MRs_c5, MRs_c6, MRs_c7, MRs_c8, MRs_c9, MRs_c10, MRs_c11)

#18:00 -- 00:38 ~6hr40
init_time <- Sys.time()
final_test_all_clust <- mclapply(1:11, function(clust) {
  final_test(res_classification = res_classification, iteration = 3, cluster = clust, nb_iterations = 1, list_MRs = list_MRs[[clust]], clustering = final_clust$clustering[patient_passed_threshold], list_VP = list_VP)
}, mc.cores = 11)
final_time <- Sys.time()
saveRDS(final_test_all_clust, './Results/final_test_all_clust_1_iter.rds') #40 minutes -- 1 iteration

################################################



feature_ranking_clinical <- function(list_VP, 
                            norm_GEM, 
                            list_networks,
                            patient_labels, 
                            test_group, ref_group, 
                            prop_training = 0.7, prop_validation = 0.3, prop_test = 0.3,
                            test_samples = NULL,
                            nb_features = 50, 
                            nb_iterations = 10){
  
  require(MASS)
  require(parallel)
  require(randomForest)
  require(viper)
  require(DESeq2)
  
  
  list_AUC <- list()
  list_AUC_neg_ctrl <- list()
  list_AUC_null <- list()
  
  all_samples <- c(test_group, ref_group)
  nb_samples <- length(all_samples)
  
  if (is.null(test_samples)){
    test_set <- sample(all_samples, round(nb_samples * prop_test), F) #left out data set  
  } else {
    test_set <- test_samples
  }
  
  training_set <- setdiff(all_samples, test_set) #training data, which will be further split into training and validation sets
  
  for (i in 1:nb_iterations){
    #write(i, './iteration_nb.txt')
    
    print(i)
    
    train <- sample(training_set, round(prop_training*length(training_set)), F)
    validation <- setdiff(training_set, train)
    
    ref <- ref_group
    ref_train <- intersect(ref, train)
    ref_val <- intersect(ref, validation)
    
    test_group_training <- intersect(test_group, train)
    test_group_validation <- intersect(test_group, validation)
    
    
    if (length(ref_train) > length(test_group_training)){
      ref_train <- sample(ref_train, length(test_group_training), F)
    }
    
    if (length(ref_train) < length(test_group_training)){
      test_group_training <- sample(test_group_training, length(ref_train), F)
    }
    
    #Feature selection - vipersignature
    nb_networks <- length(list_networks)
    vps <- viperSignature(norm_GEM[,test_group_training], norm_GEM[,ref_train], per = 100, method = 'zscore', verbose = F)
    
    if(nb_networks>1){
      res_vp <- lapply(1:nb_networks, function(i){
        res <- viper(eset = vps$signature, regulon = pruneRegulon(list_networks[[i]], 100), method = 'none', verbose = F, cores = 1)
        res <- rowSums(res)/sqrt(ncol(res))
        return(res)
      })
      
      features <- names(which(table(unlist(sapply(res_vp, function(i) names(i)))) == length(list_networks)))
      
      mat_feats <- matrix(NA, nrow = length(features), ncol = length(list_networks), dimnames = list(features, 1:length(list_networks)))
      
      
      for (iter_net in 1:length(list_networks)){
        mat_feats[features, iter_net] <- res_vp[[iter_net]][features]
      }
      mat_feats <- apply(mat_feats, 2, function(i) rank(-i))
      mat_feats <- cbind(mat_feats, rowSums(mat_feats))
      mat_feats <- mat_feats[order(mat_feats[,ncol(mat_feats)]),]
      ordered_features <- row.names(mat_feats)
    } else{
      res_vp <- viper(eset = vps$signature, regulon = pruneRegulon(list_networks[[i]], 100), method = 'none', verbose = F, cores = 1)
      res_vp <- rowSums(res_vp)/sqrt(ncol(res_vp))
      ordered_features <- names(sort(res_vp, decreasing = T))
    }
    
    
    mat_AUC <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    mat_AUC_neg_ctrl <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    mat_AUC_null <- matrix(data = NA, nrow = nb_features, ncol = length(list_VP), dimnames = list(ordered_features[1:nb_features], names(list_VP)))
    
    
    for (vp_iter in 1:length(list_VP)) {
      # write(vp_iter, './vp_iter.txt')
      
      vp <- list_VP[[vp_iter]]
      
      ordered_features_vp <- intersect(ordered_features, row.names(vp))
      
      for (iter_feat in 2:nb_features){
        # write(iter_feat, './iter_feat.txt')
        
        class_res <- classification(classification_method = 'lda', 
                                    vp_matrix = vp, 
                                    group_1_training = test_group_training, 
                                    group_2_training = ref_train, 
                                    group_1_validation = test_group_validation,
                                    group_2_validation = ref_val, 
                                    features = ordered_features_vp[1:iter_feat]
        )
        
        class_res_neg_ctrl <- lapply(1:10, function(iter_neg_ctrl){
          res <- classification(classification_method = 'lda', 
                                vp_matrix = vp, 
                                group_1_training = test_group_training, 
                                group_2_training = ref_train, 
                                group_1_validation = test_group_validation,
                                group_2_validation = ref_val, 
                                features = sample(x = row.names(vp), size = iter_feat, replace = F))
          return(res$AUC)
        })
        
        class_res_null <- lapply(1:10, function(iter_null){
          
          null_train_g1 <- sample(c(test_group_training, ref_train), size = length(test_group_training), replace = F)
          null_train_g2 <- setdiff(c(test_group_training, ref_train), null_train_g1)
          
          null_validation_g1 <- sample(c(test_group_validation, ref_val), size = length(test_group_validation), replace = F)
          null_validation_g2 <- setdiff(c(test_group_validation, ref_val), null_validation_g1)
          
          res <- classification(classification_method = 'lda', 
                                vp_matrix = vp, 
                                group_1_training = null_train_g1, 
                                group_2_training = null_train_g2, 
                                group_1_validation = null_validation_g1,
                                group_2_validation = null_validation_g2, 
                                features = ordered_features_vp[1:iter_feat])
          return(res$AUC)
        })
        
        mat_AUC[iter_feat,vp_iter] <- class_res$AUC
        mat_AUC_neg_ctrl[iter_feat,vp_iter] <- mean(as.numeric(class_res_neg_ctrl), na.rm = T)
        mat_AUC_null[iter_feat,vp_iter] <- mean(as.numeric(class_res_null), na.rm = T)
        
      }
      
    }
    
    
    list_AUC[[i]] <- mat_AUC
    list_AUC_neg_ctrl[[i]] <- mat_AUC_neg_ctrl
    list_AUC_null[[i]] <- mat_AUC_null
    
    
  }
  
  return(list(list_AUC = list_AUC,
              list_AUC_neg_ctrl = list_AUC_neg_ctrl,
              list_AUC_null = list_AUC_null,
              test_set = test_set))
  
}





final_feature_selection <- function(Top_MRs_list_per_cluster, vp, cluster_assignments, cycles = 10, verbose = T){
  list_features <- list()
  no_repeat_improvement <- NULL
  for (i in 1:cycles){
    if (verbose) {
      print(i)
    }
    features <- unlist(lapply(setdiff(unique(cluster_assignments), no_repeat_improvement), function(Top_MRs) {
      if (!is.na(row.names(Top_MRs_list_per_cluster[[Top_MRs]])[i])){
        return(row.names(Top_MRs_list_per_cluster[[Top_MRs]])[i]) 
      }
    } ))
    feat_class <- c(unlist(list_features), features)
    
    rf_res <- randomForest(x = t(vp[feat_class,names(cluster_assignments)]),
                           y = as.factor(cluster_assignments))
    
    class_error_new <- rf_res$confusion[,'class.error']
    
    if (i == 1){
      list_features[[i]] <- features
      class_error_old <- class_error_new
    } else {
      feats_improvement <- which(class_error_new < class_error_old)
      
      no_repeat_improvement_new <- setdiff(unique(cluster_assignments), feats_improvement)
      no_repeat_improvement <- c(no_repeat_improvement, no_repeat_improvement_new)
      no_repeat_improvement <- unique(no_repeat_improvement)
      feats_improvement <- setdiff(feats_improvement, no_repeat_improvement)
      
      
      if (length(feats_improvement) > 0){
        feats_improvement <- unlist(lapply(feats_improvement, function(j) row.names(Top_MRs_list_per_cluster[[j]])[i]))
        list_features[[i]] <- feats_improvement
      } else {
        list_features[[i]] <- NULL
      }
      class_error_old <- class_error_new
    }
    
  }
  return(list_features)
  
}









LogTTest <- function(x, y) {
  test.res <- t.test(x, y, alternative = 'two.sided')
  log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
  return(log.p)
}





