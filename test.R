#########
#Random GEM
GEM = matrix(data = rnorm(100000, 0, 1), 
             nrow = 1000, ncol = 100, 
             dimnames = list(sapply(1:1000, function(i) paste0(sample(letters , 8), collapse = '')), 
                             sapply(1:100, function(i) paste0(sample(letters , 8), collapse = ''))) 
             )

patients_1 <- colnames(GEM)[1:50]
patients_2 <- colnames(GEM)[51:100]
GEM = GEM - min(GEM)

network_coTF_TF <- readRDS("~/Desktop/Columbia/CalifanoLab/data/ALL/coTF_TF.rds")
network_coTF_TF <- network_coTF_TF[1:200]
names(network_coTF_TF) <- sample(row.names(GEM), size = 200, replace = F)

for (i in names(network_coTF_TF)){
  network_coTF_TF[[i]]$tfmode <- network_coTF_TF[[i]]$tfmode[1:sample(1:50, 1)]
  network_coTF_TF[[i]]$likelihood <- network_coTF_TF[[i]]$likelihood[1:length(network_coTF_TF[[i]]$tfmode)]
  names(network_coTF_TF[[i]]$tfmode) <- sample(setdiff(row.names(GEM), i), length(network_coTF_TF[[i]]$tfmode), F)
}
#Non-random GEM (LAML)
GEM <- read.delim('../../../data/AML/Adult_AML/cpm_LAML_genenames.txt', header = T, sep = '\t', row.names = 1, as.is = T) ; colnames(GEM) <- gsub(pattern = '\\.', replacement = '-', x = colnames(GEM))
network_coTF_TF <- readRDS('../../../data/AML/Adult_AML/aracne2regulon_bulk_LAML_cpm_TF_coTF.rds')
clinical <- read.delim('/Users/alexandre/Desktop/Columbia/CalifanoLab/data/AML/Adult_AML/clinical.cart.2018-12-19/clinical.tsv', header = T, sep = '\t', as.is = T, row.names = 1)
dead <- clinical$submitter_id[which(clinical$vital_status == 'dead')]
alive <- clinical$submitter_id[which(clinical$vital_status == 'alive')]

####################################################################################################
####################################################################################################
####################################################################################################
#TESTING PIPELINE
folds_res <- Initial_processing_folds(patient_group_1 = dead, 
                                      patient_group_2 = alive,
                                      CV_type = 'kfold', 
                                      nb_folds = 5, 
                                      validation_fold_size = NULL, 
                                      holdout_data = 0.2, #fraction 0 to 1, % of samples to hold out for testing 
                                      equilibrate_classes = TRUE, #bug when classes are equal in size and equilibrate_classes is true
                                      downsample_majority_class = TRUE,
                                      oversample_minority_class = FALSE, 
                                      SMOTE = TRUE, TOMEK = TRUE,
                                      GEM = GEM,
                                      savefile = TRUE)
sapply(1:5, function(j) 
  sapply(1:5, function(i) 
    length( intersect( folds_res[[paste0('iter_', i)]]$heldout_patients_group_1, folds_res[[paste0('iter_', j)]]$heldout_patients_group_1) )) )
 

sapply(1:5, function(j) sapply(1:5, function(i) length(intersect(folds_res$iter_1$folds$list_training[[j]], folds_res$iter_1$folds$list_validation[[i]]))))

sapply(1:5, function(j) sapply(1:5, function(i) length(which(folds_res[[paste0('iter_',j)]]$folds$list_training[[i]] %in% folds_res[[paste0('iter_',j)]]$patient_group_1))) )
sapply(1:5, function(j) sapply(1:5, function(i) length(which(folds_res[[paste0('iter_',j)]]$folds$list_training[[i]] %in% folds_res[[paste0('iter_',j)]]$patient_group_2))) )
sapply(1:5, function(j) sapply(1:5, function(i) length(which(folds_res[[paste0('iter_',j)]]$folds$list_validation[[i]] %in% folds_res[[paste0('iter_',j)]]$patient_group_1))) )
sapply(1:5, function(j) sapply(1:5, function(i) length(which(folds_res[[paste0('iter_',j)]]$folds$list_validation[[i]] %in% folds_res[[paste0('iter_',j)]]$patient_group_2))) )





norm_res <- Initial_normalization_clust(GEM = NULL, 
                         interactome = network_coTF_TF, 
                         Initial_processing_folds_object = readRDS('./intermediate_files/folds.rda'), 
                         compute_folds = 1, 
                         nb_clusters = 1, 
                         cluster_GEM = FALSE, 
                         run_vst = FALSE,
                         clustering_algo = 'Mclust', 
                         clustering_subset = NULL, 
                         mRNA_control = FALSE, 
                         pipeline_iter = 1, 
                         savefile = T)


val_test <- validation(Risk_group_classification_object = NULL, 
                                k = 1, 
                                c = 1, 
                                classification_algorithm = 'random_forest',
                                pipeline_iter = 1,
                                nb_iterations = 100, 
                                MR_range = 1:20, #range of ranked MRs for which to compute AUC
                                interactome = network_coTF_TF, 
                                top_and_bottom = TRUE, 
                                equilibrate_top_bottom = FALSE,
                                mRNA_control = FALSE, 
                                random_negative_control = TRUE,
                                equilibrate_classes = FALSE, #different from folds- equilibrate classes. Here will drop off extra patients from the majority class
                                save_path = './intermediate_files', savefile = TRUE)



MRs_per_fold(Risk_group_classification_object = val_test, k = 1, c = 1, pipeline_iter = 1, p_value_threshold = 1)



vp_test <- viper(GEM, network_coTF_TF, method = 'ttest', verbose = T)

test_merge <- merge_norm_clust_validation_files('./intermediate_files/', T)
a <- readRDS('./intermediate_files/norm_clust_p1_k1.rda')






plot(x = 1:20, y = val_test$iter_1$list_AUC$k1_c1, type = 'l', ylim = c(0,1))
lines(apply(val_test$iter_1$list_AUC_neg_ctrl$k1_c1, 2, mean), col = 'blue')
polygon(x = c(1:20, 20:1), y = c(apply(val_test$iter_1$list_AUC_neg_ctrl$k1_c1, 2, mean) - (2/sqrt(100))*(apply(val_test$iter_1$list_AUC_neg_ctrl$k1_c1, 2, sd)),
                                 rev(apply(val_test$iter_1$list_AUC_neg_ctrl$k1_c1, 2, mean) + (2/sqrt(100))*(apply(val_test$iter_1$list_AUC_neg_ctrl$k1_c1, 2, sd)))), 
        col = rgb(1,1,1,0.3))



test_MR_selection <- MRs_per_fold(Risk_group_classification_object = val_test, k = 1, c = 1, pipeline_iter = 1, p_value_threshold = 0.05)





