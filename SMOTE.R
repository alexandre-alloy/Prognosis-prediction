#Also implement TOMEK links
#oversampling with SMOTE followed by undersampling with TOMEK


vp_GEM <- viper(GEM, regulon = pruneRegulon(network_coTF_TF), method = 'rank', verbose = T)

#must return two vectors: list validation and list training
SMOTE <- function(input_matrix, 
                  k = 5, #bug when k == 1
                  sample_group_1,
                  sample_group_2){
  require(FNN)
  
  minority_sample_group = unlist(list(sample_group_1, sample_group_2)[which.min(list(length(sample_group_1), length(sample_group_2)))])
  majority_sample_group = unlist(list(sample_group_1, sample_group_2)[which.max(list(length(sample_group_1), length(sample_group_2)))])
  ratio_for_parity = length(majority_sample_group)/length(minority_sample_group)
  target_number_samples = length(majority_sample_group)
  knn = get.knn(data = t(input_matrix[,minority_sample_group]), k = k, algorithm = 'brute')
  
  sampling_of_k_nearest_neighbors = ceiling(target_number_samples/(k*ncol(input_matrix))) #how many k-nearest neighbors should be retained for each sample. to make sure we have enough neighbors for resampling
  list_of_k_nearest_neighbors = sapply(1:length(minority_sample_group), function(i) knn$nn.index[i , ])
  
  
  mat_res = matrix(data = NA, nrow = nrow(input_matrix), ncol = length(majority_sample_group), dimnames = list(row.names(input_matrix), paste0('SMOTE_',1:length(majority_sample_group))))
  for (i in 1:length(majority_sample_group)){
    random_minority_sample = sample(1:length(minority_sample_group), 1, F)
    random_k_neighbor = sample(1:k, 1, F)
    mat_res[,i] = (input_matrix[,minority_sample_group[random_minority_sample]] + input_matrix[,list_of_k_nearest_neighbors[random_k_neighbor,random_minority_sample]])/2
    mat_res[,i] = mat_res[,i] * runif(1, 0, 1)
  }
  
  
  res = cbind(input_matrix[,majority_sample_group], mat_res)

  return(res)
  
}

test_SMOTE <- SMOTE(input_matrix = GEM, k = 5, sample_group_1 = dead, sample_group_2 = alive)
test_SMOTE_TOMEK <- TOMEK(input_matrix = test_SMOTE, 
                          sample_group_1 = grep(pattern = 'TCGA', x = colnames(test_SMOTE), value = T), 
                          sample_group_2 = grep(pattern = 'SMOTE', x = colnames(test_SMOTE), value = T))
pca_SMOTE <- prcomp(t(test_SMOTE))
pca_SMOTE_TOMEK <- prcomp(t(test_SMOTE_TOMEK))
pca_no_SMOTE <- prcomp(t(GEM[,c(dead,alive)]))

plot(pca_no_SMOTE$x[,1], pca_no_SMOTE$x[,2], col = c(rep('blue', length(dead)), rep('red', length(alive))), pch = 21, xlab = 'PC1', ylab = 'PC2', main = 'Raw data')
plot(pca_SMOTE$x[,1], pca_SMOTE$x[,2], col = c(rep('blue', 97), rep('red', 97)), pch = 21, xlab = 'PC1', ylab = 'PC2', main = 'SMOTE')
plot(pca_SMOTE_TOMEK$x[,1], pca_SMOTE_TOMEK$x[,2], col = c(rep('blue', 78), rep('red', 78)), pch = 21, xlab = 'PC1', ylab = 'PC2', main = 'SMOTE + TOMEK')


plot(pca_SMOTE2$x[,1], pca_SMOTE2$x[,2], col = c(rep('blue', 97), rep('red', 97)), pch =21)




sapply(1:(ncol(GEM) * test_SMOTE$sampling_of_k_nearest_neighbors), 
       function(i) 
         (GEM[,i] + GEM[,test_SMOTE$list_of_k_nearest_neighbors[sample(1:test_SMOTE$sampling_of_k_nearest_neighbors, 1, F)]])/2 
)


TOMEK <- function(input_matrix, sample_group_1, sample_group_2){
  require(FNN)
  
  knn = get.knn(data = t(input_matrix[,c(sample_group_1, sample_group_2)]), k = 1, algorithm = 'brute')
  
  samples_1_to_remove = which(sapply(1:length(sample_group_1), function(i) knn$nn.index[i,] > length(sample_group_1)))
  samples_2_to_remove = knn$nn.index[samples_1_to_remove,]
  
  samples_1_to_remove = c(sample_group_1, sample_group_2)[samples_1_to_remove]
  samples_2_to_remove = c(sample_group_1, sample_group_2)[samples_2_to_remove]
  
  TOMEK_matrix = input_matrix[,setdiff(colnames(input_matrix), c(samples_1_to_remove, samples_2_to_remove))]
  
  #downsample majority 
  new_sample_group_1 = colnames(TOMEK_matrix)[which(colnames(TOMEK_matrix) %in% sample_group_1)]
  new_sample_group_2 = colnames(TOMEK_matrix)[which(colnames(TOMEK_matrix) %in% sample_group_2)]
  
  if (length(new_sample_group_1) == length(new_sample_group_2)) TOMEK_matrix = TOMEK_matrix[,c(new_sample_group_1, new_sample_group_2)]
  if (length(new_sample_group_1) > length(new_sample_group_2)) TOMEK_matrix = TOMEK_matrix[,c(sample(x = new_sample_group_1, size = length(new_sample_group_2), replace = F), new_sample_group_2)]
  if (length(new_sample_group_1) < length(new_sample_group_2)) TOMEK_matrix = TOMEK_matrix[,c(new_sample_group_1, sample(x = new_sample_group_2, size = length(new_sample_group_1), replace = F))]
  
  
  return(TOMEK_matrix)
  
}

TOMEK_test_SMOTE$new_sample_group_1
TOMEK_test_SMOTE$new_sample_group_2


TOMEK_test_SMOTE <- TOMEK(input_matrix = test_SMOTE, 
                          sample_group_1 = grep('TCGA', colnames(test_SMOTE), value = T), 
                          sample_group_2 = grep('SMOTE', colnames(test_SMOTE), value = T))



pca_TOMEK_test_SMOTE <- prcomp(t(TOMEK_test_SMOTE))
plot(pca_TOMEK_test_SMOTE$x[,1], pca_TOMEK_test_SMOTE$x[,2], col = c(rep('blue', 59), rep('red', 59)), pch =21)



