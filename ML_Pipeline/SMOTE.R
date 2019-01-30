#Also implement TOMEK links
#oversampling with SMOTE followed by undersampling with TOMEK


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


