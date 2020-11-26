#merge results of a CC_obj (create adjacency matrix and normalized adjacency matrix)
create_adjacency_matrix <- function(CC_obj, patient_labels){
  
  clust_mat <- matrix(data = NA, nrow = length(patient_labels), ncol = length(CC_obj), dimnames = list(patient_labels, 1:length(CC_obj)))
  
  # succesful_runs <- sapply(1:length(CC_obj), function(i) length(CC_obj[[i]]))
  # succesful_runs <- which(succesful_runs == 3)
  
  
  for (i in 1:length(CC_obj)){
    #print(i)
    #clust_mat[match(names(CC_obj[[i]]), row.names(clust_mat)), i] <- CC_obj[[i]]
    clust_mat[match(CC_obj[[i]]$samples, row.names(clust_mat)), i] <- CC_obj[[i]]$clusters
  }
  
  adj_mat <- matrix(data = 0, nrow = length(patient_labels), ncol = length(patient_labels), dimnames = list(patient_labels, patient_labels))
  pb <- txtProgressBar(min = 1, max = length(patient_labels), initial = 1, style = 3)
  
  for (patient in 1:length(patient_labels)){
    setTxtProgressBar(pb, patient)
    for (fold in 1:length(CC_obj)){
      clust_patient <- clust_mat[patient,fold]
      common_patients <- which(clust_mat[,fold] == clust_patient)
      adj_mat[common_patients, patient] <- adj_mat[common_patients, patient] + 1
    }
  }
  
  normalized_adj_mat <- t(t(adj_mat)/sapply(1:length(patient_labels), function(i) adj_mat[i,i]))
  
  return(list(adj_mat = adj_mat,
              normalized_adj_mat = normalized_adj_mat))
}



########################
#2- merge adj matrices 

merge_adj_mat <- function(adj_mat_obj, normalized = TRUE){
  if (normalized){
    index_mat <- 2
  } else {
    index_mat <- 1
  }
  nb_matrices <- length(adj_mat_obj)
  mat_sums <-  matrix(NA, nrow = nrow(adj_mat_obj[[1]][[index_mat]]), ncol = ncol(adj_mat_obj[[1]][[index_mat]]), dimnames = list(row.names(adj_mat_obj[[1]][[index_mat]]), colnames(adj_mat_obj[[1]][[index_mat]])))
  mat_means <- matrix(NA, nrow = nrow(adj_mat_obj[[1]][[index_mat]]), ncol = ncol(adj_mat_obj[[1]][[index_mat]]), dimnames = list(row.names(adj_mat_obj[[1]][[index_mat]]), colnames(adj_mat_obj[[1]][[index_mat]])))
  mat_sds <- matrix(NA, nrow = nrow(adj_mat_obj[[1]][[index_mat]]), ncol = ncol(adj_mat_obj[[1]][[index_mat]]), dimnames = list(row.names(adj_mat_obj[[1]][[index_mat]]), colnames(adj_mat_obj[[1]][[index_mat]])))
  pb <- txtProgressBar(min = 1, max = ncol(mat_means), initial = 1, style = 3)
  
  for (i in 1:ncol(mat_means)){
    setTxtProgressBar(pb, i)
    tmp_mat <- matrix(NA, nrow = nrow(adj_mat_obj[[1]][[index_mat]]), ncol = nb_matrices)
    for (run in 1:nb_matrices){
      tmp_mat[,run] <- adj_mat_obj[[run]][[index_mat]][,i]
    }
    sums_res <- rowSums(tmp_mat, na.rm = T)
    means_res <- rowMeans(tmp_mat, na.rm = T)
    sd_res <- rowSds(tmp_mat, na.rm = T)
    mat_means[,i] <- means_res
    mat_sds[,i] <- sd_res
    mat_sums[,i] <- sums_res
  }
  return(list(mean_adj = mat_means,
              sd_adj = mat_sds,
              sums_adj = mat_sums) )
}

###################

adj_mat_to_graph <- function(adj_mat){
  nb_samples <- dim(adj_mat)[1]
  
  edge_list <- matrix(data = NA, 
                      nrow = length(which(adj_mat !=0 )) - nrow(adj_mat), 
                      ncol = 3)
  counter = 0
  pb <- txtProgressBar(min = 0, max = nb_samples, initial = 0, style = 3)
  for (row in 1:nb_samples){
    setTxtProgressBar(pb, row)
    for (col in 1:nb_samples){
      if (row != col){
        if (adj_mat[row,col] != 0){
          counter <- counter + 1
          edge_list[counter,1] <- row.names(adj_mat)[row]
          edge_list[counter,2] <- colnames(adj_mat)[col]
          edge_list[counter,3] <- adj_mat[row,col]
        }  
      }
    }
  }
  
  graph_res <- igraph::graph_from_data_frame(edge_list)
  return(graph_res)
}

a <- adj_mat_to_graph(final_merge_CC_PHATE_LEIDEN)



calculate_euclidean_distance_from_embedding <- function(embedding, cluster_assignments){
  
  nb_clusters <- length(unique(cluster_assignments))
  
  mat_dist <- matrix(data = 0, nrow = nb_clusters, ncol = nb_clusters, dimnames = list(unique(cluster_assignments), unique(cluster_assignments)))
  
  for (i in 1:nb_clusters){
    ind_i <- which(cluster_assignments == i)
    median_1i <- median(embedding[ind_i,1])
    median_2i <- median(embedding[ind_i,2])
    median_3i <- median(embedding[ind_i,3])
    
    for (j in 1:nb_clusters){
      ind_j <- which(cluster_assignments == j)
      median_1j <- median(embedding[ind_j,1])
      median_2j <- median(embedding[ind_j,2])
      median_3j <- median(embedding[ind_j,3])
      
      dist_res <- (median_1i - median_1j)^2 + (median_2i - median_2j)^2 + (median_3i - median_3j)^2
      dist_res <- sqrt(dist_res)
      mat_dist[i,j] <- dist_res
    }
  }
  return(mat_dist)
}

