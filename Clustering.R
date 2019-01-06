#Implement consensus clustering and dirichlet clustering
#Implement UMAP + Louvain (?)

clustering <- function(input_matrix, algorithm = c('consensus clustering', 'Mclust'), nb_clusters = 3, clustering_subset = NULL){
  input_matrix = as.matrix(input_matrix)
  
  if (is.null(clustering_subset)) {
    clust_input = input_matrix
  }
  else if (clustering_subset == 'variance'){ #keeps the top 2000 features by variance.. should cluster based on dispersion instead
    clust_input = input_matrix[order(rowVars(input_matrix), decreasing = T)[1:2000], ] #if clustering on the most variable genes
  }
  else if (length(clustering_subset) > 1) {
    clust_input = input_matrix[intersect(row.names(input_matrix), clustering_subset), ] #if clustering on c(TF, coTF) or another list of features
  }
    
  
  ###
  
  if (algorithm == 'consensus clustering'){
      #c = ConsensusClusterPlus(d = as.dist(1 - cor(clust_input, method = 'spearman')), maxK = (nb_clusters + 1), reps = 4000, pItem = 0.75, pFeature = 0.75, clusterAlg = 'pam', distance = 'spearman', writeTable = F, verbose = F)
      c = ConsensusClusterPlus(d = clust_input, maxK = (nb_clusters + 1), reps = 4000, pItem = 0.75, pFeature = 0.75, clusterAlg = 'pam', distance = 'pearson', writeTable = F, verbose = F)
    return(c[[nb_clusters]]$consensusClass)
  }
  if (algorithm == 'Mclust'){
    c = Mclust(data = t(clust_input), G = 1:nb_clusters, verbose = F)
    return(c$classification)
  }
}

assign_validation_patients_to_clusters <- function(vp_training, vp_validation, cluster_assignments){
  
  rf = randomForest(x = t(vp_training), y = as.factor(cluster_assignments), keep.forest = T)
  tmp_res = predict(object = rf, newdata = t(vp_validation), type = 'response')
  c_validation = as.numeric(tmp_res)
  names(c_validation) = names(tmp_res)
  return(c_validation)
}























