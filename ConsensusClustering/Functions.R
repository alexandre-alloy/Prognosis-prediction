reticulate::use_python(python = '~/.local/share/r-miniconda/envs/r-reticulate/bin/python', required = T)
library(Seurat)
library(DESeq2)
library(tradeSeq)
library(plotfunctions)
library(mclust)
library(SingleCellExperiment)
library(phateR)

consensus_clustering <- function(dataset, 
                                 size_of_fold_prop = 0.8,
                                 with_replacement = FALSE,
                                 vary_top_features_range = c(0.1,0.3), 
                                 ndim = 2, #can also be two values (c(10,15)) will select a random integer in the range at each iteration
                                 nfolds = 1000, 
                                 dimensionality_red_algo = c('umap', 'pca', 'tsne', 'PHATE'),
                                 clustering_algo = c('mclust', 'pam', 'hc', 'leiden', 'spectral'),
                                 nb_clusters = 10:12,
                                 save_RDS_path = './RDS_obj/cc.rds',
                                 seeds = 1:1000,
                                 cores = cores){
  
  require(tsne)
  require(dbscan)
  require(cluster)
  require(leiden)
  require(uwot)
  require(mclust)
  require(kknn) # for spectral clustering
  require(phateR)
  require(parallel)
  
  if (length(seeds) != nfolds) stop('Number of seeds must equal number of folds.')
  size_of_fold = round(ncol(dataset) * size_of_fold_prop)
  
  clustering <- mclapply(1:nfolds, function(i){
    if (i%%25 == 0){
      write(i, './progress.txt')
    }
    seed = seeds[i]
    print(i)
    message(paste0('iteration ', i, ' of ', nfolds))
    set.seed(seed = seed)
    sampling <- sample(colnames(dataset), size_of_fold, with_replacement)
    set.seed(seed = seed)
    prop_top_genes <- runif(n = 1, min = vary_top_features_range[1], max = vary_top_features_range[2])
    prop_top_genes <- round(prop_top_genes * nrow(dataset))
    variable_genes <- order(rowVars(dataset), decreasing = T)[1:prop_top_genes]
    
    emb <- dim_reduction(dataset = t(dataset[variable_genes,sampling]) , algo = dimensionality_red_algo, k = ndim, seed = seed)
    
    clust_res <- clustering_data(dataset = emb, algo = clustering_algo, k = nb_clusters, seed = seed)
    
    return(list(samples = sampling,
                embedding = emb,
                clusters = clust_res))
  }, mc.cores = cores)
  

  return(clustering)
}

random_k <- function(k, seed){
  set.seed(seed = seed)
  if (length(k) > 1){
    k = round(runif(1, k[1], k[length(k)]))
  }
  return(k)
}


dim_reduction <- function(dataset, algo = '', k, seed){

  if (algo == 'umap'){
    
    k <- random_k(k, seed)
    set.seed(seed = seed)
    init <- sample(c('pca', 'spectral', 'random'), 1, F)
    set.seed(seed = seed)
    emb <- umap(X = dataset, 
                n_neighbors = round(runif(1, 40, 50)), 
                n_components = k, 
                init = init, 
                spread = runif(1, 8, 13), 
                min_dist = runif(1, 0.8,1.2), 
                verbose = T)  
    
    
  } else if (algo == 'pca'){
    
    emb <- prcomp(x = dataset)
    k <- random_k(k, seed)
    emb <- emb$x[,1:k]
    
    
  } else if (algo == 'tsne'){
    emb <- prcomp(x = dataset)
    k_pca <- random_k(c(10,50), seed)
    k_tsne <- random_k(k, seed)
    set.seed(seed = seed)
    emb <- tsne(X = emb$x[,1:k_pca], k = k_tsne, perplexity = round(runif(1,22,38)))
  }
  else if (algo == 'PHATE'){
    
    k_PHATE <- random_k(k, seed)
    
    set.seed(seed = seed)
    decay <- sample(c(1,2), 1, F, c(0.1, 0.9))
    if (decay == 1){
      decay <- NULL
    } else {
      set.seed(seed = seed)
      decay <- runif(n = 1, min = 30, max = 80)
    }
    
    emb <- phate(data = dataset, ndim = 2, knn = k_PHATE, decay = decay, gamma = 1, verbose = F, n.jobs = 1, seed = seed)
    emb <- emb$embedding
  }
###
  return(emb)
}


clustering_data <- function(dataset, algo, k, seed){
  
  ################################
  #Mclust
  if (algo == 'mclust'){
    
    k <- random_k(k, seed)
    
    set.seed(seed = seed)
    clust_res <- tryCatch({
      Mclust(dataset, G = 1:k, verbose = T) 
    }, error = function(e){
      write(paste0('Error, seed = ', seed), './errors_mclust.txt')
      set.seed(seed = seed*2)
      return(Mclust(dataset, G = 1:k, verbose = T) )
    })
    clust_res <- clust_res$classification
  }
  
  #############################
  #PAM
  if (algo == 'pam'){
    
    k <- random_k(k, seed)
    
    metric <- sample(c('euclidean', 'manhattan'), 1, F, c(0.7, 0.3))
    clust_res <- pam(x = dataset, k = k, metric = metric)
    clust_res <- clust_res$clustering
    
  }
  
  ###########################
  #leiden
  if (algo == 'leiden'){
    
    k <- random_k(k, seed)
    
    knn_res <- kNN(x = dataset, k = k)
    snn_res <- sNN(x = knn_res, k = k)
    
    edge_list <- matrix(data = NA, 
                        nrow = (dim(snn_res$id)[1] * dim(snn_res$id)[2]) - length(which(is.na((snn_res$id)))), 
                        ncol = 2)
    counter = 0
    for (row in 1:dim(snn_res$id)[1]){
      for (col in 1:dim(snn_res$id)[2]){
        if (!is.na(snn_res$id[row,col])){
          counter <- counter + 1
          edge_list[counter,1] <- row
          edge_list[counter,2] <- snn_res$id[row,col]
        }
      }
    }
    
    graph_res <- igraph::graph_from_data_frame(edge_list)
    
    clust_res <- leiden(graph_res, 
                        partition_type = 'ModularityVertexPartition', 
                        resolution_parameter = runif(1,0.8,1.2), 
                        seed = seed)
  }
  
  ###########################
  #Hierarchical clustering
  
  if (algo == 'hc'){
    
    k <- random_k(k, seed)
    
    set.seed(seed = seed)
    method <- sample(c('ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'), 1, F)
    
    clust_res <- hclust(d = dist(dataset), method = method)
    clust_res <- cutree(clust_res, k = k)
    
  }
  
  
  ###########################
  #Spectral clustering
  if (algo == 'spectral'){
    
    k <- random_k(k, seed)
    
    require(kknn)
    clust_res <- specClust(data = dataset, centers = k)$cluster
    
  }
  ###########################
  #Self-organizing maps
  
  return(clust_res)
}

############################
#merge results
create_adjacency_matrix <- function(CC_obj, patient_labels){
  clust_mat <- matrix(data = NA, nrow = length(patient_labels), ncol = length(CC_obj), dimnames = list(patient_labels, 1:1000))
  for (i in 1:length(CC_obj)){
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

# adj_mat <- lapply(1:100, function(i) {
# message(paste0('merging adjacency matrix # ', i, 'of 100'))
# create_adjacency_matrix(CC_obj = c(CC_PCA_Mclust , CC_PCA_PAM, CC_UMAP_Mclust, CC_UMAP_PAM)[[i]], patient_labels = row.names(metadata_AML_FH_ALL_normal)[1:1080])
# })


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



#################
#find patients that are not confidently classified in a particular clustering solution
find_silhouette_score_CC_obj <- function(norm_adj_mat, silhouette_threshold = 0.25){
  k_sweep <- sapply(2:40, function(i){
    print(i)
    pam_k_sweep <- pam(x = 1-norm_adj_mat, k = i, diss = T)
    return(mean(silhouette(pam_k_sweep)[,'sil_width']))
  })
  optimal_k <- which.max(k_sweep)
  optimal_k <- optimal_k + 1 #start counting at 2
  
  pam_optimal_k <- pam(x = 1-norm_adj_mat, k = optimal_k, diss = T)
  low_confidence_samples <- names(which(silhouette(pam_optimal_k)[,3] < silhouette_threshold))
  high_confidence_samples <- names(which(silhouette(pam_optimal_k)[,3] >= silhouette_threshold))
  
  return(list(silhouette_scores = silhouette(pam_optimal_k)[,3],
              low_confidence_samples = low_confidence_samples,
              high_confidence_samples = high_confidence_samples,
              optimal_k = optimal_k))
  
}

#########################3
#correct adjacency matrix by removing samples that have a low silhouette score in each clustering solution
#set adjacency matrix values all to NA for the patients with low silhouette scores

correct_adj_mat_by_sil_score <- function(adj_mat_obj, sil_score_obj){
  low_conf_samples <- sil_score_obj$low_confidence_samples
  adj_mat_obj <- adj_mat_obj
  
  pb <- txtProgressBar(min = 1, max = length(low_conf_samples), initial = 1, style = 3)
  for (low_conf_sample in low_conf_samples){
    setTxtProgressBar(pb = pb, value = which(low_conf_sample == low_conf_samples))
    for (sample in row.names(adj_mat_obj)){
        adj_mat_obj[sample,low_conf_sample] <- NA
    }
  }
  normalized_adj_mat <- t(t(adj_mat_obj)/sapply(1:ncol(adj_mat_obj), function(i) adj_mat_obj[i,i]))
  
  return(list(adj_mat = adj_mat_obj,
              normalized_adj_mat = normalized_adj_mat))
}

corr_adj_mat_CC_PHATE_LEIDEN <- correct_adj_mat_by_sil_score(adj_mat_obj = adj_mat_CC_PHATE_LEIDEN[[1]] , sil_score_obj = sil_score_CC_PHATE_LEIDEN[[1]])

corr_adj_mat_CC_PHATE_LEIDEN[1:5,sil_score_CC_PHATE_LEIDEN[[1]]$high_confidence_samples[1:5]]
colSums(corr_adj_mat_CC_PHATE_LEIDEN[,sil_score_CC_PHATE_LEIDEN[[1]]$low_confidence_samples])
colSums(corr_adj_mat_CC_PHATE_LEIDEN[,sil_score_CC_PHATE_LEIDEN[[1]]$high_confidence_samples])


corr_adj_mat_CC_PHATE_LEIDEN <- t(t(corr_adj_mat_CC_PHATE_LEIDEN)/sapply(1:ncol(corr_adj_mat_CC_PHATE_LEIDEN), function(i) corr_adj_mat_CC_PHATE_LEIDEN[i,i]))
pam_low_conf <- pam(x = 1-corr_adj_mat_CC_PHATE_LEIDEN, k = 15, diss = T)
table(pam_low_conf$clustering[sil_score_CC_PHATE_LEIDEN[[1]]$low_confidence_samples])/length(sil_score_CC_PHATE_LEIDEN[[1]]$low_confidence_samples)
table(pam_low_conf$clustering[sil_score_CC_PHATE_LEIDEN[[1]]$high_confidence_samples])/length(sil_score_CC_PHATE_LEIDEN[[1]]$high_confidence_samples)
pca_low_conf <- prcomp(corr_adj_mat_CC_PHATE_LEIDEN)
plot(pca_low_conf$x[,1:2], pch =19, col = ifelse(row.names(pca_low_conf$x) %in% sil_score_CC_PHATE_LEIDEN[[1]]$low_confidence_samples, 'red', 'green'))



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


