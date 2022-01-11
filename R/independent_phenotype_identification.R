


# Compute number of shared markers
shared_markers <- function(phen1,phen2){
  phen1 <- as.numeric(phen1)
  phen2 <- as.numeric(phen2)
  return(sum(phen1==phen2 & phen1 > -1 & phen2 > -1))
}

# Compute (shared_markers)/(len(phen1)+len(phen2))
compensated_shared_markers <- function(phen1,phen2){
  phen1 <- as.numeric(phen1)
  phen2 <- as.numeric(phen2)
  p1_have_maker <- phen1 > -1
  p2_have_maker <- phen2 > -1
  return(sum(phen1==phen2 & p1_have_maker & p2_have_maker )/(sum(p1_have_maker)+sum(p2_have_maker)))
}


#' Identify independent relevant phenotypes
#' 
#' It uses a network clustering process to identify most representative
#' phenotypes when several similar phenotypes are relevant.
#' 
#' @param phen_data A Data.Frame with Marker columns, sample columns, and (effect_size, p_value, log2foldChange) columns.
#' @param channel_data Data.Frame containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param n_phenotypes maximum number of phenotypes to be considered from \code{phen_data} filtered by lowest p-values. Default: 1000.
#' @param min_confidence Minimal confidence threshold to filter output. Default: 0.5.
#' @param n_threads Number of threads to be used. Default: 1.
#' 
#' @export
get_independent_relevant_phenotypes <- function(phen_data,
                                                channel_data,
                                                n_phenotypes = 1000,
                                                min_confidence = 0.5,
                                                n_threads = 1
){
  
  markers <- channel_data$Marker
  
  sample_ids <- colnames(phen_data)[(length(markers)+1):(ncol(phen_data)-3)]
  

  if(n_phenotypes < nrow(phen_data)){
    # Sort phenotypes by abs(effect size) and p-val and keep first n_phenotype ones
    print_log("Getting first ",n_phenotypes, " most relevant phenotypes...")
    phen_order_factor <- -(normalize_vals(unlist(phen_data[,"log2foldChange"]))*normalize_vals(unlist(phen_data[,"effect_size"])))
    phen_data <- phen_data[order(phen_order_factor),]
    phen_data <- phen_data[1:n_phenotypes,]
    rownames(phen_data) <- 1:n_phenotypes
    
  }
  
  
  
  phen_data$n_markers <- unlist(parallel::mclapply(1:nrow(phen_data), function(i) count_markers(as.numeric(phen_data[i, markers])), mc.cores = n_threads))
  
  phen_order_factor <- -(normalize_vals(unlist(phen_data[,"log2foldChange"]))*normalize_vals(unlist(phen_data[,"effect_size"])))
  
  
  # Generate adjacency matrix based on marker distance
  print_log("Generating network...")
  phen_dist_list <- expand.grid(rownames(phen_data),rownames(phen_data))
  phen_dist_list[,1] <- as.numeric(phen_dist_list[,1])
  phen_dist_list[,2] <- as.numeric(phen_dist_list[,2])
  phen_dist_list$compensated_dist <- parallel::mclapply(1:nrow(phen_dist_list), function(i) compensated_shared_markers(phen_data[as.numeric(phen_dist_list[i,1]), markers],phen_data[as.numeric(phen_dist_list[i,2]), markers]), mc.cores = n_threads)
  phen_dist_matrix <- matrix(unlist(phen_dist_list$compensated_dist), nrow = nrow(phen_data))
  diag(phen_dist_matrix) <- 0
  
  rm(phen_dist_list)
  
  
  print_log("Finding independent phenotypes...")
  thresholds <- (1:11)
  thresholds <- thresholds/max(thresholds)
  thresholds <- head(thresholds*max(phen_dist_matrix),10)
  
  cluster_top_phenotypes <- unlist(parallel::mclapply(thresholds, function(t){
    
    local_phen_dist_matrix <- phen_dist_matrix
    local_phen_dist_matrix[phen_dist_matrix<t] <- 0
    
    
    #Generate network and cluster
    
    g <- igraph::graph_from_adjacency_matrix(local_phen_dist_matrix,mode = "undirected",weighted = TRUE,diag = FALSE)
    cluster <- igraph::cluster_louvain(g)
    
    
    
    
    # Assign clusters
    cluster <- cluster$membership
    
    # Order clusters by effect size and n markers
    phen_order <- order(cluster,phen_order_factor,unlist(phen_data$n_markers))
    
    
    # Count phenotypes in each cluster
    cluster_count <- as.data.frame(table(cluster))
    cluster_count$cumFreq <- cumsum(cluster_count[,"Freq"])
    #cluster_count <- cluster_count[cluster_count$Freq>1,]
    
    if(nrow(cluster_count)){
      final_phenotypes <- phen_order[c(1,cluster_count[1:(nrow(cluster_count)-1),"cumFreq"]+1)] 
      return(as.list(final_phenotypes))
    }else{
      return(list())
    }
    
  }, mc.cores = n_threads))
  
  
  cluster_top_phenotypes <- as.data.frame(table(cluster_top_phenotypes))
  colnames(cluster_top_phenotypes) <- c("Phenotype_ID","Count")
  
  cluster_top_phenotypes$Confidence <- cluster_top_phenotypes$Count / length(thresholds)
  
  cluster_top_phenotypes$Phenotype_ID <- as.numeric(levels(cluster_top_phenotypes$Phenotype_ID))[cluster_top_phenotypes$Phenotype_ID]
  
  
  
  if(min_confidence > 0.0){
    print_log("Independent phenotypes found: ", nrow(cluster_top_phenotypes))
    print_log("Filtering phenotypes with confidence equal or above ",min_confidence)
    cluster_top_phenotypes <- cluster_top_phenotypes[cluster_top_phenotypes$Confidence >= min_confidence,]  
  }
  
  
  
  cluster_top_phenotypes <- cluster_top_phenotypes[order(-cluster_top_phenotypes$Confidence),]
  
  print_log("Final independent phenotypes: ", nrow(cluster_top_phenotypes))
  
  
  cluster_top_phenotypes$Phenotype <- unlist(parallel::mclapply(cluster_top_phenotypes$Phenotype_ID, function(i) make_phenotype_name(as.numeric(phen_data[i, markers]),markers), mc.cores = n_threads))
  
  
  final_phenotypes <- phen_data[as.numeric(cluster_top_phenotypes$Phenotype_ID),(length(markers)+1):ncol(phen_data)]
  final_phenotypes$Confidence <- cluster_top_phenotypes$Confidence
  final_phenotypes$Phenotype <- cluster_top_phenotypes$Phenotype
  
  final_phenotypes <- final_phenotypes[, c("Phenotype","n_markers","effect_size","log2foldChange","p_value","Confidence",sample_ids)]
  
  return(final_phenotypes)
  
  
}

