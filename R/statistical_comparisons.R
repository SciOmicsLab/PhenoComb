

# Transform cell count to frequency in locus
count_to_frequency <- function(samples, total_counts, n_threads = 1){
  
  sample_names <- colnames(samples)
  
  samples <- parallel::mclapply(sample_names, function(s) samples[,s] / as.numeric(total_counts[1 , s]), mc.cores = n_threads)
  
  samples <- as.data.frame(do.call(cbind,samples))
  
  colnames(samples) <- sample_names
  
  return(samples)

}

# Performs Mann-Whitney U test, returns vector with transformed effect size and p-value
mann_whitney_u_test <- function(g1,g2,n_comparisons){
  g1 <- as.numeric(g1)
  g2 <- as.numeric(g2)
  ut <- wilcox.test(g1,g2)
  return(c(((ut$statistic/n_comparisons)-0.5)*2,ut$p.value))
}

# Performs statistical comparision for each row of g1 and g2
statistical_test <- function(g1, g2, n_threads = 1){
  
  n_comparisons <- ncol(g1)*ncol(g2)
  
  st_test <- as.matrix(do.call(rbind,parallel::mclapply(1:nrow(g1), function(i) mann_whitney_u_test(g1[i,], g2[i,], n_comparisons), mc.cores = n_threads)))
  
  colnames(st_test) <- c("effect_size","p_value")
  
  st_test[is.nan(st_test)] <- 1.
  
  return(st_test)
}

# Compute the log2foldChange between two vectors
log2foldChange <- function(a,b){
  a <- as.numeric(a)
  b <- as.numeric(b)
  return(log2(mean(a)/mean(b)))
}

# Compute the log2foldChange for each row of g1 and g2
get_log2foldChange <- function(g1, g2, n_threads = 1){

  l2fc <- unlist(parallel::mclapply(1:nrow(g1), function(i) log2foldChange(g1[i,],g2[i,]), mc.cores = n_threads))
  
  return(l2fc)
}


#' Filter statistically relevant phenotypes by comparing two sample groups
#' 
#' \code{compute_statistically_relevant_phenotypes} filters statistically
#' relevant phenotypes by comparing two sample groups defined in \code{sample_data}.
#' 
#' @param phenotype_cell_counts Data.Frame object with marker and sample columns. Each row
#' is a unique phenotype.
#' @param channel_data Data.Frame containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param sample_data Data.Frame containing a Sample_ID column and additional grouping columns for the samples.
#' @param groups_column Column name in \code{sample_data} where group identifications are stored.
#' @param g1 Group label or vector with group labels for first group.
#' @param g2 Group label or vector with group labels for second group.
#' @param max_pval Maximum p-value. Used to filter phenotypes in final output.
#' @param parent_phen Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.
#' @param n_threads Number of threads to be used. Default: 1.
#' 
#' @export
compute_statistically_relevant_phenotypes <- function(phenotype_cell_counts,
                                                      channel_data,
                                                      sample_data,
                                                      groups_column,
                                                      g1,
                                                      g2,
                                                      max_pval = 0.05,
                                                      parent_phen = NULL,
                                                      n_threads = 1
){
  
  print_log("Starting statistical filtering...")
  
  print_log("Getting sample groups...")
  
  g1_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g1,"Sample_ID"])
  g2_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g2,"Sample_ID"])
  
  markers <- channel_data[,"Marker"]
  n_markers <- length(markers)
  
  sample_dt <- phenotype_cell_counts[ , c(g1_IDs,g2_IDs)]
  
  per_sample_total_counts <- sample_dt[1,]
  
  #Filter for parent phenotype
  if(!is.null(parent_phen)){
    
    print_log("Looking for parent population ", parent_phen)
    
    # must be in numeric encoding
    if(is.character(parent_phen)){
      parent_phen <- phenotype_to_numbers(parent_phen,markers)  
    }
    
    
    parent_phen_data <- find_phenotype(phenotype_cell_counts, parent_phen, markers, n_threads)
    
    if(is.null(parent_phen_data)){
      print_log("Parent penotype not found in data. Not filtering...")
    }else{
      
      print_log("Parent penotype found. Filtering phenotypes...")
      
      has_parent_phen <- unlist(parallel::mclapply(1:nrow(phenotype_cell_counts), function(i) has_phenotype(phenotype_cell_counts[i, markers], parent_phen), mc.cores = n_threads))
      
      sample_dt <- sample_dt[has_parent_phen , ]
      per_sample_total_counts <- parent_phen_data[ , c(g1_IDs,g2_IDs)]
      
      phenotype_cell_counts <- phenotype_cell_counts[has_parent_phen, ]
      
    }
    
  }
  
  print_log("Computing frequencies...")
  sample_dt <- count_to_frequency(sample_dt, per_sample_total_counts, n_threads)
  
  print_log("Computing statistical test...")
  st_test <- statistical_test(sample_dt[, g1_IDs],sample_dt[, g2_IDs],n_threads)
  
  final_significant_phens <- cbind(phenotype_cell_counts[, markers],sample_dt,st_test)
  
  print_log("Filtering statistically relevant phenotypes with p-value <= ",max_pval)
  
  pval_filter <- c(FALSE,unlist(parallel::mclapply(2:nrow(final_significant_phens), function(i) final_significant_phens[i,"p_value"] <= max_pval, mc.cores = n_threads)))
  
  final_significant_phens <- final_significant_phens[pval_filter, ]
  
  if(nrow(final_significant_phens)){
    print_log("Computing log2foldChanges...")
    final_significant_phens$log2foldChange <- get_log2foldChange(final_significant_phens[ , g1_IDs], final_significant_phens[ , g2_IDs], n_threads)
  }
  
  print_log("Final statistically relevant phenotypes: ",nrow(final_significant_phens))
  
  rm(phenotype_cell_counts,st_test,sample_dt)
  gc(verbose = F, full = T)
  
  
  return(final_significant_phens)
}



