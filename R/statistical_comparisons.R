

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
statistical_test_group <- function(counts_data, sample_data, groups_column, g1, g2, n_threads = 1){
  
  g1_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g1,"Sample_ID"])
  g2_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g2,"Sample_ID"])
  
  g1 <- counts_data[,g1_IDs]
  g2 <- counts_data[,g2_IDs]
  
  n_comparisons <- ncol(g1)*ncol(g2)
  
  st_test <- as.matrix(do.call(rbind,parallel::mclapply(1:nrow(g1), function(i) mann_whitney_u_test(g1[i,], g2[i,], n_comparisons), mc.cores = n_threads)))
  
  colnames(st_test) <- c("effect_size","p_value")
  
  st_test[is.nan(st_test)] <- 1.
  
  st_test <- as.data.frame(st_test)
  
  counts_data <- cbind(counts_data,st_test)
  
  counts_data$log2foldChange <- get_log2foldChange(counts_data[, g1_IDs], counts_data[, g2_IDs], n_threads)
  
  rm(st_test)
  gc(full = TRUE,verbose = FALSE)
  
  return(counts_data)
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


kendall_correlation_test <- function(x,y){
  if(sd(x) == 0) return(c(0.,1.))
  x <- as.numeric(x)
  y <- as.numeric(y)
  cor_test <- cor.test(x,y,method = "kendall")
  return(c(cor_test$estimate,cor_test$p.value))
}

statistical_test_correlation <- function(counts_data, sample_data, correlation_column, n_threads = 1){
  
  sample_ids <- as.character(sample_data$Sample_ID)
  
  correlation_data <- sample_data[,correlation_column]
  
  counts_data_list <- as.matrix(counts_data[,sample_ids])
  
  counts_data_list <- parallel::mclapply(seq_len(nrow(counts_data_list)), function(i) counts_data_list[i,], mc.cores = n_threads)
  
  st_test <- as.matrix(do.call(rbind,parallel::mclapply(counts_data_list, function(row) kendall_correlation_test(row, correlation_data), mc.cores = n_threads)))
  
  rm(counts_data_list)
  
  colnames(st_test) <- c("correlation","p_value")
  
  st_test[is.nan(st_test)] <- 1.
  
  st_test <- as.data.frame(st_test)
  
  counts_data <- cbind(counts_data,st_test)
  
  rm(st_test)
  gc(full = TRUE,verbose = FALSE)
  
  return(counts_data)
  
}


survival_test <- function(correlates, survival_times, status){
  
  correlates <- as.numeric(correlates)
  
  test_data <- list(time = survival_times, 
                status = status, 
                x = correlates)
  
  surv_test <- summary(survival::coxph(survival::Surv(time, status) ~ x, test_data))

  rm(test_data)
  
  return(c(surv_test$coefficients[1],surv_test$sctest[3]))
  
}


statistical_test_survival <- function(counts_data, sample_data, survival_time_column, survival_status_column, n_threads = 1){
  
  sample_ids <- as.character(sample_data$Sample_ID)
  
  survival_time_data <- as.numeric(sample_data[,survival_time_column])
  
  survival_status_data <- as.numeric(sample_data[,survival_status_column])
  
  counts_data_list <- as.matrix(counts_data[,sample_ids])
  
  counts_data_list <- parallel::mclapply(seq_len(nrow(counts_data_list)), function(i) counts_data_list[i,], mc.cores = n_threads)
  
  st_test <- as.matrix(do.call(rbind,parallel::mclapply(counts_data_list, function(row) survival_test(row, survival_time_data, survival_status_data), mc.cores = n_threads)))
  
  rm(counts_data_list)
  
  colnames(st_test) <- c("coxph_coefficient","p_value")
  
  st_test[is.nan(st_test)] <- 1.
  
  st_test <- as.data.frame(st_test)
  
  counts_data <- cbind(counts_data,st_test)
  
  rm(st_test)
  gc(full = TRUE,verbose = FALSE)
  
  return(counts_data)
  
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
#' @param test_type Type of statistical test to be performed. Value can be "group", "correlation", or "survival". Default: "group".
#' Additional parameters should be provided accordingly, such as [groups_column, g1, g2] for "group", [correlation_column] for "correlation", 
#' and [survival_time_column, survival_status_column] for "survival". Parameters not used in the test will be ignored.
#' @param groups_column Column name in \code{sample_data} where group identifications are stored.
#' @param g1 Group label or vector with group labels for first group.
#' @param g2 Group label or vector with group labels for second group.
#' @param correlation_column Column name in \code{sample_data} where data to be correlated is stored.
#' @param survival_time_column Column name in \code{sample_data} where survival time is stored.
#' @param survival_status_column Column name in \code{sample_data} where survival status is stored.
#' @param max_pval Maximum p-value. Used to filter phenotypes in final output.
#' @param parent_phen Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.
#' @param n_threads Number of threads to be used. Default: 1.
#' 
#' @export
compute_statistically_relevant_phenotypes <- function(phenotype_cell_counts,
                                                      channel_data,
                                                      sample_data,
                                                      test_type = "group", # Options: "group", "correlation", "survival"
                                                      groups_column = NULL,
                                                      g1 = NULL,
                                                      g2 = NULL,
                                                      correlation_column = NULL,
                                                      survival_time_column = NULL,
                                                      survival_status_column = NULL,
                                                      max_pval = 0.05,
                                                      parent_phen = NULL,
                                                      n_threads = 1
){
  
  print_log("Starting statistical filtering...")
  
  markers <- channel_data[,"Marker"]
  n_markers <- length(markers)
  
  sample_ids <- as.character(sample_data$Sample_ID)
  
  per_sample_total_counts <- phenotype_cell_counts[1,]
  
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
      
      phenotype_cell_counts <- phenotype_cell_counts[has_parent_phen, ]
      
      per_sample_total_counts <- parent_phen_data
      
    }
    
  }
  
  print_log("Computing frequencies...")
  phenotype_cell_counts[,sample_ids] <- count_to_frequency(phenotype_cell_counts[,sample_ids], per_sample_total_counts[,sample_ids], n_threads)
  
  
  print_log("Computing statistical test...")
  
  if(test_type == "group"){
    
    phenotype_cell_counts <- statistical_test_group(phenotype_cell_counts, sample_data, groups_column, g1, g2, n_threads = n_threads)
    
  }else if(test_type == "correlation"){
    
    phenotype_cell_counts <- statistical_test_correlation(phenotype_cell_counts, sample_data, correlation_column, n_threads = n_threads)
    
  }else if(test_type == "survival"){
    
    phenotype_cell_counts <- statistical_test_survival(phenotype_cell_counts, sample_data, survival_time_column, survival_status_column, n_threads = n_threads)
    
  }else{
    print_log("Statistical test type not valid. test_type should be in c(group, correlation, survival).")
    print_log("Aborting...")
    stop("Invalid value for test_type.")
  }
  
  
  print_log("Filtering statistically relevant phenotypes with p-value <= ",max_pval)

  pval_filter <- unlist(parallel::mclapply(1:nrow(phenotype_cell_counts), function(i) phenotype_cell_counts[i,"p_value"] <= max_pval, mc.cores = n_threads))

  phenotype_cell_counts <- phenotype_cell_counts[pval_filter, ]

  print_log("Final statistically relevant phenotypes: ",nrow(phenotype_cell_counts))


  return(phenotype_cell_counts)
}







