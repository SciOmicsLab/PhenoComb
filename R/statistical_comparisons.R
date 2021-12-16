

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
#' @param n_threads Number of threads to be used. Default: 1
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
  
  g1_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g1,"Sample_ID"])
  g2_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g2,"Sample_ID"])
  
  markers <- channel_data[,"Marker"]
  n_markers <- length(markers)
  
  sample_dt <- phenotype_cell_counts[ , c(g1_IDs,g2_IDs)]
  
  per_sample_total_counts <- sample_dt[1,]
  
  #Filter for parent phenotype
  if(!is.null(parent_phen)){
    
    # must be in numeric encoding
    if(is.character(parent_phen)){
      parent_phen <- phenotype_to_numbers(parent_phen,markers)  
    }
    
    
    parent_phen_data <- find_phenotype(phenotype_cell_counts, parent_phen, markers, n_threads)
    
    if(is.null(parent_phen_data)){
      warning("Parent penotype not found in data. Not filtering...")
    }else{
      
      
      has_parent_phen <- unlist(parallel::mclapply(1:nrow(phenotype_cell_counts), function(i) has_phenotype(phenotype_cell_counts[i, markers], parent_phen), mc.cores = n_threads))
      
      sample_dt <- sample_dt[has_parent_phen , ]
      per_sample_total_counts <- parent_phen_data[ , c(g1_IDs,g2_IDs)]
      
      phenotype_cell_counts <- phenotype_cell_counts[has_parent_phen, ]
      
    }
    
  }
  
  sample_dt <- count_to_frequency(sample_dt, per_sample_total_counts, n_threads)
  
  st_test <- statistical_test(sample_dt[, g1_IDs],sample_dt[, g2_IDs],n_threads)
  
  final_significant_phens <- cbind(phenotype_cell_counts[, markers],sample_dt,st_test)
  
  pval_filter <- c(FALSE,unlist(parallel::mclapply(2:nrow(final_significant_phens), function(i) final_significant_phens[i,"p_value"] <= max_pval, mc.cores = n_threads)))
  
  final_significant_phens <- final_significant_phens[pval_filter, ]
  
  if(nrow(final_significant_phens)){
    final_significant_phens$log2foldChange <- get_log2foldChange(final_significant_phens[ , g1_IDs], final_significant_phens[ , g2_IDs], n_threads)
  }
  
  
  rm(phenotype_cell_counts,st_test,sample_dt)
  gc(verbose = F, full = T)
  
  
  return(final_significant_phens)
}




# Filter statistically relevant phenotypes given a comparison between two groups.
# Must be called from "statistically_relevant_phenotypes_server"
memory_safe_compute_statistically_relevant_phenotypes <- function(output_folder,
                                                                  sample_data,
                                                                  channel_data,
                                                                  groups_column,
                                                                  g1,
                                                                  g2,
                                                                  max_pval = 0.05,
                                                                  parent_phen = NULL,
                                                                  start_from = 0,
                                                                  n_threads = 1
){
  
  
  print_log("Starting statistical filtering...")
  
  print_log("Getting sample groups...")
  
  g1_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g1,"Sample_ID"])
  g2_IDs <- as.character(sample_data[sample_data[,groups_column] %in% g2,"Sample_ID"])
  
  markers <- channel_data[,"Marker"]
  n_markers <- length(markers)
  
  phenotype_counts_file <- "combinatorial_phenotype_counts.csv"
  
  phenotype_counts_file_path <- file.path(output_folder,phenotype_counts_file)
  
  filtered_phenotypes_file <- paste("significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,".csv",collapse = "",sep="")
  
  filtered_phenotypes_file_path <- file.path(output_folder,filtered_phenotypes_file)
  
  
  print_log("Getting number of phenotypes in file...")
  total_lines <- count_csv_lines(phenotype_counts_file_path)
  print_log("Number of phenotypes found: ",total_lines)  
  
  total_cells <- data.table::fread(phenotype_counts_file_path, nrows = 2, header= T)
  header <- colnames(total_cells)
  
  column_types <- as.vector(sapply(total_cells, typeof))
  column_types[column_types=="character"] <- "string"
  
  total_cells <- total_cells[1 , c(g1_IDs,g2_IDs)]
  
  
  #Filter for parent phenotype
  parent_phen_data <- NULL
  if(!is.null(parent_phen)){
    
    print_log("Looking for parent population ", parent_phen)
    
    parent_phen <- phenotype_to_numbers(parent_phen, markers)
    
    parent_phen_data <- find_phenotype_in_file(phenotype_counts_file_path, parent_phen, markers, n_threads)
    
    if(is.null(parent_phen_data)){
      print_log("Parent penotype not found in data. Not filtering...")
    }else{
      print_log("Parent penotype found.")
      total_cells <- parent_phen_data[ , c(g1_IDs,g2_IDs)]
    }
    
  }
  
  
  chunk_size <- n_threads*10000
  
  if(chunk_size > total_lines){
    chunk_size <- total_lines
  }
  
  n_chunks <- ceiling((total_lines-start_from)/chunk_size)
  
  
  print_log("Dividing phenotypes into ",format(n_chunks, scientific = F)," chunks of ",format(chunk_size, scientific = F)," phenotypes each using ",n_threads," thread(s)...")
  
  if(start_from > 0){
    print_log("Finding last phenotype computed in file...")
  }
  
  laf <- LaF::laf_open_csv(phenotype_counts_file_path, column_types = column_types, skip = start_from+2)
  
  print_log("--------------------------------------------")
  
  for(i in 1:n_chunks){
    
    
    
    last_line <- (i)*chunk_size+start_from
    if(last_line > total_lines){
      last_line <- total_lines
    }
    
    start_line <- (i-1)*chunk_size+start_from+2
    
    print_log("Reading phenotypes from ",format(start_line, scientific = F), " to ",format(last_line, scientific = F))
    
    pheno_counts <- LaF::next_block(laf,nrows=chunk_size)
    colnames(pheno_counts) <- header
    
    sample_dt <- pheno_counts[ , c(g1_IDs,g2_IDs)]
    
    
    if(!is.null(parent_phen_data)){
      
      print_log("Filtering for parent population...")
      
      has_parent_phen <- unlist(parallel::mclapply(1:nrow(pheno_counts), function(i) has_phenotype(pheno_counts[i, markers],parent_phen), mc.cores = n_threads))
      
      sample_dt <- sample_dt[has_parent_phen , ]
      pheno_counts <- pheno_counts[has_parent_phen, ]
    }
    
    
    print_log("Computing frequencies...")
    
    sample_dt <- count_to_frequency(sample_dt,total_cells,n_threads)
    
    print_log("Computing statistical test...")
    
    st_test <- statistical_test(sample_dt[, g1_IDs],sample_dt[, g2_IDs],n_threads)
    
    final_significant_phens <- cbind(pheno_counts[, markers], sample_dt, st_test)
    
    print_log("Filtering statistically relevant phenotypes...")
    
    pval_filter <- unlist(parallel::mclapply(1:nrow(final_significant_phens), function(i) final_significant_phens[i,"p_value"] <= max_pval, mc.cores = n_threads))
    
    final_significant_phens <- final_significant_phens[pval_filter, ]
    
    if(nrow(final_significant_phens)){
      final_significant_phens$log2foldChange <- get_log2foldChange(final_significant_phens[ , g1_IDs], final_significant_phens[ , g2_IDs], n_threads)
    }
    
    append_file <- i > 1 | start_from > 0
    
    print_log("Writing ",nrow(final_significant_phens), " significant phenotypes to file...")
    
    if(nrow(final_significant_phens)){
      data.table::fwrite(final_significant_phens,filtered_phenotypes_file_path,append = append_file)  
    }
    
    print_log(format((last_line/total_lines)*100, digits = 2),"% done..." )
    
    rm(final_significant_phens,pheno_counts,st_test,sample_dt)
    gc(verbose = F, full = T)
    
    print_log("--------------------------------------------")
    
  }
  
}


# Find last phenotype computed to resume computation
find_last_phenotype_filtered <- function(log_file){
  
  last_phenotype <- 0
  
  con = file(log_file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(grepl("Reading phenotypes from ", line, fixed = TRUE)){
      last_phenotype <- as.integer(stringr::str_extract(line, "(?<= Reading phenotypes from )[0-9]*"))
    }
  }
  
  close(con)
  
  return(last_phenotype)
}


#' Filter statistically relevant phenotypes by comparing two sample groups (Server Version)
#' 
#' \code{compute_statistically_relevant_phenotypes} filters statistically
#' relevant phenotypes by comparing two sample groups defined in \code{sample_data}.
#' 
#' This version reads the "output_folder/combinatorial_phenotype_counts.csv" by chunks
#' and save the partial results to a file. It is intended to be used with large datasets
#' which inputs wouldn't fit in the available memory.
#' 
#' @param output_folder Path to folder where output files from this and previous steps should be saved.
#' @param channel_file Path to a ".csv" file containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param sample_file Path to a ".csv" file containing a Sample_ID column and additional grouping columns for the samples.
#' @param groups_column Column name in \code{sample_data} where group identifications are stored.
#' @param g1 Group label or vector with group labels for first group.
#' @param g2 Group label or vector with group labels for second group.
#' @param max_pval Maximum p-value. Used to filter phenotypes in final output.
#' @param parent_phen Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.
#' @param continue If TRUE, look for files to resume execution. Also needed to save necessary "continuing" files.
#' @param n_threads Number of threads to be used. Default: 1.
#' @param verbose If TRUE, print outputs from log to stdout.
#'  
#' @return Output is saved to file.
#' 
#' @export
statistically_relevant_phenotypes_server <- function(output_folder,
                                                     channel_file,
                                                     sample_file,
                                                     groups_column,
                                                     g1,
                                                     g2,
                                                     max_pval = 0.05,
                                                     parent_phen = NULL,
                                                     continue = FALSE,
                                                     n_threads = 1,
                                                     verbose = FALSE){
  
  continued <- FALSE
  
  phenotype_counts_file <- "combinatorial_phenotype_counts.csv"
  log_file <- paste("significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,"_log.txt",collapse = "",sep="")
  filtered_phenotypes_file <- paste("significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,".csv",collapse = "",sep="")
  
  channel_data <- as.data.frame(data.table::fread(channel_file,check.names = FALSE))
  sample_data <- as.data.frame(data.table::fread(sample_file,check.names = FALSE))
  
  
  if(continue){
    
    if(dir.exists(output_folder)){
      
      
      # Get files in output folder
      folder_files <- list.files(output_folder, full.names = F, recursive = F)
      
      if(length(folder_files) > 0){# If there are files in the output_folder
        
        # Check if all continuation files exist
        log_file_exists <- log_file %in% folder_files
        filtered_phenotypes_file_exists <- filtered_phenotypes_file %in% folder_files
        
        
        if(all(c(log_file_exists,
                 filtered_phenotypes_file_exists)))
        { # If all continuation files exist
          
          # Check for chunk information on log file
          last_phenotype_filtered <- find_last_phenotype_filtered(file.path(output_folder,log_file))
          
          
          if(last_phenotype_filtered){ #If chunk information found
            
            continued <- TRUE
            
            start_log(log_path = file.path(output_folder,log_file), append = TRUE, verbose = verbose)
            
            print_log("Continuation files found. Resuming...")
            
            
            
            memory_safe_compute_statistically_relevant_phenotypes(output_folder,
                                                                  sample_data,
                                                                  channel_data,
                                                                  groups_column,
                                                                  g1,
                                                                  g2,
                                                                  max_pval,
                                                                  parent_phen,
                                                                  last_phenotype_filtered,
                                                                  n_threads
            )
            
            print_log("Statistical filtering of phenotypes done.")
            
            stop_log() 
            
          }else{
            
            message("No chunk information found. Starting from scratch...")
            
          }
          
          
          
          # Do continuing stuff
          
        }else{ # If log file does not exist
          message("Continuation files not found. Starting from scratch...")
        }
        
      }else{# if no files in output_folder
        message("Continuation files not found. Starting from scratch...")
      }
      
    }else{ # if output_dir does not exist
      continued <- TRUE # To skip next if
      message("Folder ", output_folder, " does not exist. Aborting...")
    }
    
  }
  
  if(!continued){
    
    
    start_log(log_path = file.path(output_folder,log_file), append = FALSE, verbose = verbose)
    
    memory_safe_compute_statistically_relevant_phenotypes(output_folder,
                                                          sample_data,
                                                          channel_data,
                                                          groups_column,
                                                          g1,
                                                          g2,
                                                          max_pval,
                                                          parent_phen,
                                                          0,
                                                          n_threads
    )
    
    print_log("Statistical filtering of phenotypes done.")
    
    stop_log() 
    
  }
  
}


