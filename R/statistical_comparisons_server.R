

# Filter statistically relevant phenotypes given a comparison between two groups.
# Must be called from "statistically_relevant_phenotypes_server"
memory_safe_compute_statistically_relevant_phenotypes <- function(output_folder,
                                                                  sample_data,
                                                                  channel_data,
                                                                  input_phenotype_counts,
                                                                  output_file,
                                                                  test_type = "group", # Options: "group", "correlation", "survival"
                                                                  groups_column = NULL,
                                                                  g1 = NULL,
                                                                  g2 = NULL,
                                                                  correlation_column = NULL,
                                                                  survival_time_column = NULL,
                                                                  survival_status = NULL,
                                                                  max_pval = 0.05,
                                                                  parent_phen = NULL,
                                                                  start_from = 0,
                                                                  n_threads = 1
){
  
  
  print_log("Starting statistical filtering...")
  
  markers <- channel_data[,"Marker"]
  n_markers <- length(markers)
  
  sample_ids <- as.character(sample_data$Sample_ID)
  
  phenotype_counts_file_path <- file.path(output_folder,input_phenotype_counts)
  
  filtered_phenotypes_file <- output_file
  
  filtered_phenotypes_file_path <- file.path(output_folder,filtered_phenotypes_file)
  
  
  print_log("Getting number of phenotypes in file...")
  total_lines <- count_csv_lines(phenotype_counts_file_path)
  print_log("Number of phenotypes found: ",total_lines)  
  
  total_cells <- data.table::fread(phenotype_counts_file_path, nrows = 2, header= T)
  header <- colnames(total_cells)
  
  column_types <- as.vector(sapply(total_cells, typeof))
  column_types[column_types=="character"] <- "string"
  
  total_cells <- total_cells[1 ,]
  
  
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
      total_cells <- parent_phen_data[1,]
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
    
    
    if(!is.null(parent_phen_data)){
      
      print_log("Filtering for parent population...")
      
      has_parent_phen <- unlist(parallel::mclapply(1:nrow(pheno_counts), function(i) has_phenotype(pheno_counts[i, markers],parent_phen), mc.cores = n_threads))
      
      pheno_counts <- pheno_counts[has_parent_phen, ]
    }
    
    
    print_log("Computing frequencies...")
    
    pheno_counts[,sample_ids] <- count_to_frequency(pheno_counts[,sample_ids], total_cells[,sample_ids], n_threads)
    
    print_log("Computing statistical test...")
    
    if(test_type == "group"){
      
      pheno_counts <- statistical_test_group(pheno_counts, sample_data, groups_column, g1, g2, n_threads = n_threads)
      
    }else if(test_type == "correlation"){
      
      pheno_counts <- statistical_test_correlation(pheno_counts, sample_data, correlation_column, n_threads = n_threads)
      
    }
    
    print_log("Filtering statistically relevant phenotypes...")
    
    pval_filter <- unlist(parallel::mclapply(1:nrow(pheno_counts), function(i) pheno_counts[i,"p_value"] <= max_pval, mc.cores = n_threads))
    
    pheno_counts <- pheno_counts[pval_filter, ]
    
    append_file <- i > 1 | start_from > 0
    
    print_log("Writing ",nrow(pheno_counts), " significant phenotypes to file...")
    
    if(nrow(pheno_counts)){
      data.table::fwrite(pheno_counts,filtered_phenotypes_file_path,append = append_file)  
    }
    
    print_log(format((last_line/total_lines)*100, digits = 2),"% done..." )
    
    rm(pheno_counts,pval_filter)
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
                                                     input_phenotype_counts = "combinatorial_phenotype_counts.csv",
                                                     output_file = "significant_phenotypes.csv",
                                                     log_file = "significant_phenotypes.log",
                                                     test_type = "group", # Options: "group", "correlation", "survival"
                                                     groups_column = NULL,
                                                     g1 = NULL,
                                                     g2 = NULL,
                                                     correlation_column = NULL,
                                                     survival_time_column = NULL,
                                                     survival_status = NULL,
                                                     max_pval = 0.05,
                                                     parent_phen = NULL,
                                                     continue = FALSE,
                                                     n_threads = 1,
                                                     verbose = FALSE){
  
  continued <- FALSE
  
  phenotype_counts_file <- input_phenotype_counts
  filtered_phenotypes_file <- output_file
  
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
                                                                  input_phenotype_counts,
                                                                  output_file,
                                                                  test_type, # Options: "group", "correlation", "survival"
                                                                  groups_column,
                                                                  g1,
                                                                  g2,
                                                                  correlation_column,
                                                                  survival_time_column,
                                                                  survival_status,
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
                                                          input_phenotype_counts,
                                                          output_file,
                                                          test_type, # Options: "group", "correlation", "survival"
                                                          groups_column,
                                                          g1,
                                                          g2,
                                                          correlation_column,
                                                          survival_time_column,
                                                          survival_status,
                                                          max_pval,
                                                          parent_phen,
                                                          0,
                                                          n_threads
    )
    
    print_log("Statistical filtering of phenotypes done.")
    
    stop_log() 
    
  }
  
}

