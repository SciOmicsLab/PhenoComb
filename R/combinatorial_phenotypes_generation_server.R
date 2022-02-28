

# Split sequence of numbers in chunks and return chunks in a list
split_in_chunks <- function(n_rows, chunk_size){
  
  if(chunk_size >= n_rows) return(list(c(1:n_rows)))
  
  n_chunks <- ceiling(n_rows / chunk_size)
  
  chunk_list <- list()
  
  for (i in 1:(n_chunks - 1)){
    chunk_list[[i]] <- (((i - 1) * chunk_size) + 1):(i * chunk_size)
  }
  
  chunk_list[[n_chunks]] <- (((n_chunks - 1) * chunk_size) + 1):n_rows
  
  return(chunk_list)
  
}


# Compute all combinations of markers in a memory-safe way by
# splitting all the combinations in chunks and saving the results to a file
# must be called from "combinatorial_phenotype_counts_server"
memory_safe_combinatorial_phenotype_counts <- function(unique_phenotype_counts,
                                                       n_markers,
                                                       dump_file,
                                                       parent_phen = NULL,
                                                       max_phenotype_length = 0,
                                                       min_count = 10,
                                                       sample_fraction_min_counts = 0.0,
                                                       start_from = 1,
                                                       max_ram = 0,
                                                       efficient = TRUE,
                                                       max_empty_chunks = 1000,
                                                       n_threads = 1
){
  
  
  
  # Get markers names and sample IDs
  markers <- head(colnames(unique_phenotype_counts),n_markers)
  samples_id <- tail(colnames(unique_phenotype_counts),ncol(unique_phenotype_counts)-n_markers)
  
  n_samples <- length(samples_id)
  
  if(sample_fraction_min_counts <= 1/n_samples){
    sample_fraction_min_counts <- 1/n_samples
  }
  
  total_combinations <- 2^n_markers
  
  # Select parent population
  parent_markers <- NULL
  if(!is.null(parent_phen)){
    
    print_log("Filtering for parent population ", parent_phen)
    
    parent_phen <- phenotype_to_numbers(parent_phen, markers)
    
    has_parent_phen <- unlist(parallel::mclapply(1:nrow(unique_phenotype_counts), function(i) has_phenotype(unique_phenotype_counts[i,markers], parent_phen), mc.cores = n_threads))
    
    unique_phenotype_counts <- unique_phenotype_counts[has_parent_phen,]
    
    if(nrow(unique_phenotype_counts) == 0) stop("Parent phenotype not found... Aborting...")
    
    parent_markers <- markers[parent_phen > -1]
    
  }
  
  # Split job into safe memory chunks
  
  unique_phen_size <- object.size(unique_phenotype_counts)
  
  if(!max_ram){
    max_ram <- memuse::Sys.meminfo()[[1]]
  }else{
    max_ram <- memuse::as.memuse(max_ram,"MiB")
  }
  
  chunk_size <- (n_threads*100)
  
  n_chunks <- ceiling(total_combinations/chunk_size)
  
  worst_case_mem <- unique_phen_size * 2.5 * chunk_size
  
  if(as.numeric(worst_case_mem) >= as.numeric(max_ram)){
    chunk_size <- floor(total_combinations/ceiling(n_chunks*(2*worst_case_mem/max_ram)))
    n_chunks <- ceiling(total_combinations/chunk_size)
  }
  
  if(n_chunks>total_combinations){
    n_chunks <- total_combinations
    chunk_size <- 1
    print_log("WARNING: chunk memory estimation higher than memory available!!!!")
  }
  
  print_log("Estimation of worst-case scenario: ",format(worst_case_mem, units = "auto")," of memory to use...")
  print_log("Memory available: ", capture.output(print(max_ram)))
  print_log("Spliting combinations into ",n_chunks," memory-safe chunks to be stored in ", dump_file)
  
  
  append_output <- start_from > 1
  
  
  
  last_comb <- 0
  first_comb <- start_from
  
  print_log("--------------------------------------------")
  
  no_combinations_left_counter <- 0
  
  while(TRUE){
    
    if(no_combinations_left_counter >= max_empty_chunks){
      print_log("No combinations found on last ",max_empty_chunks," chunks, terminating.")  
      break
    }
    
    if(last_comb == total_combinations) break
    
    last_comb <- first_comb+chunk_size-1
    
    if(last_comb > total_combinations){
      last_comb <- total_combinations
    }
    
    print_log("Computing marker combinations from ",format(first_comb, scientific=F), " to ",format(last_comb, scientific=F)," using ",n_threads," thread(s)...")
    
    if(max_phenotype_length > 0 & max_phenotype_length < n_markers) print_log("Filtering combinations by max phenotype length of ", max_phenotype_length)
    
    current_marker_combinations <- generate_marker_combinations(n_markers,
                                                                max_phenotype_length,
                                                                lower = first_comb,
                                                                upper = last_comb,
                                                                n_threads = n_threads)
    
    first_comb <- last_comb + 1
    
    # Skip chunk if no phenotype has less than max_phenotype_length
    if(nrow(current_marker_combinations)){
      
      no_combinations_left_counter <- 0
      
      colnames(current_marker_combinations) <- markers
      
      current_marker_combinations <- as.matrix(current_marker_combinations)
      
      # Select parent population
      if(!is.null(parent_markers) & nrow(current_marker_combinations) > 0){
        
        print_log("Filtering combinations by parent phenotype...")
        
        current_marker_combinations[, parent_markers] <- 0
        
        current_marker_combinations <- unique(current_marker_combinations)
        
      }
    
    
      
      n_combinations <- nrow(current_marker_combinations)
      
      print_log("Counting cells for ", n_combinations , " combinations...")
      
      
      
      if(n_combinations == 1){
        combinatorial_phenotypes <- group_and_reduce_phenotypes(unique_phenotype_counts,current_marker_combinations[1,])
      }else{
        combinatorial_phenotypes <- do.call(rbind,parallel::mclapply(1:n_combinations, function(j) group_and_reduce_phenotypes(unique_phenotype_counts,current_marker_combinations[j,]), mc.cores = n_threads))
      }
      
      print_log(format(nrow(combinatorial_phenotypes), scientific=F)," phenotypes generated...")
      
      if(!efficient){
        # Remove phenotypes where all samples have less then min_count cells
        if(min_count>0){
          print_log("Removing phenotypes where ",format(round(sample_fraction_min_counts*100, 2), nsmall = 2),"% of the samples have at least ",min_count," cell(s) using ",n_threads," thread(s)...")
          
          counts_list <- as.matrix(combinatorial_phenotypes[,samples_id])
          
          counts_list <- parallel::mclapply(seq_len(nrow(counts_list)), function(i) counts_list[i,], mc.cores = n_threads)
          
          cell_filter <- unlist(parallel::mclapply(counts_list, function(row) (sum(row >= min_count)/n_samples) >= sample_fraction_min_counts, mc.cores = n_threads))
          
          combinatorial_phenotypes <- combinatorial_phenotypes[cell_filter,]
          
          rm(counts_list, cell_filter)
          
          print_log(nrow(combinatorial_phenotypes)," phenotypes left...")
          
        }
      }
      
      
      print_log("Writing ",format(nrow(combinatorial_phenotypes), scientific=F)," phenotypes to file...")
      
      data.table::fwrite(data.table::as.data.table(combinatorial_phenotypes), dump_file, append = append_output, nThread = n_threads)
      
      #TODO change order of logging to print number of written phenotypes after finishing writing to avoid misslogging if process break during writing.
      
      
      # To continue appending to file
      append_output <- TRUE
      
    }else{
      print_log("No combinations left...")
      no_combinations_left_counter <- no_combinations_left_counter + 1
    }
    
    
    print_log("Marker combinations computed: ",format(last_comb, scientific=F), " out of ",format(total_combinations, scientific=F), " (",format(round((last_comb/total_combinations)*100, 2), nsmall = 2),"%)")
    
    print_log("--------------------------------------------")
    
    rm(combinatorial_phenotypes)
    gc(full = TRUE,verbose = FALSE)
    
    
  }
  
  
}



# Find last marker combination computed in log in order to continue from stalled execution
find_last_marker_combination_computed <- function(log_file){
  
  last_marker_combination <- 0
  
  con = file(log_file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(grepl("Marker combinations computed: ", line, fixed = TRUE)){
      last_marker_combination <- as.integer(stringr::str_extract(line, "(?<= Marker combinations computed: )[0-9]*"))
    }
  }
  
  close(con)
  
  return(last_marker_combination)
}



#' Count cells for each phenotype for each sample (Server Version)
#'
#' \code{combinatorial_phenotype_counts_server} generates phenotypes considering all
#' possible combination of markers and count the number of cells for each phenotype
#' for each sample.
#' 
#' This server version does all the process from raw cell data, channel data, and
#' sample data writing the results to a file. It is intended to be used with large datasets
#' which outputs wouldn't fit in the available memory. It divides all the possible
#' marker combinations into smaller chunks and compute them separatedly saving results
#' to the output file called "combinatorial_phenotype_counts.csv".
#'
#' @param cell_file Path to file containing cell data. Must be a ".csv" or ".fcs" file.
#' @param channel_file Path to a ".csv" file containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param sample_file Path to a ".csv" file containing a Sample_ID column and additional grouping columns for the samples.
#' @param output_folder Path to folder where outputs and temporary files should be saved.
#' @param parent_phen Parent phenotype to filter for. All phenotypes generated will contain the parent phenotype.
#' @param min_count Minimum number of cells that a phenotype must have for at least one sample.
#' @param sample_fraction_min_counts Fraction of samples that must have at least \code{min_count} cells. Value forced to minimum of 1/n_samples.
#' @param max_phenotype_length Maximum length of markers to compose a phenotype.
#' @param sample_ID_col Name of the column in \code{cell_data} where the Sample IDs are stored. Default: "Sample_ID".
#' @param save_cell_data If TRUE, processed cell data is saved to "output_folder/cell_data.csv".
#' @param continue If TRUE, look for files to resume execution. Also needed to save necessary "continuing" files.
#' @param max_ram Maximum ram memory in Gb to be used by the function.
#' @param efficient If TRUE, filter full-length phenotypes for \code{min_count} condition before generating all marker combinations.
#' It is less sensitive for very rare phenotypes but yields a great boost in performance.
#' @param max_empty_chunks Maximum number of chunks with no valid combinations when filtering by \code{max_phenotype_length}.
#' @param n_threads Number of threads to be used. Default: 1.
#' @param verbose If TRUE, print outputs from log to stdout.
#'  
#' @return Output is saved to file.
#' 
#' @export
combinatorial_phenotype_counts_server <- function(cell_file,
                                                  channel_file,
                                                  sample_file,
                                                  output_folder,
                                                  parent_phen = NULL,
                                                  min_count = 10,
                                                  sample_fraction_min_counts = 0.0,
                                                  max_phenotype_length = 0,
                                                  sampleID_col = "Sample_ID",
                                                  save_cell_data = TRUE,
                                                  continue = TRUE,
                                                  max_ram = 0,
                                                  efficient = TRUE,
                                                  max_empty_chunks = 1000,
                                                  n_threads = 1,
                                                  verbose = TRUE
){
  
  continued <- FALSE
  log_file <- "combinatorial_phenotypes.log"
  unique_phen_file <- "unique_phen.tmp"
  cell_data_file <- "processed_cell_data.csv"
  phenotype_counts_file <- "combinatorial_phenotype_counts.csv"
  
  options(PhenoComb.verbose = verbose)
  
  
  if(continue){
    
    if(dir.exists(output_folder)){
      
      
      # Get files in output folder
      folder_files <- list.files(output_folder, full.names = F, recursive = F)
      
      if(length(folder_files)>0){# If there are files in the output_folder
        
        # Check if all continuation files exist
        log_file_exists <- log_file %in% folder_files
        unique_phen_file_exists <- unique_phen_file %in% folder_files
        phenotype_counts_file_exists <- phenotype_counts_file %in% folder_files
        
        
        if(all(c(log_file_exists,
                 unique_phen_file_exists,
                 phenotype_counts_file_exists)))
        { # If all continuation files exist
          
          # Check for chunk information on log file
          last_marker_combination <- find_last_marker_combination_computed(file.path(output_folder,log_file))
          
          if(last_marker_combination){ #If chunk information found
            
            continued <- TRUE
            
            start_log(log_path = file.path(output_folder,log_file), append = TRUE, verbose = verbose)
            
            print_log("Continuation files found. Continuing...")
            
            print_log("Reading unique phenotypes counts from ",file.path(output_folder,unique_phen_file))
            unique_phen_data <- readRDS(file.path(output_folder,unique_phen_file))
            print_log("Done.")
            
            memory_safe_combinatorial_phenotype_counts( unique_phen_data$unique_phenotypes,
                                                        unique_phen_data$n_markers,
                                                        file.path(output_folder,phenotype_counts_file),
                                                        parent_phen = parent_phen,
                                                        max_phenotype_length = max_phenotype_length,
                                                        min_count = min_count,
                                                        sample_fraction_min_counts = sample_fraction_min_counts,
                                                        start_from = last_marker_combination+1,
                                                        max_ram = max_ram,
                                                        efficient = efficient,
                                                        max_empty_chunks = max_empty_chunks,
                                                        n_threads = n_threads
            )
            
            print_log("Combinatorial phenotype cell counting done.")
            
            print_log("Removing temporary files...")
            file.remove(file.path(output_folder,unique_phen_file))
            print_log("Done.")
            
            stop_log() 
            
          }else{
            
            print_log("No chunk information found. Starting from scratch...")
            
          }
          
          
          
          # Do continuing stuff
          
        }else{ # If log file does not exist
          print_log("Continuation files not found. Starting from scratch...")
        }
        
      }else{# if no files in output_folder
        print_log("Continuation files not found. Starting from scratch...")
      }
      
    }else{ # if output_dir does not exist
      print_log("Folder ", output_folder, " does not exist. Starting from scratch...")
    }
    
  }
  
  
  if(!continued){
    
    
    
    if(!dir.exists(output_folder)){
      
      # Create output folder
      print_log("Creating output folder at ", output_folder)
      dir.create(output_folder)
      
    }
    
    start_log(file.path(output_folder,log_file), verbose = verbose, append = FALSE)
    
    # Read fcs file
    print_log("Reading cell data from ", cell_file)
    file_extension <- tools::file_ext(cell_file)
    if(!(file_extension %in% c("csv","fcs"))) stop("Cell data file format not supported. Must be .csv or .fcs.")
    if(file_extension == "csv"){
      cell_data <- as.data.frame(data.table::fread(cell_file,check.names = FALSE))
    }else if(file_extension == "fcs"){
      cell_data <- as.data.frame(flowCore::read.FCS(cell_file,truncate_max_range = FALSE)@exprs)  
    }
    
    # Read channel data
    print_log("Reading Channel data from ", channel_file)
    channel_data <- as.data.frame(data.table::fread(channel_file,check.names = FALSE))
    
    # Read sample data
    print_log("Reading Sample data from ", sample_file)
    sample_data <- as.data.frame(data.table::fread(sample_file,check.names = FALSE))
    
    
    print_log("Processing FCS data...")
    
    processed_cell_data <- process_cell_data(cell_data, channel_data, sample_data, sampleID_col = sampleID_col, n_threads = n_threads)
    
    rm(cell_data)
    
    if(save_cell_data){
      print_log("Saving cell data to ",file.path(output_folder,cell_data_file))
      data.table::fwrite(processed_cell_data,file.path(output_folder,cell_data_file))
      print_log("Writing done.")
    }
    
    # Get unique phenotypes for each sample and their counts
    print_log("Getting unique phenotypes from ",nrow(processed_cell_data)," cells...")
    
    unique_phen <- get_unique_phenotype_counts(processed_cell_data, min_count, sample_fraction_min_counts, efficient, n_threads)
    
    
    
    # Saving unique phenotypes for eventual continuation
    if(continue){
      print_log("Saving ",nrow(unique_phen)," unique phenotype counts to ", file.path(output_folder,unique_phen_file))
      saveRDS(list(n_markers = nrow(channel_data), unique_phenotypes =  unique_phen), file.path(output_folder,unique_phen_file))
      print_log("Saving done.")  
    }
    
    rm(processed_cell_data)
    gc(full = TRUE,verbose = FALSE)
    
    # Create phenotype counting file
    file.create(file.path(output_folder,phenotype_counts_file), overrite = TRUE)
    
    print_log("Starting cell counting for all combinatorial phenotypes...")
    memory_safe_combinatorial_phenotype_counts(unique_phen,
                                               nrow(channel_data),
                                               file.path(output_folder,phenotype_counts_file),
                                               parent_phen = parent_phen,
                                               max_phenotype_length = max_phenotype_length,
                                               min_count = min_count,
                                               sample_fraction_min_counts = sample_fraction_min_counts,
                                               max_ram = max_ram,
                                               efficient = efficient,
                                               max_empty_chunks = max_empty_chunks,
                                               n_threads = n_threads
    )
    print_log("Combinatorial phenotype cell counting done.")
    
    if(continue){
      print_log("Removing temporary files...")
      file.remove(file.path(output_folder,unique_phen_file))
    }
    print_log("Done.")
    
    stop_log()  
  }
  
}





