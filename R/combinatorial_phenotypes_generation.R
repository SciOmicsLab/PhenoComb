

#' Process raw cell data
#'
#' \code{process_cell_data} uses the channel_data and sample_data to filter
#' columns and cells from cell_data and discretize values based on thresholds
#' defined on channel_data.
#'
#' @param cell_data Data.Frame where each column is a channel of measurement plus one column with sample IDs. Each row is a cell.
#' @param channel_data Data.Frame containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param sample_data Data.Frame containing a Sample_ID column and additional grouping columns for the samples.
#' @param sample_ID_col Name of the column in \code{cell_data} where the Sample IDs are stored. Default: "Sample_ID"
#' @param n_threads Number of threads to be used (Only useful if there is an OOB column in channel_data). Default: 1
#'  
#' @return Data.Frame with filtered and discretized cells.
#' 
#' @export
process_cell_data <- function(cell_data, channel_data, sample_data, sampleID_col = "Sample_ID", n_threads = 1){
  
  # Must be dataframes
  cell_data <- as.data.frame(cell_data)
  channel_data <- as.data.frame(channel_data)
  sample_data <- as.data.frame(sample_data)
  
  if(!any(colnames(cell_data) == sampleID_col)){
    stop(paste("Column", sampleID_col," not found in cell data.", sep = ""))
  }
  
  # Rename SampleID and filter samples if in sample data
  colnames(cell_data)[colnames(cell_data) == sampleID_col] <- 'Sample_ID'
  cell_data[,"Sample_ID"] <- as.character(cell_data[,"Sample_ID"])
  sample_data[,"Sample_ID"] <- as.character(sample_data[,"Sample_ID"])
  
  
  # Select channel columns and rename to markers
  print_log("Filtering cell data based on channel data file...")
  cell_data = cell_data[,c(channel_data[,"Channel"],"Sample_ID")]
  colnames(cell_data) <- c(channel_data[,"Marker"],"Sample_ID")
  cell_data <- cell_data[cell_data[,"Sample_ID"] %in% sample_data[,"Sample_ID"],]
  
  # Set NA values to Inf
  channel_data[is.na(channel_data)] <- Inf
  
  
  # Remove events if any value is higher then OBB if OBB columns is present
  have_obb_col <- "OOB" %in% colnames(channel_data)
  if(have_obb_col){
    
    print_log("Removing cells with values over OOB using ",n_threads," thread(s)...")
    
    cell_data <- cell_data[unlist(parallel::mclapply(1:nrow(cell_data),function(i) !any(cell_data[i,channel_data[,"Marker"]] > channel_data[,"OOB"]), mc.cores = n_threads)),]
    
  }
  
  # Apply thresholds to discretize data
  
  markers <- channel_data[,"Marker"]
  n_markers <- nrow(channel_data)
  threshold_cols <- colnames( channel_data )[grepl("T", colnames( channel_data ) )]
  for(m in 1:n_markers){
    thresholds <- as.numeric(channel_data[m,threshold_cols])
    thresholds <- thresholds[!is.infinite(thresholds)]
    
    print_log("Applying ",length(thresholds)," threshold(s) to marker ",markers[m])
    
    thresholds <- c(-Inf,thresholds,Inf)
    cell_data[,markers[m]] <- as.numeric(cut(cell_data[,markers[m]],breaks = thresholds))-1
  }
  
  
  # Return the processed dataframe
  return(cell_data)
  
}


# Return unique phenotypes with cell count by Sample_ID
get_unique_phenotype_counts <- function(processed_cell_data, min_count = 0, efficient = TRUE, n_threads = 1){
  
  unique_phen <- as.data.frame(dplyr::ungroup(dplyr::summarise(dplyr::group_by_all(processed_cell_data), Count = dplyr::n())))
  
  sample_ids <- unique(unique_phen[,"Sample_ID"])
  markers <- head(colnames(unique_phen),length(unique_phen)-2)
  

  unique_phen <- parallel::mclapply(sample_ids, function(id) {
      sample_phenotypes <- unique_phen[unique_phen[,"Sample_ID"] == id, c(markers,"Count")]
      colnames(sample_phenotypes) <- c(markers,id)
      return(sample_phenotypes)
    }, mc.cores = n_threads)
  
  # Join all sample ID dataframes by Phenotypes
  unique_phen <- Reduce(function(dtf1, dtf2) dplyr::full_join(dtf1, dtf2, by = markers), unique_phen)
  
  #Set NA to zero
  unique_phen[is.na(unique_phen)] <- 0
  
  print_log(nrow(unique_phen)," unique phenotypes generated...")
  
  if(efficient){
    
    # Remove phenotypes where all samples have less then min_count cells
    if(min_count > 0){
      
      print_log("Removing phenotypes where all samples have less than ",min_count," cell(s) using ",n_threads," thread(s)...")
      unique_phen <- unique_phen[unlist(parallel::mclapply(1:nrow(unique_phen), function(j) !all(unique_phen[j,sample_ids] < min_count), mc.cores = n_threads)),]
      print_log(nrow(unique_phen)," unique phenotypes left...")
      
    }
  }
  
  return(as.matrix(unique_phen))
}






#' Count cells for each phenotype for each sample
#'
#' \code{combinatorial_phenotype_counts} generates phenotypes considering all
#' possible combination of markers and count the number of cells for each phenotype
#' for each sample.
#'
#' @param processed_cell_data Data.Frame containing filtered and thresholded cell data.
#' Use function \code{process_cell_data} to generate this input.
#' @param parent_phen Parent phenotype to filter for. All phenotypes generated will contain the parent phenotype.
#' @param min_count Minimum number of cells that a phenotype must have for at least one sample.
#' @param max_phenotype_length Maximum length of markers to compose a phenotype.
#' @param efficient If TRUE, filter full-length phenotypes for \code{min_count} condition before generating all marker combinations.
#' It is less sensitive for very rare phenotypes but yields a great boost in performance.
#' @param n_threads Number of threads to be used. Default: 1.
#'  
#' @return Data.Frame with all possible phenotypes and cell counts for each sample.
#' 
#' @export
combinatorial_phenotype_counts <- function(processed_cell_data,
                                           parent_phen = NULL,
                                           min_count = 10,
                                           max_phenotype_length = 0,
                                           efficient = FALSE,
                                           n_threads = 1
                                           ){
  
  
  # Get markers
  markers <- head(colnames(processed_cell_data),length(processed_cell_data)-1)
  
  # Number of markers
  n_markers <- length(markers)
  
  
  # Get unique phenotypes for each sample and their counts
  print_log("Getting unique phenotypes from ",nrow(processed_cell_data)," cells...")
  unique_phen <- get_unique_phenotype_counts(processed_cell_data, min_count, efficient, n_threads)
  
  # Create all combination of markers
  
  print_log("Generating all ", format(2^n_markers, scientific = FALSE), " marker combinations...")
  
  marker_combinations <- expand.grid(rep(list(c(T,F)), n_markers))
  colnames(marker_combinations) <- markers
  
  if(max_phenotype_length > 0 & max_phenotype_length < n_markers){
    
    print_log("Filtering combinations with maximum number of markers of ", max_phenotype_length, ".")
    marker_combinations <- marker_combinations[unlist(parallel::mclapply(1:nrow(marker_combinations), function(i) (n_markers-sum(marker_combinations[i,])) <= max_phenotype_length, mc.cores = n_threads)),]

  }
  
  # Select parent population
  if(!is.null(parent_phen)){
    
    print_log("Filtering for parent population ", parent_phen)
    
    parent_phen <- phenotype_to_numbers(parent_phen, markers)
    
    has_parent_phen <- unlist(parallel::mclapply(1:nrow(unique_phen), function(i) has_phenotype(unique_phen[i,markers], parent_phen), mc.cores = n_threads))
    
    unique_phen <- unique_phen[has_parent_phen,]
    
    parent_markers <- markers[parent_phen > -1]
    
    marker_combinations[, parent_markers] <- FALSE
    
    marker_combinations <- dplyr::distinct(marker_combinations)
    
  }
  
  marker_combinations <- as.matrix(marker_combinations)
  
  
  print_log("Counting cells for each of the ",nrow(marker_combinations)," marker combinations using ",n_threads," thread(s)...")
  # Compute counts for all phenotype combinations
  combinatorial_phenotypes <- do.call(rbind,parallel::mclapply(1:nrow(marker_combinations), function(i) group_and_reduce_phenotypes(unique_phen, marker_combinations[i,]),mc.cores = n_threads))

  print_log(nrow(combinatorial_phenotypes)," unique phenotypes generated...")

  samples_id <- tail(colnames(combinatorial_phenotypes), length(combinatorial_phenotypes)-length(markers))
  
  if(!efficient){
    # Remove phenotypes where all samples have less then min_count cells
    if(min_count>0){
      print_log("Removing phenotypes where all samples have less than ",min_count," cell(s) using ",n_threads," thread(s)...")
      
      combinatorial_phenotypes <- combinatorial_phenotypes[!unlist(mclapply(1:nrow(combinatorial_phenotypes), function(i) all(combinatorial_phenotypes[i,samples_id] < min_count), mc.cores = n_threads)),]
      
      print_log(nrow(combinatorial_phenotypes)," phenotypes left...")
    }  
    
  }
  
  return(as.data.frame(combinatorial_phenotypes))
  
}




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
                                                       start_from = 0,
                                                       max_ram = 0,
                                                       efficient = TRUE,
                                                       n_threads = 1
                                                      ){
  
  
  
  # Get markers names and sample IDs
  markers <- head(colnames(unique_phenotype_counts),n_markers)
  samples_id <- tail(colnames(unique_phenotype_counts),ncol(unique_phenotype_counts)-n_markers)
  
  
  # Create all combination of markers
  
  print_log("Generating all ",format(2^n_markers, scientific = FALSE)," marker combinations.")
  
  marker_combinations <- expand.grid(rep(list(c(T,F)),n_markers))
  colnames(marker_combinations) <- markers
  
  
  # Filtering combinations based on the maximum number of markers considered
  if(max_phenotype_length > 0 & max_phenotype_length < n_markers){
    print_log("Filtering combinations with maximum number of markers of ",max_phenotype_length,".")
    
    marker_combinations <- marker_combinations[unlist(parallel::mclapply(1:nrow(marker_combinations), function(i) (n_markers-sum(marker_combinations[i,])) <= max_phenotype_length, mc.cores = n_threads)),]
    
  }
  
  # Select parent population
  if(!is.null(parent_phen)){
    
    print_log("Filtering for parent population ", parent_phen)
    
    parent_phen <- phenotype_to_numbers(parent_phen, markers)
    
    has_parent_phen <- unlist(parallel::mclapply(1:nrow(unique_phen), function(i) has_phenotype(unique_phen[i,markers], parent_phen), mc.cores = n_threads))
    
    unique_phen <- unique_phen[has_parent_phen,]
    
    parent_markers <- markers[parent_phen > -1]
    
    marker_combinations[, parent_markers] <- FALSE
    
    marker_combinations <- dplyr::distinct(marker_combinations)
    
  }
  
  marker_combinations <- as.matrix(marker_combinations)
  n_marker_combinations <- nrow(marker_combinations)
  rownames(marker_combinations) <- 1:n_marker_combinations
  
  if(start_from > 1){
    marker_combinations <- marker_combinations[start_from:nrow(marker_combinations), ]
    print_log("Continuing to compute the ", n_marker_combinations, " marker combinations from combination #", start_from," to be stored in ", dump_file)
    print_log("Still ", nrow(marker_combinations), " to go...")
  }
  
  current_n_marker_combinations <- nrow(marker_combinations)
  
  # Split job into safe memory chunks
  
  unique_phen_size <- object.size(unique_phenotype_counts)
  
  if(!max_ram){
    max_ram <- memuse::Sys.meminfo()[[1]]
  }else{
    max_ram <- memuse::as.memuse(max_ram,"MiB")
  }
  
  chunk_size <- (n_threads*100)
  
  n_chunks <- ceiling(current_n_marker_combinations/chunk_size)
  
  worst_case_mem <- unique_phen_size * 2.5 * chunk_size
  
  if(as.numeric(worst_case_mem) >= as.numeric(max_ram)){
    chunk_size <- floor(current_n_marker_combinations/ceiling(n_chunks*(2*worst_case_mem/max_ram)))
    n_chunks <- ceiling(current_n_marker_combinations/chunk_size)
  }
  
  if(n_chunks>current_n_marker_combinations){
    n_chunks <- current_n_marker_combinations
    chunk_size <- 1
    print_log("WARNING: chunk memory estimation higher than memory available!!!!")
  }
  
  print_log("Estimation of worst-case scenario: ",format(worst_case_mem, units = "auto")," of memory to use...")
  print_log("Memory available: ", capture.output(print(max_ram)))
  print_log("Spliting combinations into ",n_chunks," memory-safe chunks to be stored in ", dump_file)
  
  
  split_marker_combinations <- split_in_chunks(current_n_marker_combinations, chunk_size)
  
  combinations_id <- as.numeric(rownames(marker_combinations))
  
  append_output <- start_from > 1
  
  print_log("--------------------------------------------")
  
  
  for(current_combinations in split_marker_combinations){
    
    n_combinations <- length(current_combinations)
    
    
    print_log("Computing marker combinations from ",format(combinations_id[current_combinations[1]], scientific=F), " to ",format(combinations_id[current_combinations[n_combinations]], scientific=F))
    
    print_log("Counting cells for each of the ",format(n_combinations, scientific=F)," marker combinations using ",n_threads," thread(s)...")
    
    
    if(n_combinations == 1){
      combinatorial_phenotypes <- group_and_reduce_phenotypes(unique_phenotype_counts,marker_combinations[current_combinations[1]])
    }else{
      combinatorial_phenotypes <- do.call(rbind,parallel::mclapply(current_combinations, function(j) group_and_reduce_phenotypes(unique_phenotype_counts,marker_combinations[j,]), mc.cores = n_threads))
    }
    
    print_log(format(nrow(combinatorial_phenotypes), scientific=F)," phenotypes generated...")
    
    if(!efficient){
      # Remove phenotypes where all samples have less then min_count cells
      if(min_count>0){
        print_log("Removing phenotypes where all samples have less than ",min_count," cell(s) using ",n_threads," thread(s)...")
        combinatorial_phenotypes <- combinatorial_phenotypes[unlist(parallel::mclapply(1:nrow(combinatorial_phenotypes), function(j) !all(combinatorial_phenotypes[j,samples_id] < min_count), mc.cores = n_threads)),]
        print_log(nrow(combinatorial_phenotypes)," phenotypes left...")
        
      }
    }
    
    
    print_log("Writing ",format(nrow(combinatorial_phenotypes), scientific=F)," phenotypes to file...")
    
    data.table::fwrite(data.table::as.data.table(combinatorial_phenotypes), dump_file, append = append_output, nThread = n_threads)
    
    # To continue appending to file
    append_output <- TRUE
    
    print_log("Marker combinations computed: ",format(combinations_id[tail(current_combinations,1)], scientific=F), " out of ",format(n_marker_combinations, scientific=F), " (",format(round((combinations_id[tail(current_combinations,1)]/n_marker_combinations)*100, 2), nsmall = 2),"%)")
    
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
#' @param max_phenotype_length Maximum length of markers to compose a phenotype.
#' @param sample_ID_col Name of the column in \code{cell_data} where the Sample IDs are stored. Default: "Sample_ID".
#' @param save_cell_data If TRUE, processed cell data is saved to "output_folder/cell_data.csv".
#' @param continue If TRUE, look for files to resume execution. Also needed to save necessary "continuing" files.
#' @param max_ram Maximum ram memory in Gb to be used by the function.
#' @param efficient If TRUE, filter full-length phenotypes for \code{min_count} condition before generating all marker combinations.
#' It is less sensitive for very rare phenotypes but yields a great boost in performance.
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
                                                  max_phenotype_length = 0,
                                                  sampleID_col = "Sample_ID",
                                                  save_cell_data = TRUE,
                                                  continue = TRUE,
                                                  max_ram = 0,
                                                  efficient = TRUE,
                                                  n_threads = 1,
                                                  verbose = TRUE
                                                  ){
  
  continued <- FALSE
  log_file <- "combinatorial_phenotypes_log.txt"
  unique_phen_file <- "unique_phen.tmp"
  cell_data_file <- "processed_cell_data.csv"
  phenotype_counts_file <- "combinatorial_phenotype_counts.csv"
  
  
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
                                                        start_from = last_marker_combination+1,
                                                        max_ram = max_ram,
                                                        efficient = efficient,
                                                        n_cores = n_cores
                                                        )
            
            print_log("Combinatorial phenotype cell counting done.")
            
            print_log("Removing temporary files...")
            file.remove(file.path(output_folder,unique_phen_file))
            print_log("Done.")
            
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
      message("Folder ", output_folder, " does not exist. Starting from scratch...")
    }
    
  }
  
  
  if(!continued){
    
    
    
    if(!dir.exists(output_folder)){
      
      # Create output folder
      message("Creating output folder at ", output_folder)
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
    
    cell_data <- process_cell_data(cell_data, channel_data, sample_data, sampleID_col = sampleID_col, n_threads = n_threads)
    
    if(save_cell_data){
      print_log("Saving cell data to ",file.path(output_folder,cell_data_file))
      data.table::fwrite(cell_data,file.path(output_folder,cell_data_file))
      print_log("Writing done.")
    }
    
    # Get unique phenotypes for each sample and their counts
    print_log("Getting unique phenotypes from ",nrow(cell_data)," cells...")
    
    unique_phen <- get_unique_phenotype_counts(cell_data, min_count, efficient, n_threads)
    
    rm(cell_data)
    
    # Saving unique phenotypes for eventual continuation
    if(continue){
      print_log("Saving ",nrow(unique_phen)," unique phenotype counts to ", file.path(output_folder,unique_phen_file))
      saveRDS(list(n_markers = nrow(channel_data), unique_phenotypes =  unique_phen), file.path(output_folder,unique_phen_file))
      print_log("Saving done.")  
    }
    
    # Create phenotype counting file
    file.create(file.path(output_folder,phenotype_counts_file), overrite = TRUE)
    
    print_log("Starting cell counting for all combinatorial phenotypes...")
    memory_safe_combinatorial_phenotype_counts(unique_phen,
                                               nrow(channel_data),
                                               file.path(output_folder,phenotype_counts_file),
                                               parent_phen = parent_phen,
                                               max_phenotype_length = max_phenotype_length,
                                               min_count = min_count,
                                               max_ram = max_ram,
                                               efficient = efficient,
                                               n_threads = n_threads
                                               )
    print_log("Combinatorial phenotype cell counting done.")
    
    print_log("Removing temporary files...")
    file.remove(file.path(output_folder,unique_phen_file))
    print_log("Done.")
    
    stop_log()  
  }
  
}


