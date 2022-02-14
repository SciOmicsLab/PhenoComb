



check_sample_data <- function(sample_data){
  
  if(!("Sample_ID" %in% colnames(sample_data))){
    
    stop("'Sample_ID' column not found in sample data. IDs must be listed in the Sample_ID column. ")
    
  }
  
}


check_cell_data <- function(cell_data, sampleID_col){
  
  if(!any(colnames(cell_data) == sampleID_col)){
    
    stop(paste("Column", sampleID_col," not found in cell data.", sep = ""))
    
  }
  
}

check_channel_data <- function(channel_data){
  
  if(!("Channel" %in% colnames(channel_data))){
    
    stop("'Channel' column not found in channel data. The 'Channel' column matching exactly the channel names in cell data is required.")
    
  }
  
  if(!("Marker" %in% colnames(channel_data))){
    
    stop("'Marker' column not found in channel data. The 'Marker' column with marker names is required.")
    
  }
  
  if(!("T1" %in% colnames(channel_data))){
    
    stop("'T1' column not found in channel data. At least one threshold column, named 'T1', must be provided. More thresholds can be indicated as 'T2', 'T3'...")
    
  }
  
}


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
  
  # Check data format required
  check_channel_data(channel_data)
  check_sample_data(sample_data)
  check_cell_data(cell_data, sampleID_col)
  
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
    
    oob_data <- as.numeric(channel_data[,"OOB"])
    
    cell_data_list <- as.matrix(cell_data[,channel_data[,"Marker"]])
    
    cell_data_list <- parallel::mclapply(seq_len(nrow(cell_data_list)), function(i) cell_data_list[i,], mc.cores = n_threads)
    
    oob_filter <- unlist(parallel::mclapply(cell_data_list, function(row) !any(row > oob_data), mc.cores = n_threads))
    
    cell_data <- cell_data[oob_filter,]
    
    rm(cell_data_list, oob_filter)
    gc(full = TRUE,verbose = FALSE)
    
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
get_unique_phenotype_counts <- function(processed_cell_data, min_count = 0, sample_fraction_min_counts = 0.0, efficient = TRUE, n_threads = 1){
  
  unique_phen <- as.data.frame(dplyr::ungroup(dplyr::summarise(dplyr::group_by_all(processed_cell_data), Count = dplyr::n())))
  
  sample_ids <- unique(unique_phen[,"Sample_ID"])
  markers <- head(colnames(unique_phen),length(unique_phen)-2)
  
  n_samples <- length(sample_ids)
  
  if(sample_fraction_min_counts <= 1/n_samples){
    sample_fraction_min_counts <- 1/n_samples
  }
  
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
      
      print_log("Removing phenotypes where ",format(round(sample_fraction_min_counts*100, 2), nsmall = 2),"% of the samples have at least ",min_count," cell(s) using ",n_threads," thread(s)...")
      
      counts_list <- as.matrix(unique_phen[,sample_ids])
      
      counts_list <- parallel::mclapply(seq_len(nrow(counts_list)), function(i) counts_list[i,], mc.cores = n_threads)
      
      cell_filter <- unlist(parallel::mclapply(counts_list, function(row) (sum(row >= min_count)/n_samples) >= sample_fraction_min_counts, mc.cores = n_threads))
      
      unique_phen <- unique_phen[cell_filter,]
      
      rm(counts_list, cell_filter)
      gc(full = TRUE,verbose = FALSE)
      
      print_log(nrow(unique_phen)," unique phenotypes left...")
      
    }
  }
  
  return(as.matrix(unique_phen))
}



# Generates combinatorial matrix with markers to drop
generate_marker_combinations <- function(n_markers, max_phenotype_length = 0, lower = NULL, upper = NULL, n_threads = 1){
  
  local_marker_comb <- as.data.frame(RcppAlgos::permuteGeneral(c(1,0), n_markers, repetition = TRUE,
                                      lower = lower,
                                      upper = upper,
                                      ))
  
  # KNOWN BUG
  # if last command returns only one combination
  # will probably produce a "argument of length 0" error
  
  if(max_phenotype_length > 0 & max_phenotype_length < n_markers){
    local_marker_comb <- local_marker_comb[unlist(parallel::mclapply(1:nrow(local_marker_comb),
                                                                     function(i) (n_markers-sum(local_marker_comb[i,])) <= max_phenotype_length,
                                                                     mc.cores = n_threads)),]  
  }
  
  return(local_marker_comb)
  
}



#' Count cells for each phenotype for each sample
#'
#' \code{combinatorial_phenotype_counts} generates phenotypes considering all
#' possible combination of markers and count the number of cells for each phenotype
#' for each sample.
#'
#' @param processed_cell_data Data.Frame containing filtered and thresholded cell data.
#' 
#' Use function \code{process_cell_data} to generate this input.
#' @param parent_phen Parent phenotype to filter for. All phenotypes generated will contain the parent phenotype.
#' @param min_count Minimum number of cells that a phenotype must have for at least one sample.
#' @param sample_fraction_min_counts Fraction of samples that must have at least \code{min_count} cells (value between 0 ~ 1). Value forced to minimum of 1/n_samples.
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
                                           sample_fraction_min_counts = 0.0,
                                           max_phenotype_length = 0,
                                           efficient = FALSE,
                                           n_threads = 1
                                           ){
  
  
  # Get markers
  markers <- head(colnames(processed_cell_data),length(processed_cell_data)-1)
  
  # Number of markers
  n_markers <- length(markers)
  
  # Get Sample IDs
  sample_ids <- unique(processed_cell_data[,"Sample_ID"])
  
  
  
  # Number of samples
  n_samples <- length(sample_ids)
  
  # Get unique phenotypes for each sample and their counts
  print_log("Getting unique phenotypes from ",nrow(processed_cell_data)," cells...")
  unique_phen <- get_unique_phenotype_counts(processed_cell_data, min_count, sample_fraction_min_counts, efficient, n_threads)
  
  # Create all combination of markers
  
  print_log("Generating all ", format(2^n_markers, scientific = FALSE), " marker combinations...")
  
  marker_combinations <- generate_marker_combinations(n_markers, max_phenotype_length, n_threads = n_threads)
  colnames(marker_combinations) <- markers
  
  # Select parent population
  if(!is.null(parent_phen)){
    
    print_log("Filtering for parent population ", parent_phen)
    
    parent_phen <- phenotype_to_numbers(parent_phen, markers)
    
    has_parent_phen <- unlist(parallel::mclapply(1:nrow(unique_phen), function(i) has_phenotype(unique_phen[i,markers], parent_phen), mc.cores = n_threads))
    
    unique_phen <- unique_phen[has_parent_phen,]
    
    parent_markers <- markers[parent_phen > -1]
    
    marker_combinations[, parent_markers] <- FALSE
    
    marker_combinations <- unique(marker_combinations)
    
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
      
      if(sample_fraction_min_counts <= 1/n_samples){
        sample_fraction_min_counts <- 1/n_samples
      }
      
      counts_list <- as.matrix(combinatorial_phenotypes[,samples_id])
      
      counts_list <- parallel::mclapply(seq_len(nrow(counts_list)), function(i) counts_list[i,], mc.cores = n_threads)
      
      cell_filter <- unlist(parallel::mclapply(counts_list, function(row) (sum(row >= min_count)/n_samples) >= sample_fraction_min_counts, mc.cores = n_threads))
      
      combinatorial_phenotypes <- combinatorial_phenotypes[cell_filter,]
      
      rm(counts_list, cell_filter)
      
      combinatorial_phenotypes <- combinatorial_phenotypes[!unlist(parallel::mclapply(1:nrow(combinatorial_phenotypes), function(i) all(combinatorial_phenotypes[i,samples_id] < min_count), mc.cores = n_threads)),]
      
      print_log(nrow(combinatorial_phenotypes)," phenotypes left...")
    }  
    
  }
  
  return(as.data.frame(combinatorial_phenotypes))
  
}



