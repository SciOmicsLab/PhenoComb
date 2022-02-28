


#' Encode phenotype name
#'
#' \code{phenotype_to_numbers} generates a numeric representation of a human-readable phenotype name
#' using the following encoding:
#' -1 : Droped marker
#' 0  : - (negative)
#' 1  : + (positive)
#' 2  : ++ (double positive)
#' More positive states can also be represented (e.g. 3, 4... : +++, ++++...)
#'
#' @param phenotype Phenotype name.
#' @param markers Names of the markers.
#' @return Numeric vector with encoded phenotype.
#'
#' @examples
#' phenotype_to_numbers("Marker1-Marker2+Marker3++Marker4-",c("Marker1","Marker2","Marker3","Marker4"))
#' phenotype_to_numbers("Marker1-Marker2+",c("Marker1","Marker2","Marker3"))
#' 
#' @export
phenotype_to_numbers <- function(phenotype, markers){
  
  unlist(lapply(markers, function(m){
    
    has_marker <- stringr::str_extract(phenotype, paste(m,"[\\+\\-]+", sep = ""))
    if( is.na(has_marker) ) return(-1)
    if( lengths(regmatches(has_marker, gregexpr("\\-", has_marker))) ) return(0)
    return(lengths(regmatches(has_marker, gregexpr("\\+", has_marker))))
    
  }))
  
}


# Compares two phenotypes (numeric vectors)
phenotypes_are_equal <- function(phen1, phen2){
  
  return(all(phen1 == phen2))
  
}

# Checks if phen1 contains phen2
has_phenotype <- function(phen1, phen2){
  
  has_marker <- phen2 > -1
  
  return(all(phen1[has_marker] == phen2[has_marker]))
  
}



#' Find phenotype in data.frame
#' 
#' @export
find_phenotype <- function(phen_table, phen, markers, n_cores = 1){
  
  # phen has to be numeric vector form
  if(is.character(phen)){
    phen <- phenotype_to_numbers(phen,markers)  
  }
  
  comparison <- unlist(parallel::mclapply(1:nrow(phen_table), function(i) phenotypes_are_equal(phen_table[i,markers],phen), mc.cores = n_cores))
  
  return(phen_table[comparison,])
  
}

#' Find phenotype in csv file
#' 
#' @export
find_phenotype_in_file <- function(phen_file, phen, markers, n_cores = 1, chunk_size = 10000){
  
  # phen must be numeric encoded
  if(is.character(phen)){
    phen <- phenotype_to_numbers(phen,markers)  
  }
  
  #Get header and col types
  file_sample <- data.table::fread(phen_file, nrows = 2, header= T)
  header <- colnames(file_sample)
  
  column_types <- as.vector(sapply(file_sample, typeof))
  column_types[column_types=="character"] <- "string"
  
  # Use LaF to read big files by chunks
  laf <- LaF::laf_open_csv(phen_file, column_types = column_types, skip = 1)
  
  while(TRUE){
    
    phenotypes <- LaF::next_block(laf, nrows=chunk_size)
    if (nrow(phenotypes) == 0) break;
    colnames(phenotypes) <- header
    
    found_phen <- find_phenotype(phenotypes, phen, markers, n_cores)
    
    if(nrow(found_phen) > 0) return(found_phen)
    
  }
  
  return(NULL)
  
  
}

count_csv_lines <- function(file_path){
  
  extension <- tools::file_ext(file_path)
  
  if(extension == "csv"){
    return(wc_count_lines(file_path) - 1)
  }else if(extension == "log"){
    return(get_n_phenotypes_log(file_path))
  }else{
    stop("File provided must be .csv or .log")
  }
  
}


get_n_phenotypes_log <- function(file_path){
  
  txt <- readLines(file_path)
  
  phenotypes_generated <- 0
  
  for(i in 1:length(txt)){
    
    if(grepl( "Writing", txt[i], fixed = TRUE)){
      
      phenotypes_generated <- phenotypes_generated + as.integer(sub(".*?Writing .*?(\\d+).*", "\\1", txt[i]))
      
    }
    
  }
  
  return(phenotypes_generated)
  
}

# Fast line counting
wc_count_lines <- function(file_path){
  as.numeric(system(paste("cat ",file_path," | wc -l", collapse = ""), intern = TRUE))
}

# Normalize numeric vector
normalize_vals <- function(x){
  x <- abs(x)
  x <- x-min(x)
  x <- x/max(x)
  return(x)
}