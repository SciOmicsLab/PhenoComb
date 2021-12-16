

#' Identify independent relevant phenotypes (Server Version)
#' 
#' It uses a network clustering process to identify most representative
#' phenotypes when several similar phenotypes are relevant.
#' 
#' This server version is just a convenience to find correct inputs and
#' outputs results to file. No performance difference.
#' 
#' Parameters \code{groups_column, g1, g2, max_pval,} and \code{parent_phen} should
#' match the parameter used in \code{statistically_relevant_phenotypes_server} in order
#' to find correct input files.
#' 
#' @param output_folder Path to folder where output files from this and previous steps should be saved.
#' @param channel_file Path to a ".csv" file containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].
#' @param sample_file Path to a ".csv" file containing a Sample_ID column and additional grouping columns for the samples.
#' @param groups_column Column name in \code{sample_data} where group identifications are stored.
#' @param g1 Group label or vector with group labels for first group.
#' @param g2 Group label or vector with group labels for second group.
#' @param max_pval Maximum p-value. Used to filter phenotypes in final output.
#' @param parent_phen Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.
#' @param n_phenotypes maximum number of phenotypes to be considered from \code{phen_data} filtered by lowest p-values. Default: 1000.
#' @param min_confidence Minimal confidence threshold to filter output. Default: 0.5.
#' @param n_threads Number of threads to be used. Default: 1.
#' @param verbose If TRUE, print outputs from log to stdout.
#' 
#' @return Output is saved to file.
#' 
#' @export
get_independent_relevant_phenotypes_server <- function(output_folder,
                                                       channel_file,
                                                       sample_file,
                                                       groups_column,
                                                       g1,
                                                       g2,
                                                       max_pval = 0.05,
                                                       parent_phen = NULL,
                                                       n_phenotypes = 1000,
                                                       min_confidence = 0.5,
                                                       n_threads = 1,
                                                       verbose = FALSE
){
  
  
  
  
  filtered_phenotypes_file <- paste("significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,".csv",collapse = "",sep="")
  final_phenotypes_file <- paste("independent_significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,".csv",collapse = "",sep="")
  log_file <- paste("independent_significant_phenotypes_",ifelse(!is.null(parent_phen),paste("parent_",parent_phen,"_",sep=""),""),groups_column,"_",g1,"_vs_",g2,"_pval_",max_pval,".log",collapse = "",sep="")
  
  # Start logging
  
  start_log(log_path = file.path(output_folder,log_file), append = FALSE, verbose = verbose) 
  
  print_log("Reading data...")
  
  channel_data <- as.data.frame(data.table::fread(channel_file,check.names = FALSE))
  sample_data <- as.data.frame(data.table::fread(sample_file,check.names = FALSE))
  
  phen_data <- as.data.frame(data.table::fread(file.path(output_folder,filtered_phenotypes_file),nThread = n_threads))
  
  final_phenotypes <- get_independent_relevant_phenotypes(phen_data,
                                                          channel_data,
                                                          n_phenotypes = n_phenotypes,
                                                          min_confidence = min_confidence,
                                                          n_threads = n_threads
                                                          )
  
  print_log("Independent phenotypes found: ", nrow(final_phenotypes_file))
  
  print_log("Writing ",nrow(final_phenotypes), " phenotypes to file: ", file.path(output_folder,final_phenotypes_file))
  
  data.table::fwrite(final_phenotypes,file.path(output_folder,final_phenotypes_file))
  
  print_log("Done.")
  
  stop_log()
  
  
}