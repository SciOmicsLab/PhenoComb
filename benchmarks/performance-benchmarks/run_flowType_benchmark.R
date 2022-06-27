#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Uncomment if using RStudio
#library("rstudioapi")
#setwd(dirname(getActiveDocumentContext()$path)

require(flowType, quietly=TRUE) # warn.conflicts = F,

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied for n_markers", call.=FALSE)
  # how to add optional parameters  
  #} else if (length(args)==1) {
  # default output file
  #  args[2] = "out.txt"
} else {
  args
}

N <- as.numeric(args[1])
print ("Number of markers", N)

LOGFILE=paste0("./outputs/flowType_benchmark_", N, ".log")
print( paste("Processing flowType benchmark for", N, "markers"))


print_log <- function(...){
  # Modified to only echo to stdout, not a logfile, but with timestamp/log format  
    time <- paste(c("[", format(Sys.time(), "%y/%m/%d %X"),"] "), collapse = "")
    txt <- paste(c(time, "[INFO] ", ...), collapse = "")
    
    message(txt)
}


time <- paste(c("[", format(Sys.time(), "%y/%m/%d %X"),"] "), collapse = "")
print_log("Starting flowType benchmark for", N, "markers")



source("synthetic_dataset_generation.R")
plot_res_bar <- function (res) {
  ###################################################
  ### code chunk number 5: BarPlot
  ###################################################
  MFIs=res@MFIs;
  Proportions=res@CellFreqs;
  Proportions <- Proportions / max(Proportions)
  names(Proportions) <- unlist(lapply(res@PhenoCodes, 
                                      function(x){return(decodePhenotype(
                                        x,res@MarkerNames[PropMarkers],
                                        res@PartitionsPerMarker))}))
  rownames(MFIs)=names(Proportions)
  index=order(Proportions,decreasing=TRUE)[1:20]
  bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
  text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0,
       labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
  axis(2);
  axis(1, at=bp, labels=FALSE);
  title(xlab='Phenotype Names', ylab='Cell Proportion')
} 


# how many repeats to do, to average the times
replicates <- 3

# Creates temporary folders for generated inputs and temporary outputs
#input_folder <- tempdir()
input_folder <- "inputs"
output_folder <- "outputs"
tmp_output_folder <- tempdir()

dir.create(input_folder, showWarnings = FALSE)
dir.create(output_folder, showWarnings = FALSE)

# Run experiments

#for(i in 1:nrow(experiments)){
#for(i in 1:2){
  i=1
  print(paste("Running experiment",i))
  
  # experiment_id <- experiments[i,"experiment_ID"]
  # n_markers <- experiments[i,"markers"]
  # max_phen_len <- experiments[i,"max_phen_len"] 
  # n_cells <- experiments[i,"cells"]
  # n_samples <- experiments[i,"samples"]/2
  # n_threads <- experiments[i,"threads"]
  
  experiment_id <- N-4
  n_markers <- N
  max_phen_len <- N
  n_cells <- 40000
  n_samples <- 10/2
  n_threads <- 50
  
  
  generate_data_set(input_folder,
                    n_markers = n_markers,
                    marker_sd = 0.5,
                    marker_states = rep(2,n_markers),
                    cells_per_sample = n_cells,
                    group1_n_samples = n_samples,
                    group2_n_samples = n_samples,
                    n_threads = 50
  )
  gc(full = TRUE,verbose = FALSE)

  
  experiment_id <- N-4
  n_markers <- N
  max_phen_len <- N
  n_cells <- 40000
  n_samples <- 10/2
  n_threads <- 50
  
  PropMarkers <- 1:n_markers
  MFIMarkers <- PropMarkers
  Thresholds <- as.list(rep(1.5, n_markers))

  MarkerNames <- sapply(1:n_markers, FUN=function(x) paste0("Marker",x))
  # e.g. MarkerNames= c('Marker1','Marker2','Marker3','Marker4')
  print_log("Marker Names ", MarkerNames)
  
  fcsData <- read.FCS("inputs/synthetic_data.fcs")

  for(j in 1:replicates){

    print_log("Starting flowType, replicate ", j)
    time1 <- format(Sys.time(), "%s")
    print(as.integer(time1)) 
     # This needs a wrapper to handle logging, times, and files
     res <- flowType(fcsData, PropMarkers, MFIMarkers, Methods='thresholds',
              MarkerNames=MarkerNames, Thresholds=Thresholds, MarkerNames,
              MaxMarkersPerPop=max_phen_len, 
              PartitionsPerMarker=2, verbose=TRUE)  
   
     head(res@Partitions)
     ##########################

     print_log("Starting flowType, replicate ", j)
     time2 <- format(Sys.time(), "%s")

     print_log("Elapsed Time ", N, " markers:", as.integer(time2) - as.integer(time1))
    
    gc(full = TRUE,verbose = FALSE)
    
  }
#}







