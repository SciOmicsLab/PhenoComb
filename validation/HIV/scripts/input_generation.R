
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(flowCore)
library(flowMeans)
library(data.table)
library(matrixStats)
library(Biobase)

input_folder <- "../data"
output_folder <- "../PhenoCombAnalysis/input"

# Get all files and file name as sample ID
fcs_files <- list.files(input_folder, pattern = ".fcs")
sample_ids <- unlist(lapply(fcs_files, function(i) as.integer(substr(i,1,nchar(i)-4))))




# Read samples and add sample ID
for(i in 1:length(sample_ids)){
  
  print(paste("Sample",i))
  
  # Read data
  fcs_sample_data <- flowCore::read.FCS(file.path(input_folder,fcs_files[i]),truncate_max_range = FALSE)
  
  # Get flow data and add Sample_ID col
  n_cells <- nrow(fcs_sample_data@exprs)
  sample_id_col <- matrix(rep(sample_ids[i],n_cells),nrow = n_cells, ncol = 1, dimnames = list(NULL,"Sample_ID"))
  fcs_sample_data <- fr_append_cols(fcs_sample_data,sample_id_col)
  
  # Write data back
  invisible(write.FCS(fcs_sample_data,file.path(input_folder,fcs_files[i])))
  #fwrite(as.data.table(fcs_flow_data),file.path(output_folder,"concat_data.csv"), showProgress = F, append = i>1)
  
  
  # Extract sample data
  sample_data <- data.frame(Sample_ID = sample_ids[i],
                            experiment_name = fcs_sample_data@description$`EXPERIMENT NAME`,
                            survival_time_from_seroconversion = fcs_sample_data@description$`CD Survival time from seroconversion`,
                            seroconversion_date = fcs_sample_data@description$`CD Seroconversion Datae`,
                            time_seroconv_to_sample = fcs_sample_data@description$`CD Time from seroc to sample`,
                            first_viral_load = fcs_sample_data@description$`CD First Viral Load`,
                            first_viral_load_date = fcs_sample_data@description$`CD First Viral Load Date`,
                            death = fcs_sample_data@description$`CD Event Censor`
  )
  
  # Write sample data
  fwrite(as.data.table(sample_data),file.path(output_folder,"sample_data.csv"), showProgress = F, append = i>1)
  
  rm(fcs_sample_data,fcs_flow_data,sample_data)
  gc(verbose = F, full = T)
  
}


fcs_sample_data <- flowCore::read.FCS(file.path(input_folder,fcs_files[1]),truncate_max_range = FALSE)
# Extract channel and marker names from description
marker_channels <- unlist(fcs_sample_data@description[paste("$P",4:16,"N",sep = "")])
marker_names <- unlist(fcs_sample_data@description[paste("$P",4:16,"S",sep = "")])

rm(fcs_sample_data)
gc(verbose = F, full = T)


# Generate thresholds

channel_data <- data.frame(
  Channel = marker_channels,
  Marker = marker_names
)

write.csv(channel_data,file.path(output_folder,"threshold_data.csv"))


# Read data back (I know, it's silly but in this way I can run the script in steps)
concatenated_data <- fread(file.path(output_folder,"concat_data.csv"))
concatenated_data <- as.matrix(concatenated_data)
gc(verbose = F, full = T)



# Generate FCS object
meta <- data.frame(name=colnames(concatenated_data),
                   desc=c(colnames(concatenated_data)[1:3],marker_names,"Sample_ID")
)


meta$minRange <- colMins(concatenated_data)
meta$maxRange <- colMaxs(concatenated_data)
meta$range <- meta$maxRange - meta$minRange


# a flowFrame is the internal representation of a FCS file
ff <- new("flowFrame",
          exprs=concatenated_data,
          parameters=AnnotatedDataFrame(meta)
)

# Write FCS file

invisible(write.FCS(ff,file.path(output_folder,"concat_data.fcs")))


