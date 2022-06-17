
# Uncomment if running on RStudio
#library("rstudioapi")
#setwd(dirname(getActiveDocumentContext()$path))

library(rlang)
library(flowCore)
library(PhenoComb)

# Dataset info
input_folder <- "../PhenoCombAnalysis/input"
output_folder <- "../PhenoCombAnalysis/output"
cell_file <- "concat_1.fcs" # or "cell_data.csv"
channel_file <- "threshold_data.csv"
sample_file <- "ncell_filtered_sample_data.csv"

#Bootstraped dataset file name
resampled_cell_file <- "subsampled_concat_1.fcs"

# Output folder with final results
output_folder <- "../BootstrapAnalysis/outputs"

# Number of replicates
replicates <- 100

dir.create(output_folder, showWarnings = FALSE)
full_logfile <- file.path(output_folder, "bootstrap_full_record.log")
file.create(full_logfile)


# Read cell data
cell_data_ff <- flowCore::read.FCS(file.path(input_folder,cell_file),truncate_max_range = FALSE)

# Resampling function
resample_flow_data <- function(flow_frame, output_file){
  
  resampled_cell_data_ff <- duplicate(flow_frame, shallow = FALSE)
  
  flowCore::exprs(resampled_cell_data_ff) <- flowCore::exprs(resampled_cell_data_ff)[base::sample(1:nrow(flowCore::exprs(resampled_cell_data_ff)), 
                                                                                                    size = nrow(flowCore::exprs(resampled_cell_data_ff)), replace = TRUE), ]
  
  write.FCS(resampled_cell_data_ff,output_file)
}


# Start bootstrap loop
for(j in 1:replicates){

  cat(paste('Bootstrap iteration', j, '\n'), file = full_logfile, append = TRUE)

  # Resample data
  resample_flow_data(cell_data_ff,file.path(input_folder,resampled_cell_file))
  
  # Process raw data and generates all combinatorial phenotypes
  combinatorial_phenotype_counts_server(file.path(input_folder,resampled_cell_file),
                                        file.path(input_folder,channel_file),
                                        file.path(input_folder,sample_file),
                                        output_folder,
                                        parent_phen = NULL,
                                        min_count = 10,
                                        sample_fraction_min_counts = 0.5,
                                        max_phenotype_length = 0,
                                        sampleID_col = "Sample_ID",
                                        save_cell_data = TRUE,
                                        continue = FALSE,
                                        verbose = TRUE,
                                        max_ram = 0,
                                        efficient = TRUE,
                                        n_threads = 50
  )
  # Collect logs
  file.append(full_logfile, file.path(output_folder,"combinatorial_phenotypes.log" ) )
  
  # Filter statistically relevant phenotypes
  
  statistically_relevant_phenotypes_server(output_folder,
                                           file.path(input_folder,channel_file),
                                           file.path(input_folder,sample_file),
                                           output_file = "significant_phenotypes_CD3+_parent.csv",
                                           log_file = "significant_phenotypes_CD3+_parent.log",
                                           test_type = "survival",
                                           survival_time_column = "survival_time_from_seroconversion",
                                           survival_status_column = "death",
                                           parent_phen = "CD3+",
                                           max_pval = 1.0,
                                           continue = FALSE,
                                           n_threads = 50,
                                           verbose = TRUE)
  # Collect logs
  file.append(full_logfile, file.path(output_folder,"significant_phenotypes_CD3+_parent.log" ) )
  
  # Compute independent statistically relevant phenotypes
  get_independent_relevant_phenotypes_server(output_folder,
                                             file.path(input_folder,channel_file),
                                             file.path(input_folder,sample_file),
                                             input_significant_phenotypes = "significant_phenotypes_CD3+_parent.csv",
                                             output_file = "independent_phenotypes.csv",
                                             log_file = "independent_phenotypes.log",
                                             n_phenotypes = 5000,
                                             min_confidence = 0.0,
                                             max_pval = 0.0000005,
                                             n_threads = 50,
                                             verbose = TRUE
  )
  # Collect logs
  file.append(full_logfile, file.path(output_folder,"independent_phenotypes.log" ) )
  
  # Copy output from temporary folder to output folder
  file.copy(file.path(output_folder,"independent_phenotypes.csv"),file.path(output_folder,paste("independent_phenotypes_",j,".csv",sep = "")))
  

  
  gc(full = TRUE,verbose = FALSE)
  
}



