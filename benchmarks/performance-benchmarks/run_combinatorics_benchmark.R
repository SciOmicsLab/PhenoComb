
# Uncomment if using RStudio
# library("rstudioapi")
# setwd(dirname(getActiveDocumentContext()$path))

library(PhenoComb)

source("synthetic_dataset_generation.R")

# Read experiment definitions
experiments <- read.csv("experiments.csv")

replicates <- 3

# Creates temporary folders for generated inputs and temporary outputs
input_folder <- tempdir()

tmp_output_folder <- tempdir()

output_folder <- "outputs"

dir.create(output_folder)

# Run experiments

for(i in 1:nrow(experiments)){
  
  print(paste("Running experiment",i))
  
  experiment_id <- experiments[i,"experiment_ID"]
  n_markers <- experiments[i,"markers"]
  max_phen_len <- experiments[i,"max_phen_len"] 
  n_cells <- experiments[i,"cells"]
  n_samples <- experiments[i,"samples"]/2
  n_threads <- experiments[i,"threads"]
  
  
  generate_data_set(input_folder,
                    n_markers = n_markers,
                    marker_sd = 0.5,
                    marker_states = rep(2,n_markers),
                    cells_per_sample = n_cells,
                    group1_n_samples = n_samples,
                    group2_n_samples = n_samples,
                    n_threads = 50
  )
  
  for(j in 1:replicates){
    
    combinatorial_phenotype_counts_server(file.path(input_folder,"synthetic_data.fcs"),
                                          file.path(input_folder,"channel_data.csv"),
                                          file.path(input_folder,"sample_data.csv"),
                                          tmp_output_folder,
                                          parent_phen = NULL,
                                          min_count = 1,
                                          max_phenotype_length = max_phen_len,
                                          sampleID_col = "Sample_ID",
                                          save_cell_data = FALSE,
                                          continue = FALSE,
                                          verbose = FALSE,
                                          max_ram = 0,
                                          efficient = TRUE,
                                          n_threads = n_threads
    )
    
    
    file.copy(file.path(tmp_output_folder,"combinatorial_phenotypes.log"),file.path(output_folder,paste("experiment_",i,"_replicate_",j,".log",sep = "")))
    
    gc(full = TRUE,verbose = FALSE)
    
  }
}






