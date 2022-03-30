
library(PhenoComb)
library(flowCore)

options(PhenoComb.verbose = TRUE)

# Paths to input files
cell_file <- "path/to/cell/file" # .csv or .fcs
channel_file <- "path/to/channel_file" # .csv
sample_file <- "path/to/sample_file" # .csv


# Reading data
# Uncoment if it's a .fcs file
cell_data <- as.data.frame(flowCore::read.FCS(cell_file,truncate_max_range = FALSE)@exprs)
#cell_data <- read.csv(cell_file) # comment if it's a .fcs file

channel_data <- read.csv(channel_file)
sample_data <- read.csv(sample_file)



# Process raw dataset
cell_data_processed <- process_cell_data(cell_data,
                                         channel_data,
                                         sample_data,
                                         sampleID_col = "Sample_ID",
                                         n_threads = 10
                                         )

# Generate combinatorial phenotypes
comb_phenotypes <- combinatorial_phenotype_counts(cell_data_processed,
                                                  parent_phen = NULL,
                                                  max_phenotype_length = 4,
                                                  min_count = 10,
                                                  efficient = TRUE,
                                                  n_threads = 10
                                                  )



# Filter statistically relevant phenotypes (Change groups_colum, g1, g2, and parent_phen accordingly to your data)
relevant_phenotypes <- compute_statistically_relevant_phenotypes(comb_phenotypes,
                                                                   channel_data,
                                                                   sample_data,
                                                                   test_type = "group",
                                                                   groups_column = "Group",
                                                                   g1 = "g1",
                                                                   g2 = "g2",
                                                                   max_pval = 0.05,
                                                                   parent_phen = NULL,
                                                                   n_threads = 10
)

# Example for correlation
# relevant_phenotypes <- compute_statistically_relevant_phenotypes(comb_phenotypes,
#                                                                  channel_data,
#                                                                  sample_data,
#                                                                  test_type = "correlation",
#                                                                  correlation_column = "Corr_data",
#                                                                  max_pval = 0.05,
#                                                                  parent_phen = NULL,
#                                                                  n_threads = 10
# )
# 
# Example for survival analysis
# relevant_phenotypes <- compute_statistically_relevant_phenotypes(comb_phenotypes,
#                                                                    channel_data,
#                                                                    sample_data,
#                                                                    test_type = "survival",
#                                                                    survival_time_column = "Survival_time",
#                                                                    survival_status_column = "Survival_status",
#                                                                    max_pval = 0.05,
#                                                                    parent_phen = NULL,
#                                                                    n_threads = 10
# )

# Get independent statistically relevant phenotypes
final_phenotypes <- get_independent_relevant_phenotypes(relevant_phenotypes,
                                                        channel_data,
                                                        n_phenotypes = 1000,
                                                        min_confidence = 0.0,
                                                        n_threads = 10
)

