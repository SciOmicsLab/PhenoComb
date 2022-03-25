
library(PhenoComb)


# input_folder <- "/Users/burke/Projects/PhenoComb/validation/Real_Datasets/COVIDome/PhenoCombAnalysis/input"
# output_folder <- "/Users/burke/Projects/PhenoComb/validation/Real_Datasets/COVIDome/PhenoCombAnalysis/output"
input_folder <- "../PhenoCombAnalysis/input"
output_folder <- "../PhenoCombAnalysis/output"
cell_file <- "concat_1.fcs" # or "cell_data.csv"
channel_file <- "Thresholds.csv"
sample_file <- "Sample_Info.csv"



# Process raw data and generates all combinatorial phenotypes
combinatorial_phenotype_counts_server(file.path(input_folder,cell_file),
                                      file.path(input_folder,channel_file),
                                      file.path(input_folder,sample_file),
                                      output_folder,
                                      parent_phen = NULL,
                                      min_count = 10,
                                      sample_fraction_min_counts = 0.5,
                                      max_phenotype_length = 6,
                                      sampleID_col = "Sample_ID",
                                      save_cell_data = TRUE,
                                      continue = TRUE,
                                      verbose = TRUE,
                                      max_ram = 0,
                                      efficient = TRUE,
                                      n_threads = 54
)


# Filter statistically relevant phenotypes
#(Change groups_colum, g1, g2, and parent_phen accordingly to your data)
# statistically_relevant_phenotypes_server(output_folder,
#                                          file.path(input_folder,channel_file),
#                                          file.path(input_folder,sample_file),
#                                          test_type = "group",
#                                          groups_column = "COVID_status",
#                                          g1 = "Positive",
#                                          g2 = "Negative",
#                                          max_pval = 1.0,
#                                          continue = TRUE,
#                                          n_threads = 54,
#                                          verbose = TRUE)
# 
# 
# statistically_relevant_phenotypes_server(output_folder,
#                                          file.path(input_folder,channel_file),
#                                          file.path(input_folder,sample_file),
#                                          output_file = "significant_phenotypes_CD3+_parent.csv",
#                                          log_file = "significant_phenotypes_CD3+_parent.log",
#                                          test_type = "survival",
#                                          survival_time_column = "survival_time_from_seroconversion",
#                                          survival_status_column = "death",
#                                          parent_phen = "CD3+",
#                                          max_pval = 1.0,
#                                          continue = TRUE,
#                                          n_threads = 50,
#                                          verbose = TRUE)

# # Compute independent statistically relevant phenotypes
# # (Change groups_colum, g1, g2, and parent_phen accordingly to your data)
get_independent_relevant_phenotypes_server(output_folder,
                                           file.path(input_folder,channel_file),
                                           file.path(input_folder,sample_file),
                                           n_phenotypes = 5000,
                                           min_confidence = 0.0,
                                           max_pval = 0.0000005,
                                           n_threads = 50,
                                           verbose = TRUE
)
