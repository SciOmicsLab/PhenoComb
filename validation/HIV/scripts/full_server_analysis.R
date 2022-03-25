
library(PhenoComb)


#input_folder <- "/Users/burke/Projects/PhenoComb/validation/Real_Datasets/FR-FCM-ZZZK/PhenoCombAnalysis/input"
#output_folder <- "/Users/burke/Projects/PhenoComb/validation/Real_Datasets/FR-FCM-ZZZK/PhenoCombAnalysis/output"
input_folder <- "../PhenoCombAnalysis/input"
output_folder <- "../PhenoCombAnalysis/output"
cell_file <- "concat_1.fcs" # or "cell_data.csv"
channel_file <- "threshold_data.csv"
sample_file <- "ncell_filtered_sample_data.csv"



# Process raw data and generates all combinatorial phenotypes
combinatorial_phenotype_counts_server(file.path(input_folder,cell_file),
                                      file.path(input_folder,channel_file),
                                      file.path(input_folder,sample_file),
                                      output_folder,
                                      parent_phen = NULL,
                                      min_count = 10,
                                      sample_fraction_min_counts = 0.5,
                                      max_phenotype_length = 0,
                                      sampleID_col = "Sample_ID",
                                      save_cell_data = TRUE,
                                      continue = TRUE,
                                      verbose = TRUE,
                                      max_ram = 0,
                                      efficient = TRUE,
                                      n_threads = 50
)

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
                                         continue = TRUE,
                                         n_threads = 50,
                                         verbose = TRUE)

# Compute independent statistically relevant phenotypes
get_independent_relevant_phenotypes_server(output_folder,
                                           file.path(input_folder,channel_file),
                                           file.path(input_folder,sample_file),
                                           #input_significant_phenotypes = "significant_phenotypes_CD3+_parent.csv",
                                           output_file = "independent_phenotypes_pval_0.0000005_5000_phenotypes.csv",
                                           log_file = "independent_phenotypes_pval_0.0000005_5000_phenotypes.log",
                                           n_phenotypes = 5000,
                                           min_confidence = 0.0,
                                           max_pval = 0.0000005,
                                           n_threads = 50,
                                           verbose = TRUE
)
1