
library(PhenoComb)


input_folder <- "/path/to/input/folder"
output_folder <- "path/to/output/folder"
cell_file <- "cell_data.fcs" # or "cell_data.csv"
channel_file <- "channel_data.csv"
sample_file <- "sample_data.csv"


# Process raw data and generates all combinatorial phenotypes
combinatorial_phenotype_counts_server(file.path(input_folder,cell_file),
                                      file.path(input_folder,channel_file),
                                      file.path(input_folder,sample_file),
                                      output_folder,
                                      parent_phen = NULL,
                                      min_count = 10,
                                      max_phenotype_length = 0,
                                      sampleID_col = "Sample_ID",
                                      save_cell_data = TRUE,
                                      continue = TRUE,
                                      verbose = TRUE,
                                      max_ram = 0,
                                      efficient = TRUE,
                                      n_threads = 10
)

# Filter statistically relevant phenotypes
#(Change groups_colum, g1, g2, and parent_phen accordingly to your data)
statistically_relevant_phenotypes_server(output_folder,
                                         file.path(input_folder,channel_file),
                                         file.path(input_folder,sample_file),
                                         test_type = "group",
                                         groups_column = "Group",
                                         g1 = "g1",
                                         g2 = "g2",
                                         max_pval = 0.05,
                                         parent_phen = NULL,
                                         continue = TRUE,
                                         n_threads = 10,
                                         verbose = TRUE
)


# Compute independent statistically relevant phenotypes
get_independent_relevant_phenotypes_server(output_folder,
                                           file.path(input_folder,channel_file),
                                           file.path(input_folder,sample_file),
                                           n_phenotypes = 1000,
                                           min_confidence = 0.0,
                                           n_threads = 10,
                                           verbose = TRUE
)
