
# To run using Rscript

library(PhenoComb)


input_folder <- "../PhenoCombAnalysis/input"
output_folder <- "../PhenoCombAnalysis/output_lineages"
cell_file <- "concat_1.fcs" # or "cell_data.csv"
sample_file <- "Sample_Info.csv"


parent_populations <- c("CD45+CD19+","CD45+CD56+","CD45+CD14+","CD45+CD3-CD19-CD56-CD14-","CD45+CD3+")


for(pp in parent_populations){

  lineage_channel_file <- paste("Thresholds_",pp,".csv",sep = "")

  # Process raw data and generates all combinatorial phenotypes
  combinatorial_phenotype_counts_server(file.path(input_folder,cell_file),
                                        file.path(input_folder,lineage_channel_file),
                                        file.path(input_folder,sample_file),
                                        output_folder,
                                        parent_phen = pp,
                                        min_count = 100,
                                        sample_fraction_min_counts = 0.5,
                                        max_phenotype_length = 22,
                                        sampleID_col = "Sample_ID",
                                        save_cell_data = TRUE,
                                        continue = TRUE,
                                        verbose = TRUE,
                                        max_ram = 0,
                                        efficient = TRUE,
                                        n_threads = 54
  )

  file.rename(file.path(output_folder,"combinatorial_phenotypes.log"),file.path(output_folder,paste("combinatorial_phenotypes_parent_",pp,".log",sep = "")))
  file.rename(file.path(output_folder,"combinatorial_phenotype_counts.csv"),file.path(output_folder,paste("combinatorial_phenotype_counts_parent_",pp,".csv",sep = "")))


}



# Filter statistically relevant phenotypes
#(Change groups_colum, g1, g2, and parent_phen accordingly to your data)
for(pp in parent_populations){
  
  lineage_channel_file <- paste("Thresholds_",pp,".csv",sep = "")
  
  statistically_relevant_phenotypes_server(output_folder,
                                           file.path(input_folder,lineage_channel_file),
                                           file.path(input_folder,sample_file),
                                           input_phenotype_counts = paste("combinatorial_phenotype_counts_parent_",pp,".csv",sep = ""),
                                           input_phenotype_counts_log = paste("combinatorial_phenotypes_parent_",pp,".log",sep = ""),
                                           output_file = paste("significant_phenotypes_parent_",pp,".csv",sep = ""),
                                           log_file = paste("significant_phenotypes_parent_",pp,".log",sep = ""),
                                           test_type = "group",
                                           groups_column = "COVID_status",
                                           g1 = "Positive",
                                           g2 = "Negative",
                                           max_pval = 0.0005,
                                           continue = TRUE,
                                           n_threads = 54,
                                           verbose = TRUE)
  
  
  get_independent_relevant_phenotypes_server(output_folder,
                                             file.path(input_folder,lineage_channel_file),
                                             file.path(input_folder,sample_file),
                                             input_significant_phenotypes = paste("significant_phenotypes_parent_",pp,".csv",sep = ""),
                                             output_file = paste("independent_phenotypes_parent_",pp,".csv",sep = ""),
                                             log_file = paste("independent_phenotypes_parent_",pp,".log",sep = ""),
                                             n_phenotypes = 10000,
                                             min_confidence = 0.0,
                                             max_pval = 0.0005,
                                             n_threads = 54,
                                             verbose = TRUE
  )
  

}


