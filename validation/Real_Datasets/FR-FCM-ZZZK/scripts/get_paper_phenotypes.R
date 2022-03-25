
library(PhenoComb)

marker_data <- read.csv("../PhenoCombAnalysis/input/threshold_data.csv")

markers <- marker_data$Marker

significant_phenotypes_file <- "../PhenoCombAnalysis/output/significant_phenotypes_CD3+_parent.csv"

paper_final_phenotypes <- c("KI67+CD127-","CD45RO-CD8+CD57+CCR5-CD27+CCR7-CD127-","CD28-CD45RO+CD57-")

paper_original_phenotypes <- c("KI67+CD4-CCR5+CD127-","CD45RO-CD8+CD4-CD57+CCR5-CD27+CCR7-CD127-","CD28-CD45RO+CD4-CD57-CD27-CD127-")


#Add parent population
parent_population <- "CD3+"
paper_final_phenotypes <- paste(parent_population, paper_final_phenotypes, sep = "")
paper_original_phenotypes <- paste(parent_population, paper_original_phenotypes, sep = "")

paper_final_phenotypes_data <- do.call(rbind,lapply(paper_final_phenotypes, function(phen) find_phenotype_in_file(significant_phenotypes_file,phen,markers,50)))

if(!is.null(paper_final_phenotypes_data)){
  paper_final_phenotypes_data$n_markers <-  unlist(lapply(1:nrow(paper_final_phenotypes_data), function(i) count_markers(as.numeric(paper_final_phenotypes_data[i,markers]))))
  paper_final_phenotypes_data$Phenotype <-  unlist(lapply(1:nrow(paper_final_phenotypes_data), function(i) make_phenotype_name(as.numeric(paper_final_phenotypes_data[i,markers]),markers)))
  
  total_cols <- ncol(paper_final_phenotypes_data)
  
  paper_final_phenotypes_data <- paper_final_phenotypes_data[,c(total_cols:(total_cols-3),(length(markers)+1):(total_cols-4))]
  
  write.csv(paper_final_phenotypes_data,"../PhenoCombAnalysis/output/paper_final_phenotypes_data_CD3+_parent.csv",row.names = F)  
}


paper_original_phenotypes_data <- do.call(rbind,lapply(paper_original_phenotypes, function(phen) find_phenotype_in_file(significant_phenotypes_file,phen,markers,50)))

if(!is.null(paper_original_phenotypes_data)){
  paper_original_phenotypes_data$n_markers <-  unlist(lapply(1:nrow(paper_original_phenotypes_data), function(i) count_markers(as.numeric(paper_original_phenotypes_data[i,markers]))))
  paper_original_phenotypes_data$Phenotype <-  unlist(lapply(1:nrow(paper_original_phenotypes_data), function(i) make_phenotype_name(as.numeric(paper_original_phenotypes_data[i,markers]),markers)))
  paper_original_phenotypes_data <- paper_original_phenotypes_data[,c(total_cols:(total_cols-3),(length(markers)+1):(total_cols-4))]
  
  write.csv(paper_original_phenotypes_data,"../PhenoCombAnalysis/output/paper_original_phenotypes_dataCD3+_parent.csv",row.names = F)  
}


