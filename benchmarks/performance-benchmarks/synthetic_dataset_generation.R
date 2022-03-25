

library(PhenoComb)
library(Biobase)
library(flowCore)
library(parallel)


generate_proportions <- function(marker_combinations,
                                 group1_perturbed_phenotypes = NULL,
                                 group2_perturbed_phenotypes = NULL,
                                 group1_phenotype_perturbations = NULL,
                                 group2_phenotype_perturbations = NULL,
                                 n_threads = 1){
  
  proportions_g1 <- rep(1,nrow(marker_combinations))
  proportions_g2 <- rep(1,nrow(marker_combinations))
  
  if(!is.null(group1_perturbed_phenotypes)){
    for(i in 1:length(group1_perturbed_phenotypes)){
      has_phen <- unlist(mclapply(1:nrow(marker_combinations), function(j) has_phenotype(marker_combinations[j,],group1_perturbed_phenotypes[i]),mc.cores = n_threads))
      proportions_g1[has_phen] <- proportions_g1[has_phen]+group1_phenotype_perturbations[i]
    }
    
  }
  
  if(!is.null(group2_perturbed_phenotypes)){
    for(i in 1:length(group2_perturbed_phenotypes)){
      has_phen <- unlist(mclapply(1:nrow(marker_combinations), function(j) has_phenotype(marker_combinations[j,],group2_perturbed_phenotypes[i]),mc.cores = n_threads))
      proportions_g2[has_phen] <- proportions_g2[has_phen]+group2_phenotype_perturbations[i]
    }
    
    
  }
  
  #normalization_factor <- sum(c(proportions_g1,proportions_g2))
  
  proportions_g1[proportions_g1 < 0] <- 0
  proportions_g2[proportions_g2 < 0] <- 0
  
  proportions_g1 <- proportions_g1 / sum(proportions_g1)
  proportions_g2 <- proportions_g2 / sum(proportions_g2)
  
  return(list(g1 = proportions_g1, g2 = proportions_g2))
  
}

generate_cells <- function(marker_combinations,
                           proportions,
                           marker_sd,
                           total_cells,
                           sample_name
){
  
  cells <- marker_combinations[sample(1:nrow(marker_combinations), size = n_cells, prob = proportions, replace = TRUE),]
  
  cells <- as.data.frame(apply(cells,c(1,2), function(i) rnorm(1,i,marker_sd)))
  
  cells$Sample_ID <- sample_name
  
  return(cells)
  
}




write_channel_data <- function(markers,marker_states,output_folder){
  
  channel_data <- data.frame(
    Channel = paste("Ch",1:length(markers),sep = ""),
    Marker = markers
  )
  
  n_thresholds <- max(marker_states) - 1
  
  for(i in 1:n_thresholds){
    col_name <- paste("T",i,sep = "")
    channel_data[,col_name] <- NA
    channel_data[marker_states-1 >= i,col_name] <- i+0.5
  }
  
  write.csv(channel_data,file.path(output_folder,"channel_data.csv"), row.names = F, na = "")
}


write_sample_data <- function(group1_n_samples,group2_n_samples,output_folder){
  total_samples <- group1_n_samples+group2_n_samples
  sample_data <- data.frame(
    Sample_ID = 1:total_samples,
    Group = c(rep("g1",group1_n_samples),rep("g2",group2_n_samples))
  )
  write.csv(sample_data,file.path(output_folder,"sample_data.csv"), row.names = F, na = "")
}


generate_data_set <- function(
  output_folder,
  n_markers = 4,
  markers = NULL,
  marker_states = c(2,2,2,2),
  marker_sd = 0.3,
  cells_per_sample = 1000,
  group1_n_samples = 3,
  group2_n_samples = 3,
  group1_perturbed_phenotypes = NULL,
  group2_perturbed_phenotypes = NULL,
  group1_phenotype_perturbations = NULL,
  group2_phenotype_perturbations = NULL,
  n_threads = 1
){
  
  # Fix marker lengths or names
  
  if(length(markers)<1){
    markers <- paste("Marker",1:n_markers,sep = "")
  }else{
    n_markers <- length(markers)
  }
  
  # Input consistency checks
  
  if(n_markers<1) stop("The number of markers must be at least 1.")
  
  if(length(marker_states) != n_markers) stop("The length of marker states must match the number of markers.")
  
  if(any(marker_states<2)) stop("Each marker must have a minimum of 2 states.")
  
  if(group1_n_samples < 1 | group2_n_samples < 1) stop("Number of samples per group must be at least 1.")
  
  
  # Generate channel file
  
  write_channel_data(markers,marker_states,output_folder)
  
  # Generate sample file
  
  write_sample_data(group1_n_samples,group2_n_samples,output_folder)
  
  
  # Generate synthetic data
  
  
  marker_combinations <- expand.grid(lapply(marker_states, function(i) 1:i))
  colnames(marker_combinations) <- markers
  
  
  proportions <- generate_proportions(marker_combinations,
                                      group1_perturbed_phenotypes,
                                      group2_perturbed_phenotypes,
                                      group1_phenotype_perturbations,
                                      group2_phenotype_perturbations,
                                      n_threads)
  
  
  group_1_cells <- do.call(rbind,mclapply(1:group1_n_samples, function(i) cells <- generate_cells(marker_combinations,
                                                                                                proportions$g1,
                                                                                                marker_sd,
                                                                                                cells_per_sample,
                                                                                                i
  ),mc.cores = n_threads))
  
  group_2_cells <- do.call(rbind,mclapply(1:group2_n_samples, function(i) cells <- generate_cells(marker_combinations,
                                                                                                proportions$g2,
                                                                                                marker_sd,
                                                                                                cells_per_sample,
                                                                                                (i+group1_n_samples)
  ),mc.cores = n_threads))
  
  synthetic_data <- rbind(group_1_cells,group_2_cells)
  colnames(synthetic_data) <- c(paste("Ch",1:length(markers),sep = ""),"Sample_ID")
  
  
  # Generate FCS object
  meta <- data.frame(name=colnames(synthetic_data),
                     desc=paste('Synthetic',colnames(synthetic_data))
  )
  meta$range <- apply(apply(synthetic_data,2,range),2,diff)
  meta$minRange <- apply(synthetic_data,2,min)
  meta$maxRange <- apply(synthetic_data,2,max)
  
  # a flowFrame is the internal representation of a FCS file
  ff <- new("flowFrame",
            exprs=as.matrix(synthetic_data),
            parameters=AnnotatedDataFrame(meta)
  )
  
  # Write FCS file
  
  invisible(write.FCS(ff,file.path(output_folder,"synthetic_data.fcs")))
  
}





# folder <- "../.."
# 
# 
# generate_data_set(folder,
#                   n_markers = 5,
#                   marker_sd = 1.5,
#                   marker_states = rep(2,5),
#                   cells_per_sample = 2000,
#                   group1_n_samples = 4,
#                   group2_n_samples = 4,
#                   group1_perturbed_phenotypes = c("Marker1+Marker2+Marker3-"),
#                   group1_phenotype_perturbations = c(10)
#                   #group2_perturbed_phenotypes = c("Marker1-Marker2-"),
#                   #group2_phenotype_perturbations = c(10)
# )
