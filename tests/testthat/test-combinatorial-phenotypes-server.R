

# Path to test data
# This test dataset comprises one cell for each sample for each combination
# of full-length phenotype.
# n_markers = 4
# n_samples = 4

cell_file <- "../testdata/combinatorial_phenotypes/combinatorics_test_data.csv"
channel_file <- "../testdata/combinatorial_phenotypes/channel_data.csv"
sample_file <- "../testdata/combinatorial_phenotypes/sample_data.csv"

n_markers <- 4
n_samples <- 4


tmp_folder <- tempdir(check = TRUE)


test_that("Generates correct number combinatorial phenotypes", {
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 1,
                                        continue = FALSE,
                                        verbose = FALSE)
  
  combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  
  expect_equal(nrow(combinatorial_phenotypes), 81)
})

test_that("Generates correct counts for combinatorial phenotypes", {
  # Counts must be equal to 2^(total_markers-n_markers)
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 1,
                                        continue = FALSE,
                                        verbose = FALSE)
  
  combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  
  row_tests <- unlist(lapply(1:nrow(combinatorial_phenotypes), function(i){
    
    considered_markers <- count_markers(as.numeric(combinatorial_phenotypes[i,1:n_markers]))
    
    return(all(combinatorial_phenotypes[i,(n_markers+1):(n_markers+n_samples)] == 2^(n_markers-considered_markers)))
    
  }))
  
  expect_true(all(row_tests))
})


test_that("Generates correct number combinatorial phenotypes using min counts", {
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 5,
                                        efficient = FALSE,
                                        continue = FALSE,
                                        verbose = FALSE)
  
  filtered_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  
  expect_equal(nrow(filtered_combinatorial_phenotypes), 9)
  
  col_tests <- unlist(lapply((n_markers+1):(n_markers+n_samples), function(i){
    
    return(all(filtered_combinatorial_phenotypes[,i] >= 5))
    
  }))
  
  expect_true(all(col_tests))
  
})

test_that("Filters correctly for parent phenotypes", {
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 1,
                                        parent_phen = "Marker1+",
                                        continue = FALSE,
                                        verbose = FALSE)
  
  parent_filtered_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  expect_true(all(parent_filtered_combinatorial_phenotypes[,"Marker1"]==1))
  expect_equal(nrow(parent_filtered_combinatorial_phenotypes),27)
  
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 1,
                                        parent_phen = "Marker1+Marker2-",
                                        continue = FALSE,
                                        verbose = FALSE)
  
  parent_filtered_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  expect_true(all(parent_filtered_combinatorial_phenotypes[,"Marker1"] == 1) & all(parent_filtered_combinatorial_phenotypes[,"Marker2"] == 0))
  expect_equal(nrow(parent_filtered_combinatorial_phenotypes),9)
})


test_that("Phenotype length limit is correct", {
  for(n in 1:4){
    combinatorial_phenotype_counts_server(cell_file,
                                          channel_file,
                                          sample_file,
                                          tmp_folder,
                                          min_count = 0,
                                          max_phenotype_length = n,
                                          continue = FALSE,
                                          verbose = FALSE)
    
    length_limited_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
    
    row_tests <- unlist(lapply(1:nrow(length_limited_combinatorial_phenotypes), function(i){
      
      considered_markers <- count_markers(as.numeric(length_limited_combinatorial_phenotypes[i,1:n_markers]))
      
      return(considered_markers <= n)
      
    }))
    
    expect_true(all(row_tests))
  }
  
})


test_that("Non-Efficient mode produces same result", {
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 1,
                                        efficient = FALSE,
                                        continue = FALSE,
                                        verbose = FALSE)
  
  non_efficient_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  expect_equal(nrow(non_efficient_combinatorial_phenotypes), 81)
})


test_that("Multi-thread works", {
  combinatorial_phenotype_counts_server(cell_file,
                                        channel_file,
                                        sample_file,
                                        tmp_folder,
                                        min_count = 0,
                                        n_threads = 2,
                                        continue = FALSE,
                                        verbose = FALSE)
  
  parallel_combinatorial_phenotypes <- read.csv(file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))
  
  expect_equal(nrow(parallel_combinatorial_phenotypes), 81)
  
  row_tests <- unlist(lapply(1:nrow(parallel_combinatorial_phenotypes), function(i){
    
    considered_markers <- count_markers(as.numeric(parallel_combinatorial_phenotypes[i,1:n_markers]))
    
    return(all(parallel_combinatorial_phenotypes[i,(n_markers+1):(n_markers+n_samples)] == 2^(n_markers-considered_markers)))
    
  }))
  
  expect_true(all(row_tests))
})





# Clean-up

unlink(tmp_folder, recursive = TRUE)

rm(cell_file,
   channel_file,
   sample_file,
   n_markers,
   n_samples)