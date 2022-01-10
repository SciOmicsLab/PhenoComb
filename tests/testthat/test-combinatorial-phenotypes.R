

# Temporary verbose false
options(PhenoComb.verbose = FALSE)


# Read test data
# This test dataset comprises one cell for each sample for each combination
# of full-length phenotype.
# n_markers = 4
# n_samples = 4

cell_data <- read.csv("../testdata/combinatorial_phenotypes/combinatorics_test_data.csv")
channel_data <- read.csv("../testdata/combinatorial_phenotypes/channel_data.csv")
sample_data <- read.csv("../testdata/combinatorial_phenotypes/sample_data.csv")

n_markers <- nrow(channel_data)
n_samples <- nrow(sample_data)


# Generate outputs
processed_cell_data <- process_cell_data(cell_data,channel_data,sample_data)

combinatorial_phenotypes <- combinatorial_phenotype_counts(processed_cell_data, min_count = 0)

filtered_combinatorial_phenotypes <- combinatorial_phenotype_counts(processed_cell_data, min_count = 5)



# Test code

test_that("Generates correct number combinatorial phenotypes", {
  expect_equal(nrow(combinatorial_phenotypes), 81)
})

test_that("Generates correct number combinatorial phenotypes using min counts", {
  expect_equal(nrow(filtered_combinatorial_phenotypes), 9)
})

test_that("Generates correct counts for combinatorial phenotypes", {
  # Counts must be equal to 2^(total_markers-n_markers)
  row_tests <- unlist(lapply(1:nrow(combinatorial_phenotypes), function(i){
    
    considered_markers <- count_markers(as.numeric(combinatorial_phenotypes[i,1:n_markers]))
    
    return(all(combinatorial_phenotypes[i,(n_markers+1):(n_markers+n_samples)] == 2^(n_markers-considered_markers)))
    
  }))
  
  expect_true(all(row_tests))
})






# Clean-up
options(PhenoComb.verbose = TRUE)

rm(cell_data,
   channel_data,
   sample_data,
   n_markers,
   n_samples,
   processed_cell_data,
   combinatorial_phenotypes)