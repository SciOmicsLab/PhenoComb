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


# Test code

test_that("Processed cell data has correct number of rows", {
  expect_equal(nrow(processed_cell_data),nrow(cell_data))
})

test_that("Processed cell data all marker values are integer", {
  expect_equal(sum(processed_cell_data[1:n_markers] - floor(processed_cell_data[1:n_markers])),0)
})

test_that("Processed cell data all marker values match number of thresholds", {
  expect_true(all(unlist(lapply(processed_cell_data[,1:n_markers], function(i) all(i %in% c(0,1))))))
})





# Clean-up
options(PhenoComb.verbose = TRUE)

rm(cell_data,
   channel_data,
   sample_data,
   n_markers,
   n_samples,
   processed_cell_data)