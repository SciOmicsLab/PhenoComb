
# Temporary verbose false
options(PhenoComb.verbose = FALSE)


# Read test data
# This test dataset has counts for 81 phenotypes
# Phenotype 1 is total cells
# Phenotypes 2:21 have equal counts for all groups
# Phenotypes 22:51 have g1 doubled
# Phenotypes 52:81 have g2 doubled
# n_markers = 4
# n_samples = 4

combinatorial_phenotypes <- read.csv("../testdata/statistical_comparison/fake_phenotype_counts.csv",check.names = FALSE)
channel_data <- read.csv("../testdata/combinatorial_phenotypes/channel_data.csv")
sample_data <- read.csv("../testdata/combinatorial_phenotypes/sample_data.csv")

n_markers <- nrow(channel_data)
n_samples <- nrow(sample_data)


test_that("Perform statistical comparison correctly", {
  
  relevant_phenotypes <- suppressWarnings(compute_statistically_relevant_phenotypes(combinatorial_phenotypes,
                                                                                    channel_data,
                                                                                    sample_data,
                                                                                    "Group",
                                                                                    "g1",
                                                                                    "g2",
                                                                                    max_pval = 1.
  ))
  
  expect_equal(nrow(relevant_phenotypes), 80)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  expect_true(all(relevant_phenotypes[1:20,"effect_size"] == 0))
  expect_true(all(relevant_phenotypes[1:20,"log2FoldChange"] == 0))
  
  expect_true(all(relevant_phenotypes[21:50,"effect_size"] == 1))
  expect_true(all(relevant_phenotypes[21:50,"log2FoldChange"] == 1))
  
  expect_true(all(relevant_phenotypes[51:80,"effect_size"] == -1))
  expect_true(all(relevant_phenotypes[51:80,"log2FoldChange"] == -1))
  
})

test_that("Filters by p-value correctly", {
  
  relevant_phenotypes <- suppressWarnings(compute_statistically_relevant_phenotypes(combinatorial_phenotypes,
                                                                                    channel_data,
                                                                                    sample_data,
                                                                                    "Group",
                                                                                    "g1",
                                                                                    "g2",
                                                                                    max_pval = 0.5
  ))
  
  expect_equal(nrow(relevant_phenotypes), 60)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  expect_true(all(relevant_phenotypes[1:30,"effect_size"] == 1))
  expect_true(all(relevant_phenotypes[1:30,"log2FoldChange"] == 1))
  
  expect_true(all(relevant_phenotypes[31:60,"effect_size"] == -1))
  expect_true(all(relevant_phenotypes[31:60,"log2FoldChange"] == -1))
  
  expect_true(all(relevant_phenotypes[,"p_value"] <= 0.5))
  
})


test_that("Filters for parent phenotype", {
  
  relevant_phenotypes <- suppressWarnings(compute_statistically_relevant_phenotypes(combinatorial_phenotypes,
                                                                                    channel_data,
                                                                                    sample_data,
                                                                                    "Group",
                                                                                    "g1",
                                                                                    "g2",
                                                                                    max_pval = 1.0,
                                                                                    parent_phen = "Marker1+"
  ))
  
  expect_equal(nrow(relevant_phenotypes), 26)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  relevant_phenotypes <- suppressWarnings(compute_statistically_relevant_phenotypes(combinatorial_phenotypes,
                                                                                    channel_data,
                                                                                    sample_data,
                                                                                    "Group",
                                                                                    "g1",
                                                                                    "g2",
                                                                                    max_pval = 1.0,
                                                                                    parent_phen = "Marker1+Marker2-"
  ))
  
  expect_equal(nrow(relevant_phenotypes), 8)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  
  
})

# Clean-up
options(PhenoComb.verbose = TRUE)

rm(combinatorial_phenotypes,
   channel_data,
   sample_data,
   n_markers,
   n_samples)