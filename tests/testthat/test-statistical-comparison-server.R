

# Read test data
# This test dataset has counts for 81 phenotypes
# Phenotype 1 is total cells
# Phenotypes 2:21 have equal counts for all groups
# Phenotypes 22:51 have g1 doubled
# Phenotypes 52:81 have g2 doubled
# n_markers = 4
# n_samples = 4

combinatorial_phenotypes <- "../testdata/statistical_comparison/fake_phenotype_counts.csv"
channel_file <- "../testdata/combinatorial_phenotypes/channel_data.csv"
sample_file <- "../testdata/combinatorial_phenotypes/sample_data.csv"

tmp_folder <- tempdir(check = TRUE)

file.copy(combinatorial_phenotypes,tmp_folder)

file.rename(file.path(tmp_folder,"fake_phenotype_counts.csv"),file.path(tmp_folder,"combinatorial_phenotype_counts.csv"))

n_markers <-   nrow(read.csv(channel_file))
n_samples <- nrow(read.csv(sample_file))




test_that("Perform statistical comparison correctly", {
  
  suppressWarnings(statistically_relevant_phenotypes_server(
    tmp_folder,
    channel_file,
    sample_file,
    test_type = "group",
    groups_column = "Group",
    g1 = "g1",
    g2 = "g2",
    max_pval = 1.
  ))
  
  relevant_phenotypes <- read.csv(file.path(tmp_folder,"significant_phenotypes.csv"))
  
  
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
  
  
  suppressWarnings(statistically_relevant_phenotypes_server(
    tmp_folder,
    channel_file,
    sample_file,
    test_type = "group",
    groups_column = "Group",
    g1 = "g1",
    g2 = "g2",
    max_pval = 0.5
  ))
  
  relevant_phenotypes <- read.csv(file.path(tmp_folder,"significant_phenotypes.csv"))
  
  expect_equal(nrow(relevant_phenotypes), 60)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  expect_true(all(relevant_phenotypes[1:30,"effect_size"] == 1))
  expect_true(all(relevant_phenotypes[1:30,"log2FoldChange"] == 1))
  
  expect_true(all(relevant_phenotypes[31:60,"effect_size"] == -1))
  expect_true(all(relevant_phenotypes[31:60,"log2FoldChange"] == -1))
  
  expect_true(all(relevant_phenotypes[,"p_value"] <= 0.5))
  
})


test_that("Filters for parent phenotype", {
  
  suppressWarnings(statistically_relevant_phenotypes_server(
    tmp_folder,
    channel_file,
    sample_file,
    test_type = "group",
    groups_column = "Group",
    g1 = "g1",
    g2 = "g2",
    max_pval = 1.0,
    parent_phen = "Marker1+"
  ))
  
  relevant_phenotypes <- read.csv(file.path(tmp_folder,"significant_phenotypes.csv"))
  
  
  expect_equal(nrow(relevant_phenotypes), 27)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  suppressWarnings(statistically_relevant_phenotypes_server(
    tmp_folder,
    channel_file,
    sample_file,
    test_type = "group",
    groups_column = "Group",
    g1 = "g1",
    g2 = "g2",
    max_pval = 1.0,
    parent_phen = "Marker1+Marker2-"
  ))
  
  relevant_phenotypes <- read.csv(file.path(tmp_folder,"significant_phenotypes.csv"))
  
  expect_equal(nrow(relevant_phenotypes), 9)
  expect_equal(ncol(relevant_phenotypes),n_markers+n_samples+3)
  
  
  
})


# Clean-up

unlink(tmp_folder, recursive = TRUE)

rm(combinatorial_phenotypes,
   channel_file,
   sample_file,
   n_markers,
   n_samples)
