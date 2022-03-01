# PhenoComb

An R package for combinatorial analysis of phenotypes

## Installing PhenoComb

To install PhenoComb you will need to have the following packages installed:

-   devtools
-   BiocManager
-   flowCore

Copy and paste the following code into R to install them.

``` r
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

Then, to install PhenoComb copy and paste the following code into R:

``` r
devtools::install_github("SciOmicsLab/PhenoComb")
```

## Tutorial

This tutorial is intended to show the basic workflow and features of PhenoComb.

It has two workflow modes:

-   Local: having all your data in-memory, useful for small datasets (Up to 14 markers recommended).
-   Server: storing all outputs into files, aimed to process large datasets and number of markers.

Regardless of the desired workflow, the inputs are the same. However, for the server version, the inputs must be accessed directly from files rather than R objects.

### Input Preparation

PhenoComb requires three inputs:

1.  Cell Data (data.frame, .fcs, or .csv)
2.  Channel Data (data.frame or .csv)
3.  Sample Data (data.frame or .csv)

#### Cell Data

Cell Data can be your concatenated FCS file with the channels as columns and fluorescence values as rows for each cell. It must also contain one column with the sample identification. As a default, we suggest a column called "Sample_ID" (not to be confounded with SampleID standard column from FlowJo, **never use that**), but can be called as whatever the user desire.

EXAMPLE

#### Channel Data

Channel Data is a table containing the columns described below. Required columns are marked in bold.

-   **Channel**: channel name respective to its column in Cell Data. Names must match.
-   **Marker**: name of the marker measured in the respective Channel. Please avoid spaces and special characters.
-   **T1**: column with a threshold value used to discretize measured values into states (-/+).
-   T2, T3, ..., Tn: additional thresholds to discretize values into more than two states (-/++/+++/...).
-   OOB: aka Out Of Bounds values, used to discard cell measurement if measured value is above this threshold.

Channels not specified in this file will be discarded.

EXAMPLE

#### Sample Data

Sample Data is a table containing the names of your samples and meta information used to perform further statistical tests.

The required column is marked in bold.

-   **Sample_ID**: column with sample IDs contained in the respective column in Cell Data. The following columns are only examples. Can be named or represent data as the user convenience.
-   Sample_Group: group label respective to the sample.
-   Correlated_Measurement: a value for the phenotypes be correlated to.
-   ...

### Local Workflow

PhenoComb provide some tools that can be conveniently chained to obtain final filtered statistically relevant phenotypes. We reccomend using the Local Workflow to datasets with up to 14 markers. For more markers, consider using the Server Workflow.

The main tools are:

1.  Data Pre-Processing.
2.  Combinatorial Phenotype Cell Counting
3.  Statistical Testing.
4.  Independent Phenotype Filtering.

The following sections will demonstrate the basic use of each tool.

#### Reading Inputs

For the "local workflow", all data must be read into memory as data.frame objects. The following code sniped shows an example of how to do that:

``` r
# Use this line if your input is a FCS file (Must have flowCore installed).
cell_data <- as.data.frame(flowCore::read.FCS("/path/to/cell_data.fcs",truncate_max_range = FALSE)@exprs)

# Use this line if your input is a CSV file.
cell_data <- read.csv("/path/to/cell_data.csv") # comment if it's a .fcs file

channel_data <- read.csv("/path/to/channel_data.csv")
sample_data <- read.csv("/path/to/sample_data.csv")
```

All input files must be in accordance with the specified formats described before.

#### Data Pre-Processing

The data pre-processing consists of getting the fluorescent intensities measured for each marker for each cell and transform them into discrete states.

To run the pre-processing, use the following code snippet:

``` r
cell_data_processed <- process_cell_data(cell_data,
                                         channel_data,
                                         sample_data,
                                         sampleID_col = "Sample_ID",
                                         n_threads = 8
                                         )
```

where `sampleID_col` is the respective column in your `cell_data` containing the sample IDs matching the ones described in `sample_data`. Set the `n_threads` parameter to the desired number of threads to be used.

The marker states will be encoded as follows:

| State | Encoding |
|:-----:|:--------:|
|  \-   |    0     |
|  \+   |    1     |
|  ++   |    2     |
|  ...  |    ..    |

#### Combinatorial Phenotypes Cell Counting

This step will take the pre-processed Cell Data and output a data.frame with all possible combinations of markers and marker states counting the respective number of cells that contain that phenotype for each sample.

``` r
comb_phenotypes <- combinatorial_phenotype_counts(cell_data_processed,
                                                  min_count = 10,
                                                  sample_fraction_min_counts = 0.5,
                                                  n_threads = 8
                                                  )
```

where `min_count` is the minimum number of cells a phenotype must have, and `sample_fraction_min_counts` is the fraction of the samples that must meet the `min_count` criterion. All phenotypes that do not match these minimum requirements will be discarded.

The resulting `data.frame` will have its first N columns with the encoding respective to the N markers considered and the next M columns respective to the M samples considered containing the cell counts.

#### Statistical Testing

PhenoComb offers, for now, three statistical tests to assess the statistical relevance of each phenotype. They are:

-   Group Comparison using Mann-Whitney U test.
-   Correlation using Kendall's Rank Correlation.
-   Survival test using Cox Proportional Hazards Model (CPHR)

The following code snippet shows an example of a group comparison test:

``` r
relevant_phenotypes <- compute_statistically_relevant_phenotypes(comb_phenotypes,
                                                                   channel_data,
                                                                   sample_data,
                                                                   test_type = "group",
                                                                   groups_column = "Group",
                                                                   g1 = "g1",
                                                                   g2 = "g2",
                                                                   max_pval = 0.05,
                                                                   n_threads = 8
)
```

where `test_type` indicates a group comparison test, the `groups_column` is the respective column in `sample_data` where your grouping information is stored, and `g1` and `g2` are the values in `groups_column` to be considered for each group.

The resulting `data.frame` will have its first N columns with the encoding respective to the N markers considered, the next M columns respective to the M samples considered containing the cell frequencies, and the last columns contain the results from the respective statistical test.

#### Independent Phenotypes Filtering

PhenoComb provides a tool to help identifying the most relevant phenotypes that are independent on each other. For example, a phenotype Marker1+Marker2+ is dependent on Marker1+, thus this tool will help identifying the most relevant ones.

The following code snippet shows an example of how use this tool:

``` r
final_phenotypes <- get_independent_relevant_phenotypes(relevant_phenotypes,
                                                        channel_data,
                                                        n_phenotypes = 1000,
                                                        min_confidence = 0.5,
                                                        n_threads = 8
)
```

where `n_phenotypes` is the maximum number of phenotypes considered in the analysis filtering fow lowest associated p-value. The `min_confidence` is a measurement of how frequent the phenotype shows up during the analysis steps.

The resulting `data.frame` will have its first column with the phenotype name generated from the encoding, the next one to three columns will have the results from the statistical test depending on which one was chosen in the previous step, the confidence obtained by the algorithm, and the next M columns are respective to the M samples considered containing the cell frequencies used in the statistical test.
