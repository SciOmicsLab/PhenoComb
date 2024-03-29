# PhenoComb

Burke PEP, Strange A, Monk E, Thompson B, Amato CM, Woods DM. PhenoComb: a discovery tool to assess complex phenotypes in high-dimensional single-cell datasets. Bioinform Adv. 2022 Aug 3;2(1):vbac052. doi: 10.1093/bioadv/vbac052. PMID: 36699375; PMCID: PMC9710698.

An R package for combinatorial analysis of phenotypes.

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
-   Server: storing all outputs into files, aimed at processing large datasets and number of markers.

Regardless of the desired workflow, the inputs are the same. However, for the server version, the inputs must be accessed directly from files rather than R objects.

### Input Preparation

PhenoComb requires three inputs:

1.  Cell Data (data.frame, .fcs, or .csv)
2.  Channel Data (data.frame or .csv)
3.  Sample Data (data.frame or .csv)

#### Cell Data

Cell Data can be your concatenated FCS file with the channels as columns and fluorescence values as rows for each cell. It must also contain one column with the sample identification. As a default, we suggest a column called "Sample_ID" (not to be confounded with SampleID standard column from FlowJo, **never use that**), but can be called as whatever the user desire.

You can find a step-by-step tutorial on how to prepare your fcs data using FlowJo [here](docs/FlowJo_tutorial/FlowJo_Tutorial.md).

Example:

| SampleID   | FSC-A  | SSC-A  | FJComp-BL1-A | FJComp-BL3-A | FJComp-RL1-A | FJComp-RL2-A | FJComp-VL1-A | Sample_ID |
|------------|--------|--------|--------------|--------------|--------------|--------------|--------------|-----------|
| 54695.4688 | 345236 | 76813  | 1743.154663  | -415.9430847 | 667.958557   | 5997.42432   | -26.10379    | 1         |
| 54883.9648 | 394485 | 89886  | 1945.750854  | -2512.404541 | 2948.14404   | 14496.3232   | -491.70593   | 1         |
| 54184.4531 | 338380 | 60629  | 1059.707153  | -966.7955933 | 3331.49146   | 7877.6333    | 14.0653687   | 2         |
| 54737.4766 | 306087 | 114762 | 1153.990112  | 726.7556152  | 142.62999    | 487.799774   | -59.632195   | 2         |
| 53959.1289 | 615779 | 145288 | 3343.247559  | 1712.121582  | 353.986755   | 16573.8262   | -1022.7864   | 3         |
| 53953.7852 | 158319 | 269510 | 1944.918457  | 851.4907227  | 149.089676   | 680.195435   | 676.559326   | 3         |

#### Channel Data

Channel Data is a table containing the columns described below. Required columns are marked in bold.

-   **Channel**: channel name respective to its column in Cell Data. Names must match.
-   **Marker**: name of the marker measured in the respective Channel. Please avoid spaces and special characters.
-   **T1**: column with a threshold value used to discretize measured values into states (-/+).
-   T2, T3, ..., Tn: additional thresholds to discretize values into more than two states (-/++/+++/...).
-   OOB: aka Out Of Bounds values, used to discard cell measurement if measured value is above this threshold.

Channels not specified in this file will be discarded.

Example:

| Channel      | Marker | T1   | T2    | OOB    |
|--------------|--------|------|-------|--------|
| FJComp-BL1-A | CD95   | 1133 |       | 63661  |
| FJComp-BL3-A | CD45RO | 4249 |       | 104384 |
| FJComp-RL1-A | CD127  | 1316 |       | 16279  |
| FJComp-RL2-A | CD44   | 4249 | 25183 | 227456 |
| FJComp-VL1-A | KLRG1  | 814  |       | 55306  |

#### Sample Data

Sample Data is a table containing the names of your samples and meta information used to perform further statistical tests.

The required column is marked in bold.

-   **Sample_ID**: column with sample IDs contained in the respective column in Cell Data. The following columns are only examples. Can be named or represent data as the user convenience.
-   Sample_Group: group label respective to the sample.
-   Correlated_Measurement: a value for the phenotypes be correlated to.
-   ...

Example:

| Sample_ID | Treatment |
|-----------|-----------|
| 1         | None      |
| 2         | None      |
| 3         | None      |
| 4         | IL10      |
| 5         | IL10      |
| 6         | IL10      |

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
library(PhenoComb)
library(flowCore)

# Use this line if your input is a FCS file (Must have flowCore installed).
cell_data <- as.data.frame(flowCore::read.FCS("/path/to/cell_data.fcs",truncate_max_range = FALSE)@exprs)

# Use this line if your input is a CSV file.
cell_data <- read.csv("/path/to/cell_data.csv")

channel_data <- read.csv("/path/to/channel_data.csv")
sample_data <- read.csv("/path/to/sample_data.csv")
```

All input files must be in accordance to the specified formats described before.

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

where `test_type` indicates a group comparison test, the `groups_column` is the respective column in `sample_data` where your grouping information is stored, and `g1` and `g2` are the values in `groups_column` to be considered for each group. The `max_pval` parameter will filter the outputs considering it as a maximum p-value.

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

### Server Workflow

The server workflow follows the same steps of the local workflow, except that the pre-processing step is incorporated into the Combinatorial Phenotypes Cell Counting step. It will save all the outputs from each step, and also log files, to the desired output folder. It means that the size of the dataset you can analyze is not limited by your computer memory, but by time and disk space for the outputs.

#### Folder Structure

Although is not obligatory, we recommend the following folder structure to perform your analyses.

``` bash
DataSetName
├── OriginalData
├── inputs
│   ├── cell_data.fcs
│   ├── channel_data.csv
│   └── sample_data.csv
├── outputs
└── scripts
    └── phenoComb_script.R
```

> Please be aware that, for high number of markers (26 markers and up), the size of the output files can be on the TeraByte scale. Make sure you will have enough disk space.


#### Combinatorial Phenotypes Cell Counting

The following code snippet is considering the folder structure previously recommended. It will perform the pre-processing and the combinatorial phenotypes cell counting. It will generate two output files ("combinatorial_phenotype_counts.csv" and "combinatorial_phenotypes.log") in the `output_folder` indicated.

``` r
# File: phenoComb_script.R

library(PhenoComb)

# Run the next two lines if you are running this script from RStudio
# to change your working directory to where this file is. Comment/remove them
# if you are running using Rscript.
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

# Process raw data and generates all combinatorial phenotypes
combinatorial_phenotype_counts_server(cell_file = "../inputs/cell_data.fcs",
                                      channel_file = "../inputs/channel_data.csv",
                                      sample_file = "../inputs/sample_data.csv",
                                      output_folder = "../outputs",
                                      min_count = 10,
                                      sample_fraction_min_counts = 0.5,
                                      sampleID_col = "Sample_ID",
                                      verbose = TRUE,
                                      n_threads = 8
)
```

The verbose option will also print the logging to your default prompt. The `sampleID_col` parameter is the respective column in your `cell_file` containing the sample IDs matching the ones described in `sample_file`. The `min_count` option is the minimum number of cells a phenotype must have, and `sample_fraction_min_counts` is the fraction of the samples that must meet the `min_count` criterion. All phenotypes that do not match these minimum requirements will be discarded. Set the `n_threads` parameter to the desired number of threads to be used.


#### Statistical Testing

Similarly to the local workflow, the statistical testing will evaluate all phenotypes computed in the last step accordingly to the specified statistical test. However, it will read the input from a file and output the results to another file. The following code snippet is an example of performing a group comparison statistical test:

``` r
statistically_relevant_phenotypes_server(output_folder = "../outputs",
                                         channel_file = "../inputs/channel_data.csv",
                                         sample_file = "../inputs/sample_data.csv",
                                         input_phenotype_counts = "combinatorial_phenotype_counts.csv",
                                         input_phenotype_counts_log = "combinatorial_phenotypes.log",
                                         output_file = "significant_phenotypes.csv",
                                         log_file = "significant_phenotypes.log",
                                         groups_column = "Group",
                                         g1 = "g1",
                                         g2 = "g2",
                                         max_pval = 0.05,
                                         verbose = TRUE,
                                         n_threads = 8
)
```

where `test_type` indicates a group comparison test, the `groups_column` is the respective column in `sample_data` where your grouping information is stored, and `g1` and `g2` are the values in `groups_column` to be considered for each group. The `max_pval` parameter will filter the outputs considering it as a maximum p-value.

If `input_phenotype_counts` and `input_phenotype_counts_log` are not provided, it will automatically search for them in the `output_folder`. However, if explicitly provided, it will make some of the processing faster. Consider using it when analyzing big datasets.

The `output_file` and `log_file` parameters ar optional, but it's a convenience when performing different group comparisons on the same dataset, for example. They will be stored into `output_folder`, **DO NOT** put the full path to a file here, only the file name.


#### Independent Phenotypes Filtering

Differently from the previous two steps, this tool **is** limited by your machine's memory. To avoid any memory issues, choose a reasonable number of phenotypes to be considered using the `n_phenotypes` option. A reasonable number is 5000 to 10000 phenotypes. They will be selected from the input based on the smallest associated p-values. It will automatically detect the statistical test used in the previous step.

The following code snippet will perform the Independent Phenotypes Filtering:

``` r
get_independent_relevant_phenotypes_server(output_folder = "../outputs",
                                           channel_file = "../inputs/channel_data.csv",
                                           sample_file = "../inputs/sample_data.csv",
                                           input_significant_phenotypes = "significant_phenotypes.csv",
                                           output_file = "independent_phenotypes.csv",
                                           log_file = "independent_phenotypes.log",
                                           n_phenotypes = 1000,
                                           min_confidence = 0.5,
                                           verbose = TRUE,
                                           n_threads = 8
)
```

The `output_file` and `log_file` parameters ar optional, but it's a convenience when trying different number of phenotypes on the same dataset, for example. They will be stored into `output_folder`, **DO NOT** put the full path to a file here, only the file name.



