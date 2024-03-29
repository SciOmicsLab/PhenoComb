% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/independent_phenotype_identification_server.R
\name{get_independent_relevant_phenotypes_server}
\alias{get_independent_relevant_phenotypes_server}
\title{Identify independent relevant phenotypes (Server Version)}
\usage{
get_independent_relevant_phenotypes_server(
  output_folder,
  channel_file,
  sample_file,
  input_significant_phenotypes = "significant_phenotypes.csv",
  output_file = "independent_phenotypes.csv",
  log_file = "independent_phenotypes.log",
  n_phenotypes = 1000,
  min_confidence = 0.5,
  max_pval = NULL,
  n_threads = 1,
  verbose = FALSE
)
}
\arguments{
\item{output_folder}{Path to folder where output files from this and previous steps should be saved.}

\item{channel_file}{Path to a ".csv" file containing columns named: Channel, Marker, T1, [T2, T3, ... , Tn], [OOB].}

\item{sample_file}{Path to a ".csv" file containing a Sample_ID column and additional grouping columns for the samples.}

\item{n_phenotypes}{maximum number of phenotypes to be considered from \code{phen_data} filtered by lowest p-values. Default: 1000.}

\item{min_confidence}{Minimal confidence threshold to filter output. Default: 0.5.
#' @param parent_phen Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.}

\item{max_pval}{Maximum p-value. Used to filter phenotypes in final output.}

\item{n_threads}{Number of threads to be used. Default: 1.}

\item{verbose}{If TRUE, print outputs from log to stdout.}

\item{parent_phen}{Parent phenotype to filter for. All phenotypes in the output will contain the parent phenotype.}
}
\value{
Output is saved to file.
}
\description{
It uses a network clustering process to identify most representative
phenotypes when several similar phenotypes are relevant.
}
\details{
This server version is just a convenience to find correct inputs and
outputs results to file. No performance difference.

Parameters \code{groups_column, g1, g2, max_pval,} and \code{parent_phen} should
match the parameter used in \code{statistically_relevant_phenotypes_server} in order
to find correct input files.
}
