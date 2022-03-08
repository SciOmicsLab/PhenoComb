
#' PhenoComb
#' 
#' Package for Phenotype Combinatorics Analysis
#' 
#' @docType package
#' @author Paulo Burke <pauloepburke@gmail.com>
#' @import data.table
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib PhenoComb
#' @name PhenoComb
NULL  

.onLoad <- function(libname, pkgname){
  
  
  op <- options()
  op.phenocomb <- list(
    PhenoComb.verbose = TRUE,
    PhenoComb.log.file = NULL
  )
  
  toset <- !(names(op.phenocomb) %in% names(op))
  if(any(toset)) options(op.phenocomb[toset])
  
  invisible()
}


.onUnload <- function(libpath){
  library.dynam.unload("PhenoComb",libpath)
}