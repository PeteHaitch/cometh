#' Analyse, manage and visualise co-methylation data.
#'
#' \pkg{cometh} provides tools for analysing, managing and visualising co-methylation 
#' data. Loosely speaking, co-methylation is the correlation structure of DNA 
#' methylation.
#'
#' Please refer to the vignettes to see how to use the  \pkg{cometh} package.
#'
#' @docType package
#' @name cometh-package
#' @useDynLib cometh
#' @import Biobase S4Vectors IRanges GenomeInfoDb GenomicRanges Rcpp   
#' @import BiocParallel
#' @import data.table
#' @import R.utils
#' @references TODO
NULL

## TODO: Refine @import to @importFrom, e.g. I only need fread from data.table, 
## gunzip from R.utils, and probably only a subset of functions from 
## BiocParallel.