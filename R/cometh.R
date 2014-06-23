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
#' @importFrom BiocParallel bplapply
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @references TODO
NULL