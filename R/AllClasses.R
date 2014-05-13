### =========================================================================
### All classes 
### =========================================================================

## TODO: Remove use of class-specific methods

### -------------------------------------------------------------------------
### MTuples 
###

.valid.MTuples.pos <- function(object){
  
  msg <- NULL
    
  # The "empty" MTuples object
  if (isTRUE(all(is.na(object@ranges@start)))){
    m <- NA_integer_
  } else if (isTRUE(all(is.na(object@extraPos)))){
    # m = 1 or 2
    if (isTRUE(all(object@ranges@start == (object@ranges@start + object@ranges@width - 1L)))){
      m <- 1L
      pos <- as.matrix(object@ranges@start)
    } else{
      m <- 2L
      pos <- cbind(object@ranges@start, object@ranges@start + object@ranges@width - 1L)
    }
  } else{
    # m > 2
    m <- ncol(object@extraPos) + 2L
    pos <- cbind(object@ranges@start, object@extraPos, object@ranges@start + object@ranges@width - 1L)
  }
  
  if (!is.na(m)){
    if (!allRowsSortedCpp(pos)){
      msg <- validMsg(msg, paste0("positions in each m-tuple must be sorted in strictly increasing order, i.e. ", sQuote('pos1'), " < ", sQuote('pos2'), " < ", sQuote('...'), " < ", sQuote('posm')))
    }
  }
  
  if (!is.na(m)){
    if (min(pos) < 0){  # min(x) < 0 is faster than any(x < 0)
      msg <- validMsg(msg, paste0("positions in each m-tuple must be positive integers."))
    }
  }
  
  return(msg)
}

.valid.MTuples <- function(object){
  msg <- c(.valid.MTuples.pos(object)) # Include all .valid.MTuples.* functions in this vector
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

## TODO: Decide whether to export the class definition; see vignette(topic = 'namespace', package = 'roxygen2').
## TODO: Currently need to do '?"MTuples-class"' to find help; would prefer '?MTuples'
#' An S4 class to represent m-tuples of genomic positions.
#' 
#' @details
#' The \code{MTuples} class extends the \code{\link[GenomicRanges]{GRanges}} class by adding 
#' the \code{extraPos} slot (see below for details). An m-tuple is a tuple of 
#' individual basepairs that are on the same chromosome, where 'm' is the number
#' of positions in the tuple. For example, (chr1:30, chr1:33, chr1:40) is 
#' a 3-tuple of the positions on chromosome 1. Note the strand of the m-tuple 
#' is optional.
#' 
#' Internally, this example 3-tuple is stored as a GRanges object with the first 
#' and last positions of the m-tuple stored as the \code{start} and \code{end} 
#' of the GRanges interval, respectively. That is,
#' \code{GRanges('chr1', IRanges(start = 30, end = 40))}. The "extra" position, 
#' chr1:33, is stored in the \code{extraPos} matrix.
#'
#' @slot extraPos A numeric matrix storing "extra" positions in m-tuples, 
#' provided m >= 3. If m = 1 or 2, \code{extraPos} is a matrix of \code{NA}. 
#' The \code{extraPos} matrix has as many rows as there are m-tuples in the 
#' \code{MTuples} object.
#' 
#' @section Constructor:
#' \strong{TODO}: Insert help for constructor method.
#' 
#' @section Coercion:
#' \strong{TODO}: Insert help for any coerction methods.
#' 
#' @section Accessors:
#' \strong{TODO}: Insert help for any accessor methods.
#' 
#' @section Splitting and combining:
#' \strong{TODO}: Insert help for any splitting and combining methods.
#' 
#' @section Subsetting:
#' \strong{TODO}: Describe any subsetting methods.
#'
#' @section Filtering:
#' \strong{TODO}: Describe any filtering methods.
#' 
#' @section Methods based on findOverlaps:
#' \strong{TODO} Insert help for any findOverlaps-based methods.
#' 
#' @section Other methods:
#' \strong{TODO}: Describe any other methods.
#' 
#' @include AllGenerics.R
setClass('MTuples',
         representation(extraPos = "matrix"),
         contains = "GRanges",
         validity = .valid.MTuples
         )

### -------------------------------------------------------------------------
### CoMeth 
###

.valid.CoMeth.counts <- function(object){
  
  msg <- NULL
  
  ## Check assay names
  if (isTRUE(all(is.na(object@ranges@start)))){
    m <- NA_integer_
  } else if (isTRUE(all(is.na(object@extraPos)))){
    # m = 1 or 2
    if (isTRUE(all(object@ranges@start == (object@ranges@start + object@ranges@width - 1L)))){
      m <- 1L
    } else{
      m <- 2L
    }
  } else{
    # m > 2
    m <- ncol(object@extraPos) + 2L
  }
  if (!identical(names(object@assays$data@listData), .make_m_tuple_names(m))){
    msg <- validMsg(msg, paste0("assay names must be: ", paste0("assay names must be: ", paste0(sQuote(.make_m_tuple_names(m)), collapse = ', ')), "."))
  }
  
  ## Check that all 'counts' are non-negative
  ## Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if (min(sapply(object@assays$data, min, na.rm = TRUE), na.rm = TRUE) < 0) {
    msg <- validMsg(msg, paste0(sQuote('counts'), " cannot have negative entries."))
  }  
  return(msg)
}

.valid.CoMeth.methylation_type <- function(object){
  
  msg <- NULL
  
  if (sum(grepl("^methylation_type$", colnames(object@colData))) != 1L){
    msg <- validMsg(msg, paste0(sQuote('colData'), " of ", sQuote('CoMeth'), " must contain column ", sQuote('methylation_type'), " once and only once."))
  }
  
  ## Can only run the next check if colData contains the 'methylation_type' column
  if (is.null(msg)){
    if (!all(object@colData[, grepl("^methylation_type$", colnames(object@colData))] %in% c('CG', 'CHG', 'CHH', 'CNN', 'CG/CHG', 'CG/CHH', 'CG/CNN', 'CHG/CHH', 'CHG/CNN', 'CHH/CNN', 'CG/CHG/CHH', 'CG/CHG/CNN', 'CHG/CHH/CNN', 'CG/CHG/CHH/CNN'))){
      msg <- validMsg(msg, paste0(sQuote('methylation_type'), " for each sample must be ", sQuote('CG'), ", ", sQuote('CHG'), ", ", sQuote('CHH'), " or ", sQuote('CNN'), ", or some combination of these, e.g., ", sQuote("CG/CHG"), ".\nCombinations must sorted alphabetically and be separated by a forward slash (", sQuote('/'), ")."))
    }
  }
  
  return(msg)
}

.valid.CoMeth.rowData <- function(object){
  
  msg <- NULL
  
  if (class(object@rowData) != "MTuples"){
    msg <- validMsg(msg, paste0(sQuote('rowData(CoMeth)'), " must be an ", sQuote('MTuples'), " object."))
  }
  return(msg)
}

.valid.CoMeth.noDuplicates <- function(object){
  
  msg <- NULL
  
  ## Check that there are no duplicate m-tuples
  if (any(duplicated(object))){
    msg <- validMsg(msg, paste0(sQuote('CoMeth'), " object cannot contain duplicate m-tuples."))
  }
  return(msg)
}

.valid.CoMeth <- function(object){
  
  # First need to check that rowData is an MTuples object.
  # Otherwise some of the .valid.CoMeth.* functions won't work
  msg <- .valid.CoMeth.rowData(object)
  if (is.null(msg)){
    msg <- c(.valid.CoMeth.counts(object), .valid.CoMeth.methylation_type(object), .valid.CoMeth.noDuplicates(object)) # Include all other .valid.CoMeth.* functions in this vector
  }
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

## TODO: Decide whether to export the class definition; see vignette(topic = 'namespace', package = 'roxygen2').
#' An S4 class to represent co-methylation patterns at m-tuples of genomic positions.
#' 
#' @details
#' The \code{CoMeth} class is based on the
#' \code{\link[GenomicRanges]{SummarizedExperiment}} class. The main difference
#' is that rather than using a \code{\link[GenomicRanges]{GRanges}} object as 
#' the \code{rowData}, a \code{CoMeth} object uses an \code{\link{MTuples}} 
#' object.
#' 
#' The assays of a \link{CoMeth} object are the counts of how many times each 
#' co-methylation pattern is observed for that m-tuple in each sample. 
#' For example, the possible co-methylation patterns of 2-tuples are 'MM', 'MU', 
#' 'UM' and 'UU' and thus there are four assays of the same names.
#' 
#' @section Constructor:
#' \strong{TODO}: Insert help for constructor method.
#' 
#' @section Coercion:
#' \strong{TODO}: Insert help for any coerction methods.
#' 
#' @section Accessors:
#' \strong{TODO}: Insert help for any accessor methods.
#' 
#' @section Splitting and combining:
#' \strong{TODO}: Insert help for any splitting and combining methods.
#' 
#' @section Subsetting:
#' \strong{TODO}: Describe any subsetting methods.
#'
#' @section Filtering:
#' \strong{TODO}: Describe any filtering methods.
#' 
#' @section Methods based on findOverlaps:
#' \strong{TODO} Insert help for any findOverlaps-based methods.
#' 
#' @section Other methods:
#' \strong{TODO}: Describe any other methods.
#' 
#' @include AllGenerics.R
setClass('CoMeth', 
         contains = "SummarizedExperiment",
         validity = .valid.CoMeth)
