### =========================================================================
### All classes 
### =========================================================================

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
    if (!.allRowsSortedCpp(pos)){
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
### CoMeth is a VIRTUAL class with concrete subclasses CoMeth1 (for 1-tuples),
### CoMeth2 (for 2-tuples) and CoMeth3Plus (for m-tuples, m >= 3).

.valid.CoMeth.rowData <- function(object){
  
  msg <- NULL
  
  if (class(object@rowData) != "MTuples"){
    msg <- validMsg(msg, paste0(sQuote('rowData(CoMeth)'), " must be an ", sQuote('MTuples'), " object."))
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

.valid.CoMeth.assayNames <- function(object){
  
  msg <- NULL
  
  # Get m
    if (isTRUE(all(is.na(object@rowData@extraPos)))){
      # m = 1 or 2
      if (isTRUE(all(object@rowData@ranges@start == 
                       (object@rowData@ranges@start 
                        + object@rowData@ranges@width - 1L)))){
        m <- 1L
      } else{
        m <- 2L
      }
    } else{
      # m > 2
      m <- ncol(object@rowData@extraPos) + 2L
    }
  
  # Check that object has 'MM..M', ..., 'UU..U', and 'EP' assays names
  if (!all(c(.make_m_tuple_names(m), "EP") %in% names(object@assays$data))){
    msg <- validMsg(msg, paste0("assay names must include: ", 
                                paste0(sQuote(c(.make_m_tuple_names(m), "EP")), 
                                       collapse = ", "), "."))
  }
}

.valid.CoMeth.noDuplicates <- function(object){
  
  msg <- NULL
  
  ## Check that there are no duplicate m-tuples
  if (any(duplicated(object))){
    msg <- validMsg(msg, paste0(sQuote('CoMeth'), " object cannot contain duplicate m-tuples."))
  }
  
  return(msg)
}

.valid.CoMeth.counts <- function(object){
  
  msg <- NULL
  
  m <- getM(rowData(object))
  assay_names <- .make_m_tuple_names(m)
  
  ## Check that all 'counts' are non-negative
  ## Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if (min(sapply(object@assays$data[assay_names], min, na.rm = TRUE), na.rm = TRUE) < 0) {
    msg <- validMsg(msg, paste0(sQuote('counts'), " cannot have negative entries."))
  }
  
  return(msg)
}

.valid.CoMeth <- function(object){
  
  # First need to check that rowData is an MTuples object.
  # Otherwise some of the .valid.CoMeth.* functions won't work
  msg <- .valid.CoMeth.rowData(object)
  if (is.null(msg)){
    msg <- c(.valid.CoMeth.methylation_type(object), 
             .valid.CoMeth.noDuplicates(object),
             .valid.CoMeth.assayNames(object)) # Include all other .valid.CoMeth.* functions in this vector
  }
  ## Can't run this check unless the assayNames are correct
  if (is.null(msg)){
    msg <- .valid.CoMeth.counts(object)
  }
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

## TODO: Decide whether to export the class definition; see vignette(topic = 'namespace', package = 'roxygen2').
#' An S4 class to store co-methylation patterns at m-tuples of genomic 
#' positions.
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
         contains = c("VIRTUAL", "SummarizedExperiment"),
         validity = .valid.CoMeth)


#### TODO: Define concrete subclasses of CoMeth
# The current ones (below) aren't quite correct.
# CoMeth1: 1-tuples
# CoMeth2: 2-tuples
# CoMeth3Plus: m-tuples, m >= 3

# CoMeth should have 'MM..M', ..., 'UU..U' and 'EP' assays (not enforced)
# CoMeth1 should have 'M', 'U', 'EP' and 'beta' as an assay 
# CoMeth2 should also have 'LOR' as an assay.
# CoMeth3Plus doesn't need any extra assays.
# User could extend CoMeth VIRTUAL class, for specific m-tuples,
# e.g. CoMeth7 for 7-tuples, which might include an additional (and, as yet, 
# unknown) assay that is specific to 7-tuples.

# TODO: Move EP to CoMeth VIRTUAL class

# What about zeta? Where should it go? Well, it depends. 
# zeta is the average level of methylation for all reads that overlap all 
# methylation loci in the m-tuple. __This is generally not equivalent to the 
# average beta values of the m-tuple__.
# __Because of this, I will not include this as part of the CoMeth VIRTUAL class
# or any of its concrete subclasses__

### -------------------------------------------------------------------------
### CoMeth1
###
### A concrete subclass of CoMeth for storing methylation patterns at 1-tuples.

.valid.CoMeth1.counts <- function(object){
  
  msg <- NULL
  
  ## Check that contains 1-tuples
  if (getM(rowData(object)) != 1L){
    msg <- validMsg(msg, paste0("Expected 1-tuples in a ", sQuote("CoMeth1"), " object."))
  } 

  ## Check assay names 
  ## M, U and EP are already checked by validity method for CoMeth 
  ## VIRTUAL class
  assay_names <- names(object@assays$data@listData)
  extra_assay_names <- assay_names[-which(assay_names %in% c('M', 'U', 'EP'))] 
  if (!extra_assay_names %in% "beta"){
    msg <- validMsg(msg, 
                    paste0(sQuote("CoMeth1"), 
                           " object must include assay: ", 
                           sQuote("beta")))
  }
    
  return(msg)
}

.valid.CoMeth1 <- function(object){
  
  msg <- c(.valid.CoMeth(object), .valid.CoMeth1.counts(object)) # Include all other .valid.CoMeth1.* functions in this vector
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

setClass("CoMeth1", 
         contains = "CoMeth",
         validity = .valid.CoMeth1)

### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### (Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R)
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html

setAs("CoMeth1", "SummarizedExperiment",
      def = function(from){
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value)
      {
        firstTime <- TRUE
        for (nm in slotNames(value)) {
          v <- slot(value, nm)
          if (firstTime) {
            slot(from, nm, FALSE) <- v
            firstTime <- FALSE
          } else {
            `slot<-`(from, nm, FALSE, v)
          }
        }
        from
      }
)

### -------------------------------------------------------------------------
### CoMeth2
###
### A concrete subclass of CoMeth for storing methylation patterns at 2-tuples.

.valid.CoMeth2.counts <- function(object){
  
  msg <- NULL
  
  ## Check that contains 2-tuples
  if (getM(rowData(object)) != 2L){
    msg <- validMsg(msg, paste0("Expected 2-tuples in a ", sQuote("CoMeth2"), " object."))
  } 
  
  ## Check assay names 
  ## MM, MU, UM, UU and EP are already checked by validity method for CoMeth 
  ## VIRTUAL class
  assay_names <- names(object@assays$data@listData)
  extra_assay_names <- assay_names[-which(assay_names %in% 
                                            c('MM', 'MU', 'UM', 'UU', 'EP'))] 
  if (!extra_assay_names %in% "LOR"){
    msg <- validMsg(msg, 
                    paste0(sQuote("CoMeth2"), 
                           " object must include assay: ", 
                           sQuote("LOR")))
  }
    
  return(msg)
}

.valid.CoMeth2 <- function(object){
  
  msg <- c(.valid.CoMeth(object), 
           .valid.CoMeth2.counts(object)) # Include all other .valid.CoMeth2.* functions in this vector
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

setClass("CoMeth2", 
         contains = "CoMeth",
         validity = .valid.CoMeth2)

### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### (Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R)
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html

setAs("CoMeth2", "SummarizedExperiment",
      def = function(from){
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value)
      {
        firstTime <- TRUE
        for (nm in slotNames(value)) {
          v <- slot(value, nm)
          if (firstTime) {
            slot(from, nm, FALSE) <- v
            firstTime <- FALSE
          } else {
            `slot<-`(from, nm, FALSE, v)
          }
        }
        from
      }
)

### -------------------------------------------------------------------------
### CoMeth3Plus
###
### A concrete subclass of CoMeth for storing methylation patterns at m-tuples,
### when m >= 3.

.valid.CoMeth3Plus.counts <- function(object){
  
  msg <- NULL
  
  ## Check that contains m-tuples, with m >= 3
  if (getM(rowData(object)) < 3L){
    msg <- validMsg(msg, paste0("Expected m-tuples (m >= 3) in a ", 
                                sQuote("CoMeth3Plus"), " object."))
  } 
  
  ## Compulsory assay names already checked by validity method for CoMeth 
  ## VIRTUAL class
    
  return(msg)
}

.valid.CoMeth3Plus <- function(object){
  
  msg <- c(.valid.CoMeth(object), 
           .valid.CoMeth3Plus.counts(object)) # Include all other .valid.CoMeth3Plus.* functions in this vector
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

setClass("CoMeth3Plus", 
         contains = "CoMeth",
         validity = .valid.CoMeth3Plus)

### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### (Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R)
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html

setAs("CoMeth3Plus", "SummarizedExperiment",
      def = function(from){
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value)
      {
        firstTime <- TRUE
        for (nm in slotNames(value)) {
          v <- slot(value, nm)
          if (firstTime) {
            slot(from, nm, FALSE) <- v
            firstTime <- FALSE
          } else {
            `slot<-`(from, nm, FALSE, v)
          }
        }
        from
      }
)

### -------------------------------------------------------------------------
### MethylationLociSet 
###

.valid.MethylationLociSet.methylation_type <- function(object){
  
  msg <- NULL
  
  if (!object@methylation_type %in% c('CG', 'CHG', 'CHH', 'CNN', 'CG/CHG', 'CG/CHH', 'CG/CNN', 'CHG/CHH', 'CHG/CNN', 'CHH/CNN', 'CG/CHG/CHH', 'CG/CHG/CNN', 'CHG/CHH/CNN', 'CG/CHG/CHH/CNN')){
    msg <- validMsg(msg, paste0("Invalid ", sQuote("methylation_type")))
  }
  return(msg)
}

.valid.MethylationLociSet <- function(object){
  
  msg <- c(.valid.MethylationLociSet.methylation_type(object)) # Include all other .valid.MethylationLociSet* functions in this vector
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

setClass("MethylationLociSet", 
         contains = "GRanges",
         slots = c(methylation_type = "character"),
         validity = .valid.MethylationLociSet)