### =========================================================================
### All classes 
### =========================================================================

## TODO: Should I define prototypes in class definitions?
## TODO: Define MTuplesList; check split(MTuples) produces an MTuplesList.

### -------------------------------------------------------------------------
### MTuples 
###

#' @rdname MTuples
#' @export
setClass('MTuples',
         representation(extraPos = "matrix"),
         contains = "GRanges",
)

setValidity("MTuples", function(object) {
  
  # Include all .valid.MTuples.* functions in this vector
  msg <- c(.valid.MTuples.pos(object))
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.MTuples.pos <- function(object) {
  # There are effectively no checks made if is.na(m) == TRUE, because this
  # corresponds to the "empty" MTUples object
  
  msg <- NULL
    
  # The "empty" MTuples object
  if (isTRUE(all(is.na(object@ranges@start)))) {
    m <- NA_integer_
  } else if (isTRUE(all(is.na(object@extraPos)))) {
    # m = 1 or 2
    if (isTRUE(all(object@ranges@start == 
                     (object@ranges@start + object@ranges@width - 1L)))) {
      m <- 1L
    } else{
      m <- 2L
      if (any(object@ranges@start >=
            (object@ranges@start + object@ranges@width - 1L))) {
        msg <- validMsg(msg, paste0("positions within each 2-tuple must be ", 
                                    "sorted in strictly increasing order, ", 
                                    "i.e., ", sQuote('pos1'), " < ", 
                                    sQuote('pos2'), '.'))
      }
    }
  } else{
    # m >= 3
    m <- ncol(object@extraPos) + 2L
    if (any(object@ranges@start >= object@extraPos[, 1]) |
          any(object@extraPos[, m - 2L] >= 
                (object@ranges@start + object@ranges@width - 1L)) |
          !.allRowsSortedCpp(object@extraPos)) {
      msg <- validMsg(msg, paste0("positions within each ", m, 
                                  "-tuple must be ", "sorted in strictly ", 
                                  "increasing order, ", 
                                  "i.e., ", sQuote('pos1'), " < ... < ", 
                                  sQuote(paste0('pos', m)), '.'))
    }
  }
  
  # Since we already have checked that pos is sorted, we need only check that
  # the first element, i.e. start, is positive.
  # We also need only check that the extraPos are integers, because the start 
  # and end are guaranteed to be integers since they are part of an IRanges.
  # And since extraPos only really exists if m >= 3, this check is only 
  # necessary if that condition is satisfied.
  if (!is.na(m)){
    if (m < 3) {
      if (any(object@ranges@start) < 0) {
        msg <- validMsg(msg, paste0("positions within each m-tuple must be ", 
                                    "positive integers."))
      }
    } else{
      if (any(object@ranges@start) < 0 | !is.integer(object@extraPos)) {
        msg <- validMsg(msg, paste0("positions within each m-tuple must be ", 
                                    "positive integers."))
      }
    }
  }
  
  return(msg)
}

#' MTuples class and constructor
#' 
#' @description
#' An m-tuple is a tuple of individual basepairs that are on the same 
#' chromosome, where \eqn{m = 1, 2, \ldots} is the number of positions in the 
#' tuple. For example, (chr1:30, chr1:33, chr1:40) is a 3-tuple of the positions 
#' 30, 33 and 40 on chromosome 1. Note that the strand of the m-tuple is optional. 
#' The \code{MTuples} class extends the \code{\link[GenomicRanges]{GRanges}} 
#' class by adding the \code{extraPos} slot (see below for details). 
#'  
#' @param seqnames \code{\link[IRanges]{Rle}} object, character vector, or 
#' factor containing the sequence names.
#' @param pos \code{matrix} where each row stores the ordered positions of an
#' m-tuple.
#' @param strand \code{\link[IRanges]{Rle}} object, character vector, or factor 
#' containing the strand information.
#' @param ... Optional metadata columns. These columns cannot be named "start", 
#' "end", "width", or "element". A named integer vector "\code{seqlength}" can 
#' be used instead of \code{seqinfo}.
#' @param seqlengths an integer vector named with the sequence names and 
#' containing the lengths (or \code{NA}) for each \code{level(seqnames)}.
#' @param seqinfo a \code{\link[GenomeInfoDb]{Seqinfo}} object containing 
#' allowed sequence names and lengths (or \code{NA}) for each 
#' \code{level(seqnames)}.
#' 
#' @section Differences to \code{\link[GenomicRanges]{GRanges}}:
#' There are several differences between 
#' \code{\link[GenomicRanges]{GRanges}} and \code{MTuples}.
#' 
#' \subsection{Conceptual}{
#' The main difference is suggested in the names; a 
#' \code{\link[GenomicRanges]{GRanges}} object is designed to store 
#' \emph{ranges} whereas an \code{MTuples} object is designed to store tuples. 
#' A range is a closed interval and so includes all points in between the 
#' "start" and "end", inclusive. By contrast, a tuple is a \emph{set} and so 
#' only includes those points listed in the tuple. See the section below called 
#' \code{findOverlaps} for some implications of this difference.}
#' 
#' \subsection{\code{extraPos}}{
#' An \code{MTuples} object contains an extra slot called \code{extraPos}. This 
#' is only relevant if creating m-tuple with \eqn{m \ge 3} when it stores 
#' positions 2, ..., m - 1 of the m-tuple as an integer matrix. When \eqn{m < 3} 
#' this slot is a 1-column matrix of \code{NA}.}
#' 
#' @details
#' Internally, an m-tuple is stored as a \code{\link[GenomicRanges]{GRanges}} 
#' with any "extra" positions stored in the \code{extraPos} slot. For example, 
#' the first and last positions of the  3-tuple (chr1:30, chr1:33, chr1:40) are 
#' stored as the \code{start} and \code{end}, respectively, of the 
#' \code{\link[GenomicRanges]{GRanges}}, i.e., 
#' \code{GRanges('chr1', IRanges(start = 30, end = 40))}. The "extra" position, 
#' chr1:33, is stored as an element in matrix in the \code{extraPos} slot.
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
#' @return An \code{MTuples} object.
#' 
#' @aliases MTuples-class
#' 
#' @examples
#' cat("TODO")
#' 
#' @include AllGenerics.R
#' @rdname MTuples
#' @export
MTuples <- function(seqnames = Rle(), pos = matrix(), 
                    strand = Rle("*", length(seqnames)), ..., 
                    seqlengths = NULL, seqinfo = NULL) {
  
  ## TODO: Any other argument checks from the old constructor
  ## Check that all required arguments are of the correct type
  ## Or try to coerce to the correct type
  if (!is.matrix(pos)){
    stop(sQuote('pos'), ' must be a matrix.') 
  }
  if (all(!is.na(pos))) {
    # min(pos) < 0 is faster than any(pos < 0)
    if (min(pos) < 0) {
      stop('Some values in ', sQuote('pos'), ' are negative.')
    }
    if (!is.integer(pos)) {
      if (isTRUE(any(round(pos) != pos))) {
        stop(paste0('Some values in ', sQuote('pos'), ' are not integers.'))
      }
      mode(pos) <- "integer"
    }
  }

  ## Other argument checks are deferred to the GRanges constructor, which
  ## is called further down.
  m <- ncol(pos)
  
  if (!isTRUE(all(is.na(pos)))){
    if (m == 1L){
      ranges <- IRanges(start = pos[, 1L], width = 1L)
      extraPos <- matrix(NA_integer_, nrow = length(ranges))
    } else if (m == 2L){
      ## When m = 2, need extra check that start != end, because this isn't 
      ## caught by the validity methods when a single m-tuple is passed.
      ## Otherwise this supposed 2-tuple is silently converted to a 1-tuple.
      ## See GitHub issue #8 (https://github.com/PeteHaitch/cometh/issues/8)
      if (!.allRowsSortedCpp(pos)){
        stop(paste0("positions within each 2-tuple must be sorted in strictly increasing order, i.e. ", sQuote('pos1'), " < ", sQuote('pos2'), '.'))
      }
      ranges <- IRanges(start = pos[, 1L], end = pos[, 2L])
      extraPos <- matrix(NA_integer_, nrow = length(ranges))
    } else{
      ranges <- IRanges(start = pos[, 1L], end = pos[, m])
      extraPos <- pos[, seq(from = 2L, to = m - 1, by = 1L), drop = FALSE]
    }
  } else{
    ranges <- IRanges()
    extraPos <- matrix()
  }
  
  new("MTuples", GRanges(seqnames = seqnames, ranges = ranges, strand = strand, seqinfo = seqinfo, ...), extraPos = extraPos) 
}

### -------------------------------------------------------------------------
### CoMeth 
###

#' @rdname CoMeth
#' @export
setClass('CoMeth', 
         contains = c("VIRTUAL", "SummarizedExperiment"))

setValidity("CoMeth", function(object) { 
  
  # First need to check that rowData is an MTuples object.
  # Otherwise some of the .valid.CoMeth.* functions won't work
  msg <- .valid.CoMeth.rowData(object)
  if (is.null(msg)){
    
    # Include all other .valid.CoMeth.* functions in this vector
    msg <- c(.valid.CoMeth.methylation_type(object), 
             .valid.CoMeth.noDuplicates(object),
             .valid.CoMeth.assayNames(object))
  }
  
  # Can't run this check unless the assayNames are correct
  if (is.null(msg)){
    msg <- .valid.CoMeth.counts(object)
  }
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.CoMeth.rowData <- function(object) {
  
  msg <- NULL
  
  if (class(object@rowData) != "MTuples"){
    msg <- validMsg(msg, paste0(sQuote('rowData(CoMeth)'), " must be an ", 
                                sQuote('MTuples'), " object."))
  }
  return(msg)
}

.valid.CoMeth.methylation_type <- function(object) {
  
  msg <- NULL
  
  if (sum(grepl("^methylation_type$", colnames(object@colData))) != 1L){
    msg <- validMsg(msg, paste0(sQuote('colData'), " of ", sQuote('CoMeth'), 
                                " must contain column ", 
                                sQuote('methylation_type'), 
                                " once and only once."))
  }
  
  # Can only run the next check if colData contains the 'methylation_type' 
  # column
  if (is.null(msg)){
    
    mts <- object@colData[, grepl("^methylation_type$", 
                                  colnames(object@colData))]
    
    if (!all(sapply(X = mts, FUN = .valid_methylation_type))){
      msg <- validMsg(msg,
                      paste0(sQuote('methylation_type'), 
                             " for each sample must be ", sQuote('CG'), ", ",
                             sQuote('CHG'), ", ", sQuote('CHH'), " or ", 
                             sQuote('CNN'), ", or some combination of these, ", 
                             "e.g., ", sQuote("CG/CHG"), 
                             ".\nCombinations must sorted alphabetically and ", 
                             "be separated by a forward slash (", 
                             sQuote('/'), ")."))
    }
  }
  return(msg)
}

.valid.CoMeth.assayNames <- function(object) {
  
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
  
  # Check that object has 'MM..M', ..., 'UU..U' assays names
  if (!all(c(.make_m_tuple_names(m)) %in% names(object@assays$data))){
    msg <- validMsg(msg, paste0("assay names must include: ", 
                                paste0(sQuote(c(.make_m_tuple_names(m))), 
                                       collapse = ", "), "."))
  }
}

## TODO: Might remove this because the duplicated method for large MTuples is
## freaking slow.
.valid.CoMeth.noDuplicates <- function(object) {
  
  msg <- NULL
  
  # Check that there are no duplicate m-tuples
  if (any(duplicated(object))){
    msg <- validMsg(msg, paste0(sQuote('CoMeth'), 
                                " object cannot contain duplicate m-tuples."))
  }
  
  return(msg)
}

.valid.CoMeth.counts <- function(object) {
  
  msg <- NULL
  
  m <- getM(rowData(object))
  assay_names <- .make_m_tuple_names(m)
  
  # Check that all 'counts' are non-negative
  # Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than 
  # any(assay() < 0) if it is false
  if (min(sapply(object@assays$data[assay_names], min, na.rm = TRUE), 
          na.rm = TRUE) < 0) {
    msg <- validMsg(msg, paste0(sQuote('counts'), 
                                " cannot have negative entries."))
  }
  
  return(msg)
}

## TODO: Profile constructor.
## TODO: Proper way to link to papers; DOI? URL?
## TODO: Re-write the CoMeth constructor and include here. Remove the old CoMeth
## constructor from methods-CoMeth-class.R.

#' CoMeth class and constructor
#' 
#' @description
#' The \code{CoMeth} class is a virtual class extended from 
#' \code{\link[GenomicRanges]{SummarizedExperiment}}. The subclasses, 
#' \code{CoMeth1} (for 1-tuples), \code{CoMeth2} (for 2-tuples) and 
#' \code{CoMeth3Plus} (for m-tuples, \eqn{m \ge 3}), are containers for 
#' storing methylation patterns at methylation loci m-tuples.
#' 
#' Most users will use the \code{\link{read.comethylation}} function to 
#' construct the appropriate \code{CoMeth} object from the output \code{.tsv} 
#' file(s) produced by the \code{comethylation} Python software.
#' 
#' @param assays A \code{\link[IRanges]{SimpleList}} of matrix elements. All 
#' elements of the list must have the same dimensions, and dimension names 
#' (if present) must be consistent across elements and with the row names of 
#' \code{rowData} and \code{colData}. See below for the required assays for 
#' each subclass of \code{CoMeth}.
#' @param rowData A \code{\link{MTuples}} instance describing the m-tuples of 
#' interest. Row names, if present, become the row names of the \code{CoMeth}. 
#' The length of the \code{CoMeth} must equal the number of rows of the 
#' matrices in assays.
#' @param colData A \code{\link[IRanges]{DataFrame}} describing the samples. 
#' Row names, which are required, become the column names of the \code{CoMeth}. 
#' See below for additional details.
#' @param exptData An optional \code{\link[IRanges]{SimpleList}} of arbitrary 
#' content describing the overall experiment.
#' @param ... S4 methods \code{list} and \code{matrix}, arguments identical to 
#' those of the \code{SimpleList} method.
#' @param verbose A \code{logical(1)} indicating whether messages about data 
#' coercion during construction should be printed.
#' 
#' @section Differences to \code{\link[GenomicRanges]{SummarizedExperiment}}:
#' There are several differences between 
#' \code{\link[GenomicRanges]{SummarizedExperiment}} and \code{CoMeth}.
#' 
#' \subsection{\code{rowData}}{
#' A \code{CoMeth} object uses an \code{\link{MTuples}} object in the 
#' \code{rowData} slot rather than a \code{\link{GRanges}} or 
#' \code{GenomicRangesList} object.
#' 
#' Furthermore, the \code{MTuples} object must include 
#' \code{\link[GenomeInfoDb]{Seqinfo}} of the reference genome. This 
#' "reference genome" can accommodate a spiked-in unmethylated genome (normally 
#' lambda phage) that is commonly used as a control in bisulfite-sequencing 
#' experiments.}
#' 
#' \subsection{\code{assays}}{
#' The assays of a \code{CoMeth} object must include the counts of how many 
#' times each methylation pattern is observed for that m-tuple in each sample. 
#' For example, there are \eqn{2^1} possible methylation patterns at 1-tuples, 
#' namely \code{M} and \code{U}; there are \eqn{2^2 = 4} possible methylation 
#' patterns at 2-tuples, namely \code{MM}, \code{MU}, \code{UM} and \code{UU}; 
#' there are \eqn{2^3} possible patterns at 3-tuples, namely \code{MMM}, ..., 
#' \code{'UUU'}. The \code{CoMeth} class enforces non-negative integer values 
#' in these "counts" assays. The assays must be a \code{\link[IRanges]{SimpleList}}.
#' 
#' There are is an additional assay for each of the \code{CoMeth1} and 
#' \code{CoMeth2} subclasses; the assays of a \code{CoMeth1} object also 
#' include the \eqn{\beta-{\text{values}}}, \code{beta}, defined as 
#' \eqn{\beta = \frac{M}{M + U}}; the assays of a \code{CoMeth2} object also 
#' include the within-fragment comethylation, which is measured as a log-odds 
#' ratio, \code{LOR}, and is defined as 
#' \eqn{\log(\frac{(MM + 0.5)(UU + 0.5)}{(MU + 0.5)(UM + 0.5)})}. The 
#' \code{CoMeth3Plus} contains no additional assays.}
#' 
#' \subsection{\code{colData}}{
#' Unlike a \code{\link[GenomicRanges]{SummarizedExperiment}}, a \code{CoMeth} 
#' object requires the specification of \code{colData}. The row names of this 
#' \code{colData}, which are required, must be the sample names and will be 
#' used as the column names of the \code{CoMeth}. The \code{colData} must 
#' include the methylation type as a column of the \code{colData} and use the 
#' column name \code{methylation_type}. The \code{methylation_type} must be one 
#' of \code{CG}, \code{CHG}, \code{CHH} or \code{CNN}, or some combination of 
#' these separated by a '/', e.g., \code{CG/CHG}.}
#' 
#' @section Coercion:
#' \strong{TODO}: Insert help for any coercion methods.
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
#' @section For developers:
#' Developers may wish to extend the \code{CoMeth} virtual class to other 
#' specific m-tuples, e.g. \code{CoMeth7} to include an additional assay 
#' that is specific to 7-tuples.
#' 
#' @return The \code{CoMeth} constructor returns the appropriate subclass 
#' depending on  the size of the m-tuples, i.e., the value of \eqn{m}.
#'
#' @aliases CoMeth-class CoMeth1-class CoMeth2-class CoMeth3Plus-class
#'  
#' @references See \url{http://www.github.com/PeteHaitch/comethylation} for the 
#' \code{comethylation} Python software.
#' 
#' @seealso \code{\link{read.comethylation}} for a function to read in the 
#' \code{comethylation} output \code{.tsv} file(s) and construct the 
#' appropriate \code{CoMeth} object.
#' @seealso \code{\link[GenomicRanges]{SummarizedExperiment}} for the class that
#' \code{CoMeth} extends.
#' 
#' @examples
#' cat("TODO")
#' @include AllGenerics.R
#' @rdname CoMeth
#' @export
CoMeth <- function(assays = SimpleList(), rowData = MTuples(), 
                   colData = DataFrame(), exptData = SimpleList(), ..., 
                   verbose = FALSE) {
  
  ## Argument checks: assay names and contents, colData names and contents,
  ## seqinfo in MTuples, rowData is MTuples
  
  m <- getM(rowData)
  
  if (is.na(m)) {
    stop("Missing ", sQuote('rowData'), '.')
  } else if (m == 1L){
    class <- "CoMeth1"
  } else if (m == 2L){
    class <- "CoMeth2"
  } else{
    class <- "CoMeth3Plus"
  }
  new(class, SummarizedExperiment(assays = assays, rowData = rowData, 
                                  colData = colData, exptData = exptData, 
                                  verbose = verbose))
}

### -------------------------------------------------------------------------
### CoMeth1
###
### A concrete subclass of CoMeth for storing methylation patterns at 1-tuples.
### CoMeth1 should have 'M', 'U' and 'beta' as assays.

setClass("CoMeth1", 
         contains = "CoMeth")

## TODO: Should .calid.CoMeth() be included in the validity check?
setValidity("CoMeth1", function(object) {
  
  # Include all other .valid.CoMeth1.* functions in this vector
  msg <- c(.valid.CoMeth1.counts(object)) 
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.CoMeth1.counts <- function(object) {
  
  msg <- NULL
  
  # Check that contains 1-tuples
  if (getM(rowData(object)) != 1L){
    msg <- validMsg(msg, paste0("Expected 1-tuples in a ", sQuote("CoMeth1"), 
                                " object."))
  } 

  # Check assay names 
  # M, and U are already checked by validity method for CoMeth 
  # VIRTUAL class
  assay_names <- names(object@assays$data@listData)
  extra_assay_names <- assay_names[-which(assay_names %in% c('M', 'U'))] 
  if (!extra_assay_names %in% "beta"){
    msg <- validMsg(msg, 
                    paste0(sQuote("CoMeth1"), 
                           " object must include assay: ", 
                           sQuote("beta")))
  }
    
  return(msg)
}

## TODO: Is this a bottleneck?
### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html
setAs("CoMeth1", "SummarizedExperiment",
      def = function(from) {
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value) { 
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
### CoMeth2 should have 'MM', 'MU', 'UM', 'UU' and 'LOR' as assays.

setClass("CoMeth2", 
         contains = "CoMeth")

## TODO: Should .calid.CoMeth() be included in the validity check?
setValidity("CoMeth2", function(object) {
  
  # Include all other .valid.CoMeth2.* functions in this vector
  msg <- c(.valid.CoMeth2.counts(object)) 
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.CoMeth2.counts <- function(object) {
  
  msg <- NULL
  
  # Check that contains 2-tuples
  if (getM(rowData(object)) != 2L){
    msg <- validMsg(msg, paste0("Expected 2-tuples in a ", sQuote("CoMeth2"), 
                                " object."))
  } 
  
  # Check assay names 
  # MM, MU, UM, and UU are already checked by validity method for CoMeth 
  # VIRTUAL class
  assay_names <- names(object@assays$data@listData)
  extra_assay_names <- assay_names[-which(assay_names %in% 
                                            c('MM', 'MU', 'UM', 'UU'))] 
  if (!extra_assay_names %in% "LOR"){
    msg <- validMsg(msg, 
                    paste0(sQuote("CoMeth2"), 
                           " object must include assay: ", 
                           sQuote("LOR")))
  }
    
  return(msg)
}

## TODO: Is this a bottleneck?
### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### (Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R)
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html
setAs("CoMeth2", "SummarizedExperiment",
      def = function(from) {
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value) { 
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
### Unlike CoMeth1 and CoMeth2, CoMeth3Plus only has the default assays of the 
### CoMeth VIRTUAL class.
setClass("CoMeth3Plus", 
         contains = "CoMeth")

## TODO: Should .calid.CoMeth() be included in the validity check?
setValidity("CoMeth3Plus", function(object) {
  
  # Include all other .valid.CoMeth3Plus.* functions in this vector
  msg <- c(.valid.CoMeth3Plus.counts(object)) 
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.CoMeth3Plus.counts <- function(object) {
  
  
  msg <- NULL
  
  # Check that contains m-tuples, with m >= 3
  if (getM(rowData(object)) < 3L){
    msg <- validMsg(msg, paste0("Expected m-tuples (m >= 3) in a ", 
                                sQuote("CoMeth3Plus"), " object."))
  } 
  
  # Compulsory assay names already checked by validity method for CoMeth 
  # VIRTUAL class
    
  return(msg)
}

## TODO: Is this a bottleneck?
### Coercion:
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.
### (Source: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/VariantAnnotation/R/AllClasses.R)
### See also https://stat.ethz.ch/pipermail/r-devel/2012-October/065028.html
setAs("CoMeth3Plus", "SummarizedExperiment",
      def = function(from) {
        if (strict) {
          force(from)
          class(from) <- "SummarizedExperiment"
        }
        from
      },
      replace = function(from, to, value) {
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
setClass("MethylationLociSet", 
         contains = "GRanges",
         slots = c(methylation_type = "character"))

setValidity("MethylationLociSet", function(object) {
  
  # Include all other .valid.MethylationLociSet* functions in this vector
  msg <- c(.valid.MethylationLociSet.methylation_type(object)) 
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
  })

.valid.MethylationLociSet.methylation_type <- function(object) {
  
  msg <- NULL
  
  if (!all(sapply(X = object@methylation_type, .valid_methylation_type))){
    msg <- validMsg(msg, paste0("Invalid ", sQuote("methylation_type")))
  }
  return(msg)
}


### -------------------------------------------------------------------------
### BetaCor 
###
setClass("BetaCor", 
         contains = "DataFrame",
         slots = c(methylation_type = "character", NIL = "integer", 
                   IPD = "integer", feature_name = "character", 
                   same_feature = "logical", ignore_strand = "logical", 
                   seqinfo = "Seqinfo", method = "character"))

## TODO: BetaCor validity
setValidity("BetaCor", function(object){
  TRUE
  })