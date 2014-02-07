#' The cometh class
setClass("CoMeth", contains = "SummarizedExperiment")

# TODO: Require seqlengths or seqinfo when constructing a CoMeth
# TODO: Make documentation more like GRanges, i.e. ?CoMeth describes the class and there is a subsection titled "Construction" that described the constructor, etc.
#' The constructor function for cometh objects
#'
#' Construct a \code{CoMeth} object. Users will generally use the \code{\link{read.comethylation}} to construct a \code{CoMeth} object. 
#' @param seqnames Rle object, character vector, or factor containing the sequence names.
#' @param pos A numeric matrix storing the positions of each m-tuple as a row. Thus nrow(pos) = the number of m-tuples and ncol(pos) = m.
#' @param counts A numeric matrix storing the number of times each co-methylation pattern is observed at each m-tuple. Thus nrow(pos) = the number of m-tuples and ncol(pos) = the number of possible co-methylation patterns = 2 ^ m.
#' @param m An integer storing the size of the m-tuples, i.e. the "m" in m-tuple.
#' @param methylation_type A character vector storing the type of methylation loci for these m-tuples. Possible values are "CG", "CHG", "CHH" or "CNN". Multiple values can be specified, e.g. c("CG", "CHG").
#' @param sample_name The name of the sample.
#' @param strand An (optional) Rle object, character vector, or factor containing the strand information.
#' @param seqlengths The sequence lengths of the reference genome of the sample. Must be an integer vector named with the sequence names and containing the lengths (or NA) for each level(seqnames).
#' @param seqinfo An (optional) \code{\link[GenomicRanges]{Seqinfo}} object containing information about the reference genome of the sample
#'
#' @export
#' @seealso \code{\link{read.comethylation}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
CoMeth <- function(seqnames, pos, counts, m, methylation_type, sample_name, strand = "*", seqlengths = NULL, seqinfo = NULL){
  if (missing(seqnames) || all(!is.character(seqnames), !is(seqnames, 'Rle'), !is.factor(seqnames))){
    stop("Need 'seqnames'. Must be an Rle object, a character vector or a factor containing the sequence names")
  }
  if (missing(pos) || !is.matrix(pos) || !is.numeric(pos)){
    stop("Need 'pos'. Must be a numeric matrix.")
  }
  if (missing(counts) || !is.matrix(counts) || !is.numeric(counts)){
    stop("Need 'counts'. Must be a numeric matrix.'")
  }
  if (missing(m) || !is.integer(m) || length(m) > 1 || m < 1){
    stop("'m' must be specified. Must be a single int and must be greater than 0")
  }
  if (missing(methylation_type) || !is.character(methylation_type) || !all(methylation_type %in% c('CG', 'CHG', 'CHH', 'CNN'))){
    stop("'methylation_type' must be specified and must be a character vector. 'CG', 'CHG', 'CHH' and 'CNN', or a vector of some combination of these, are the only valid 'methylation_type'")
  }
  if (!all(methylation_type %in% c('CG', 'CHG', 'CHH', 'CNN'))){
    stop("'methylation_type' must be a character vector")
  }
  if (ncol(pos) != m){
    stop("ncol(pos) should be equal to m")
  }
  if (ncol(counts) != (2 ^ m)){
    stop("ncol(counts) should be equal to 2^m")
  }
  if (nrow(counts) != nrow(pos)){
    stop("nrow(pos) should be equal to nrow(counts)")
  }
  if (is.null(colnames(pos)) || ! all(grepl('pos[0-9]', colnames(pos)))){
    stop("Column names of 'pos' must be: pos1, pos2, ..., posm")
  }
  if (is.null(colnames(counts)) || ! all(colnames(counts) == make_m_tuple_names(m))){
    stop_msg <- paste0("Column names of 'counts' must be: ", paste0(make_m_tuple_names(m), collapse = ', '))
    stop(stop_msg)
  }
  message("There is minimal checking of the 'pos' and 'counts' matrices. E.g. no check is made that each row of 'pos' is sorted; no check is made that all entries of 'counts' and 'pos' are positive integers")

  gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = pos[, 1], end = pos[, ncol(pos)]), strand = strand, seqlengths = seqlengths, seqinfo = seqinfo)
  # Store each column of counts as a separate assay
  counts_colnames <- colnames(counts)
  counts <- lapply(seq_len(2 ^ m), function(i, counts){
    counts[, i, drop = FALSE]
  }, counts = counts)
  names(counts) <- counts_colnames
  # Need to store the other positions if m > 2. Each additional position is stored as a separate assay
  if (m > 2){
    add_pos <- lapply(seq(2, m - 1, 1), function(i, pos){
      pos[, i, drop = FALSE]
    }, pos = pos)
    names(add_pos) <- paste0('pos', seq(2, m - 1, 1))
  } else {
    add_pos <- NULL
  }
  assays <- SimpleList(c(counts, add_pos))
  colData <- DataFrame(sample_name = sample_name, m = m, methylation_type = paste0(sort(methylation_type), collapse = '/'))
  cometh <- SummarizedExperiment(assays = assays, rowData = gr, colData = colData)
  cometh <- as(cometh, "CoMeth")

  return(cometh)
}

setMethod("show", signature(object = "CoMeth"), function(object) {
            cat("An object of type 'CoMeth' with:\n")
            cat(paste0(" ", nrow(object), " ", colData(object)$methylation_type[1], " ", colData(object)$m[1], "-tuples\n"))
            cat(paste0(" ", ncol(object), " samples\n"))
          })

setReplaceMethod("sampleNames",
                 signature = signature(
                   object = "CoMeth",
                   value = "ANY"),
                 function(object, value) {
                   colnames(object) <- value
                   object
                 })

# Set the validity of a CoMeth object
setValidity("CoMeth", function(object) {
  msg <- validMsg(NULL, .checkAssayNames(object, make_m_tuple_names(getM(object))))
  if(class(rowData(object)) != "GRanges")
    msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
  ## Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if(min(getCounts(object)) < 0) {
    msg <- validMsg(msg, "'counts' has negative entries")
  }
  if(min(getPos(object)) < 0) {
    msg <- validMsg(msg, "'pos' has negative entries")
  }
  # TODO: Check that each row of getPos(object) is sorted. The naive implementation "any(apply(X = getPos(object), FUN = function(x){is.unsorted(x)}, MARGIN = 1))" is very slow

  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

setMethod("sampleNames", "CoMeth", function(object) {
  colData(object)[, 'sample_name']
})

#' Get the \code{counts} from a \code{\link{CoMeth}} object
#'
#' Get the counts from a \code{\link{CoMeth}} object. \code{counts} is a numeric matrix storing the number of times each co-methylation pattern is observed at each m-tuple. 
#'
#' @param CoMeth a \code{\link{CoMeth}} object
#'
#' @return A numeric matrix
getCounts <- function(CoMeth) {
  if (ncol(CoMeth) > 1){
    stop('Can only get counts from a CoMeth object that contains data on a single sample')
  }
  stopifnot(is(CoMeth, "CoMeth"))

  cc <- grep('[MU]', assayNames(CoMeth))
  counts <- sapply(cc, function(cc, CoMeth){assay(CoMeth, cc)}, CoMeth = CoMeth)
  colnames(counts) <- grep('[MU]', assayNames(CoMeth), value = TRUE)

  return(counts)
}

#' Get the \code{coverage} from a \code{\link{CoMeth}} object.
#'
#' The \code{coverage} of a \code{\link{CoMeth}} object is the number of reads for each m-tuple, i.e. \code{rowSums(getCounts(CoMeth))}
#' @param CoMeth A \code{\link{CoMeth}} object
#'
#' @return A numeric vector
#' @export
getCoverage <- function(CoMeth) {
  if (ncol(CoMeth) > 1){
    stop('Can only get counts from a CoMeth object that contains data on a single sample')
  }
  stopifnot(is(CoMeth, "CoMeth"))

  rowSums(getCounts(CoMeth))
}

#' Get \code{m} from a \code{\link{CoMeth}} object
#'
#' Get the size of m-tuples, i.e. the \code{m} in m-tuples, from a \code{\link{CoMeth}} object.
#' @param CoMeth A CoMeth object
#'
#' @return An integer
#' @export
getM <- function(CoMeth) {
  if (ncol(CoMeth) > 1){
    stop('Can only get counts from a CoMeth object that contains data on a single sample')
  }
  stopifnot(is(CoMeth, "CoMeth"))

  m <- colData(CoMeth)[, 'm']

  return(m)
}

#' Obtain the positions from a CoMeth object
#' 
#' The \code{pos} of a \code{\link{CoMeth}} object are the positions of the cytosines that make up each m-tuple.
#' @param CoMeth A CoMeth object
#'
#' @return A matrix
#' @export
getPos <- function(CoMeth) {
  if (ncol(CoMeth) > 1){
    stop('Can only get counts from a CoMeth object that contains data on a single sample')
  }
  stopifnot(is(CoMeth, "CoMeth"))

  m <- getM(CoMeth)
  if (m == 1){
    pos  <- matrix(start(CoMeth), ncol = 1)
  } else if (m == 2) {
    pos <- cbind(start(CoMeth), end(CoMeth))
  } else {
    pc <- grep('pos', assayNames(CoMeth))
    pos <- cbind(start(CoMeth), sapply(pc, function(pc, CoMeth){assay(CoMeth, pc)}, CoMeth = CoMeth), end(CoMeth))
  }
  colnames(pos) <- paste0('pos', seq_len(m))

  return(pos)
}

# TODO: cbind, rbind, combine. Require that 'm', 'seqinfo', 'methylation_type' are identical when rbind-ing cometh objects