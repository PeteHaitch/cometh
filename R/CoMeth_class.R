#### TODOs ####
# TODO: Make documentation more like GRanges, i.e. ?CoMeth describes the class and there is a subsection titled "Construction" that described the constructor, etc.
# TODO: cbind, rbind, combine, comparse, %<=% (and similar). Require that 'm', 'seqinfo', 'methylation_type' are identical when rbind-ing CoMeth objects. Allow for multiple samples in a CoMeth object
# TODO: Function collapseStrand(). Should only work with symmetric methylation_type, e.g. CG, CHG. Will be hard to do for more complication methylation_type, e.g. CG/CHG, without using the reference genome (slow). Add the 'collapse_strand' option to CoMeth constructor.
# TODO: Perhaps simplify the description of the 'pos' and 'counts' parameters to note that these can be list of data.frame or data.frame-like objects (DataFrame) since a data.frame is just a list of lists.
# TODO: Check for any orphan TODOs

#### CoMeth class definition, constructor and validity function ####
#' The CoMeth class
.CoMeth <- setClass("CoMeth", representation(extra_pos = "list"), contains = "SummarizedExperiment")

#' The constructor function for CoMeth objects
#'
#' Users will generally use the \code{\link{read.comethylation}} to initialise a \code{CoMeth} object, however, this is the general method to initialise a \code{CoMeth} object. The initialisation allows for multiple samples in the same object by passing most arguments as a list where each element of the list corresponds to the arguments for a given sample. These multiple samples will then be combined into a single \code{CoMeth} object.
#' @param sample_names A character vector containing the names of the samples.
#' @param seqnames A list of character vectors, a list of factors or an \code{\link[IRanges]{RleList}} object containing the sequence names. Each element of the list should be named and match one of the \code{sample_names}.
#' @param pos A list of lists. Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the \code{sample_names}. The number of elements of each sub-list should be equal to \code{m}. For a given sample, each element of the inner-list stores the positions of each m-tuple as an integer vector. E.g. \code{pos[[2]][[1]]} contains for \code{sample2}, all \code{pos1} for each m-tuple as an integer vector.
#' @param counts Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the \code{sample_names}. The number of elements of each sub-list should be equal to \eqn{2 ^ \code{m}}. For a given sample, each element of the inner-list stores the number of times each co-methylation pattern is observed at each m-tuple. E.g. \code{pos[[2]][[1]]} contains for \code{sample2}, how many times the first co-methylation pattern was observed for each m-tuple as an integer vector.
#' @param strand An (optional) list of character vectors, list of factors or an \code{\link[IRanges]{RleList}} object containing the strand information of each m-tuple. \strong{WARNING}: If 'strand = NULL' all m-tuples in the resulting \code{\link{CoMeth}} object will have their strand set to '*' to signify that the strand is unknown or irrelevant (such as when methylation measurements have been combined across strands). \strong{WARNING}: m-tuples will not be combined across samples if they are on different strands. 
#' @param m An integer storing the size of the m-tuples, i.e. the \code{m} in m-tuple. Only a single value is accepted, and not a list, because \code{m} must be the same for all samples in a \code{CoMeth} object.
#' @param methylation_type A character vector storing the type of methylation loci for these m-tuples. Possible values are "CG", "CHG", "CHH" or "CNN". Multiple values can be specified, e.g. c("CG", "CHG"). Only a single value is accepted, and not a list, because \code{methylation_type} must be the same for all samples in a \code{CoMeth} object.
#' @param seqinfo A \code{\link[GenomicRanges]{Seqinfo}} object containing information about the reference genome of the samples. All samples must be mapped against the same reference genome, however, multiple reference genomes per sample are allowed. This is to allow for a spike-in unmethylated genome (normally lambda phage), which is a common step in a bisulfite-sequencing protocol.
#' @param sort_cometh logical. Should the cometh object be sorted by genomic co-ordinates? Note that sorting is based soley on the coordinates of the first and last methylation loci in each m-tuple. Thus, if m > 2 the order of m-tuples on the same chromosome with the same start and end co-ordinates is random, e.g. there is no guarantee on the order of chr1:(1, 8, 10) and chr1:(1, 4, 10).
#'
#' @export
#' @note The sub-lists of 'pos' and 'counts' can in fact be data.frames since a data.frame is just a list of lists
#' @seealso \code{\link{read.comethylation}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
CoMeth <- function(sample_names, seqnames, pos, counts, strand = NULL, m, methylation_type, seqinfo, sort_cometh = TRUE){
  if (missing(sample_names) || !identical(unique(sample_names), sample_names)){
    stop("Need 'sample_names'. Must be a character vector where each element is a unique sample name.")
  }
  
  # Must be checked before 'counts' is checked because it uses make_m_tuple_names(m)
  if (missing(m) || !is.integer(m) || length(m) > 1 || m < 1){
    stop("'m' must be specified and must be a single, positive integer.")
  }
  
  if (missing(seqnames) || is.null(seqnames) || length(seqnames) != length(sample_names) || !all(names(seqnames) %in% sample_names) || all( !all(sapply(seqnames, is.character)), !all(sapply(seqnames, is.factor)), !all(sapply(seqnames, FUN = function(x){is(x, 'Rle')})))){
    stop("Need 'seqnames'. Must be a list of character vectors, a list of factors or an RleList object containing the sequence names. Each element of the list should be named and match one of the 'sample_names'.")
  }
  
  if (missing(pos) || !all(names(pos) %in% sample_names) || length(pos) != length(sample_names) || !all(sapply(pos, is.list)) || !all(sapply(pos, function(x){sapply(x, is.integer)})) || !all(apply(sapply(pos, function(x){sapply(x, length)}), FUN = zero_range, MARGIN = 2)) || !all(sapply(pos, length) == m)){
    stop("Need 'pos'. Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the 'sample_names'. The number of elements of each sub-list should be equal to 'm'. For a given sample, each element of the inner-list stores the positions of each m-tuple as an integer vector. E.g. pos[[2]][[1]] contains for 'sample2', all 'pos1' for each m-tuple as an integer vector.")
  }
  if (any(is.null(sapply(pos, names))) || !all(grepl('pos[0-9]', sapply(pos, names)))){
    stop("Names of each sub-list of 'pos' must be: ", paste0('pos', seq_len(m), collapse = ', '))
  }
  
  if (missing(counts) || !all(names(counts) %in% sample_names) || length(counts) != length(sample_names) || !all(sapply(counts, is.list)) || !all(sapply(counts, function(x){
    sapply(x, is.integer)
    })) || !all(apply(sapply(counts, function(x){sapply(x, length)}), FUN = zero_range, MARGIN = 2)) || !all(sapply(counts, length) == 2 ^ m)){
    stop("Need 'counts'. Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the 'sample_names'. The number of elements of each sub-list should be equal to '2^m'. For a given sample, each element of the inner-list stores the number of times each co-methylation pattern is observed at each m-tuple. E.g. pos[[2]][[1]] contains for 'sample2', how many times the first co-methylation pattern was observed for each m-tuple as an integer vector.")
  }
  if (any(sapply(sapply(counts, names), is.null)) || !all(apply(X = sapply(counts, names), FUN = function(x, m){identical(sort(x), make_m_tuple_names(m))}, MARGIN = 2, m = m))){
    stop("Names of each sub-list of 'counts' must be: ", paste0(make_m_tuple_names(m), collapse = ', '))
  }
  
  if (!all(sapply(sample_names, FUN = function(sample_name, counts, pos){
    zero_range(c(sapply(counts[[sample_name]], length), sapply(pos[[sample_name]], length)))
    }, counts = counts, pos = pos))){
    stop("For each element of 'sample_names', the lengths of all sub-lists for its corresponding element of 'pos' should equal the lengths of sub-lists for its corresponding element of 'counts'.")
  }
  
  if (!is.null(strand) || !all(sapply(strand, function(x){x}) == '*')){
    if (!all(names(strand) %in% sample_names) || all( !all(sapply(strand, is.character)), !all(sapply(strand, is.factor)), !all(sapply(strand, FUN = function(x){is(x, 'Rle')}))) || !all(sapply(sample_names, FUN = function(sample_name, strand, seqnames){
      length(strand[[sample_name]] == length(seqnames[[sample_name]]))
      }, strand = strand, seqnames = seqnames))){
      stop("If 'strand' is not NULL then it must be a list of character vectors, a list of factors or an RleList object containing the strand of each m-tuple. Each element of the list should be named and match one of the 'sample_names' and have the length equal to the number of m-tuples for that sample.")
    } else if (("*" %in% as.vector(unlist(strand, use.names = FALSE))) && (("+" %in% as.vector(unlist(strand, use.names = FALSE))) || ("-" %in% as.vector(unlist(strand, use.names = FALSE))))){
      warning("'strand' contains '*' as well as at least one of '+' or '-'. m-tuples will not be combined across samples if they are on different strands. This means that m-tupes with the exact same genomic co-ordinates across samples, but with different strands (e.g. '*' and '+'), will not be combined.\nTo combine all m-tuples based on genomic co-ordinates regardless of strand, set 'strand = NULL'.")
    }
  } else{
    strand <- RleList(lapply(seqnames, function(x){rep('*', length(x))}))
  }
  
  if (missing(methylation_type) || !is.character(methylation_type) || !all(methylation_type %in% c('CG', 'CHG', 'CHH', 'CNN'))){
    stop("'methylation_type' must be specified. 'methylation_type' must be 'CG', 'CHG', 'CHH' and 'CNN', or a vector of some combination of these, e.g. methylation_type = c('CG', 'CHG')")
  }

  if (missing(seqinfo)){
    stop("'seqinfo' must be specified.")
  }
  
  if (!isTRUEorFALSE(sort_cometh)){
    stop("'sort_cometh' must be TRUE or FALSE")
  }
  
  warning("There is minimal checking of the 'pos' and 'counts' matrices. E.g. no check is made that each row of the matrices in 'pos' is sorted; no check is made that all entries of the matrices in 'counts' and all the entries of the matrices in 'pos' are positive integers")
  
  # Combine the list-wise data into matrix-like data whilst taking care of common and sample-specific m-tuple, i.e. filling in zeros when a sample doesn't have any observations for that m-tuple.
  combined_data <- .combine(sample_names, seqnames, pos, counts, m, strand)
  
  # Construct GRanges
  gr <- GRanges(seqnames = combined_data[['seqnames']], ranges = IRanges(start = combined_data[['pos']][[1]], end = combined_data[['pos']][[m]]), strand = combined_data[['strand']], seqinfo = seqinfo)
  # Need to store other positions if m > 2
  if (m > 2){
    extra_pos <- lapply(X = seq(from = 2, to = m - 1, by = 1), FUN = function(i, pos){
      pos[[i]]
      }, pos = combined_data[['pos']])
    names(extra_pos) <- paste0('pos', seq(from = 2, to = m - 1, by = 1))
    } else{
      extra_pos = list()
    }
  assays <- SimpleList(c(combined_data[['counts']]))
  colData <- DataFrame(m = rep(m, length(combined_data[['sample_names']])), methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(combined_data[['sample_names']])), row.names = combined_data[['sample_names']])
  cometh <- SummarizedExperiment(assays = assays, rowData = gr, colData = colData)
  cometh <- .CoMeth(cometh, extra_pos = extra_pos)

  # Sort CoMeth object if required
  if (sort_cometh){
    cometh <- sort(cometh)
  }
  
  # Return CoMeth object
  return(cometh)
}

# CoMeth validity function
validCoMeth <- function(object){
  ## Check that m is identical for all samples
  ## Can't use getM(object) because getM checks whether object is a valid CoMeth object which it isn't when it gets called by the CoMeth constructor
  m <- colData(object)$m[1] 
  if (!zero_range(colData(object)$m)){
    msg <- validMsg(msg, "'m' must be identical for all samples")
  }
  ## Check assay names
  msg <- validMsg(NULL, .checkAssayNames(object, make_m_tuple_names(m)))
  ## Check rowData is GRanges
  if (class(rowData(object)) != "GRanges"){
    msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
  }
  ## Check that all 'counts' are non-negative
  # Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if (min(sapply(assays(object), min, na.rm = TRUE), na.rm = TRUE) < 0) {
    msg <- validMsg(msg, "'counts' has negative entries")
  }
  ## Check that all 'pos' are non-negative
  if (min(sapply(getPos(object)[, -1], min, na.rm = TRUE)) < 0) { 
    msg <- validMsg(msg, "'pos' has negative entries")
  }
  ## Check for existance and validity of 'extra_pos' slot
  if (!.hasSlot(object, "extra_pos")){
    msg <- validMsg(msg, "'extra_pos' slot is missing")
  } else{
    if ((m > 2) && (!is.list(object@extra_pos) || length(object@extra_pos) != (m - 2)) || !all(sapply(object@extra_pos, length) == nrow(object))){
      msg <- validMsg(msg, "'extra_pos' must be a list of length 'm' - 2 if 'm' > 2. Each element of 'extra_pos' must be of the same length as 'nrow(object)'")
    } else if (m == 2 && length(object@extra_pos) != 0){
        msg <- validMsg(msg, "'extra_pos' must be an empty list if 'm' = 2.")
      }
    } 
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}
setValidity("CoMeth", validCoMeth)

#### Methods for CoMeth objects ####

## TODO: Define all (necessary) BiocGenerics: 

setMethod(show, "CoMeth", function(object){
  cat("An object of type 'CoMeth' with:\n")
  cat(paste0(" ", nrow(object), " ", getMethylationType(object), " ", getM(object), "-tuples\n"))
  cat(paste0(" ", ncol(object), " samples\n"))
  #callNextMethod()
})

setMethod(length, "CoMeth", function(x){
  nrow(x)
})

# TODO: This doesn't work if j is also specified, e.g. CoMeth[1:10, 1:3]@extra_pos[[1]] is not of length 10
# The "[" method for a Cometh object is almost identical to that for a SummarizedExperiment object. The only difference is the for a CoMeth object we have to also take care to also subset the 'extra_pos' field.
# i indexes m-tuples, j indexes samples
setMethod("[", c("CoMeth", "ANY", "missing"),
          function(x, i, j, ..., drop = FALSE)
            {
            if (missing(j)){
              j <- seq_len(ncol(x))
            }            
            initialize(x, 
                       as(x, "SummarizedExperiment")[i, j, ..., drop = drop],
                       # x[i, j, ..., drop = drop], # Seems to be a circular definition
                       extra_pos = lapply(x@extra_pos, FUN = function(xx, i){xx[i]}, i = i))
          })

# The cbind method for a CoMeth object differs to that for a SummarizedExperiment. cbind allows for the objects to have different 'pos' data (getPos()) provided all objects have the same seqinfo. Samples without 'counts' data for a particular m-tuple will have the corresponding value filled with NA. The 'm' of all objects must be identical. 'methylation_type' may differ between object but the CoMeth object will have the combined 'methylation_type', e.g. cbind-ing a "CG" with "CHG" will result in a new object with a "CG/CHG" methylation_type.  sample_names must be unique across all objects being combined. The way cbind works is effectively to extract all the relevant information from the objects being combined and then create a new CoMeth object with the CoMeth constructor.
setMethod("cbind", function(..., deparse.level = 1){
  args <- unname(list(...))
  
  # Ensure that all objects have the same seqinfo. Don't want to rbind/cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed/cbind-ed have the same seqinfo. While it might be nice to to rbind/cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from mm9). Rather than handling each case in rbind/cbind I will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed/cbind-ed
  try(seq_info <- do.call("merge", lapply(args, seqinfo)))
  if (is(seq_info, "try-error")){
    stop("Can only rbind 'CoMeth' objects with identical 'seqinfo'. See ?seqinfo for help with replacing the 'seqinfo' of a 'CoMeth' object.")
  }
  
  # Ensure all sample_names are unique across the objects being combined
  sample_names <- unlist(lapply(args, sampleNames))
  if (length(sample_names) != length(unique(sample_names))){
    stop("'sample_names' across all CoMeth object must be unique.")
  }
  
  # Ensure all objects have the same 'm'
  m <- sapply(args, getM)
  if (!zero_range(m)){
    stop("Can only rbind 'CoMeth' objects with the same 'm'.")
  }
  
  # If necessary, combine 'methylation_type'
  methylation_type <- sapply(args, getMethylationType)
  if (!all(methylation_type == methylation_type[1])){
    warning("Combining 'CoMeth' objects with different 'methylation_type'. The 'methylation_type' of the new CoMeth object is: ", paste0(sort(unique(methylation_type)), collapse = '/'))
  }
  
  cometh <- CoMeth(sample_names = sample_names, seqnames = lapply())
})


setMethod("rbind", function(..., deparse.level = 1){
  args <- unname(list(...))
  
  # Ensure that all CoMeth objects have the same seqinfo. Don't want to rbind/cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed/cbind-ed have the same seqinfo. While it might be nice to to rbind/cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from mm9). Rather than handling each case in rbind/cbind I will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed/cbind-ed
  try(seq_info <- do.call("merge", lapply(args, seqinfo)))
  if (is(seq_info, "try-error")){
    stop("Can only rbind 'CoMeth' objects with identical 'seqinfo'. See ?seqinfo for help with replacing the 'seqinfo' of a 'CoMeth' object.")
  }
  
  
  

})

# OLD VERSION: 
# WARNING: Only order on rowData. Therefore there is no guarantee how it will order when m > 2
# setMethod(order, "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
#   args <- lapply(list(...), rowData)
#   do.call("order", c(args, list(na.last = na.last, decreasing = decreasing)))
# })

# Orders a CoMeth object based on rowData and extra_pos.
# TODO: Remove dependence on plyr::quickdf
# WARNING: I don't understand when you would use the case where "..." contains multiple objects but I've included it anyway based on selectMethod("order", "GenomicRanges")
setMethod(order, "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
# Based on selectMethod("order", "GenomicRanges")
  if (!isTRUEorFALSE(decreasing)){
    stop("'decreasing' must be TRUE or FALSE")
  }
  
  args <- list(...)
  
  if (!zero_range(sapply(args, getM))){
    stop("All 'CoMeth' objects must have the same 'm' value")
  }
  
  m <- getM(args[[1]])
  
  if (m < 3){
    ## Just defer to the order method defined for GRanges
    do.call("order", sapply(args, rowData))
  } else{
    ## Need to use an order method that takes note of the 'extra_pos' field in a CoMeth object
    if (length(args) == 1L) {
      x <- args[[1L]]
      do.call("order", c(cbind(as.factor(seqnames(x)), as.factor(strand(x)), start(x), plyr::quickdf(x@extra_pos), end(x)), list(na.last = na.last, decreasing = decreasing)))
      } else{
      order_args <- vector("list", (m + 2L) * length(args))
      idx <- (m + 2L) * seq_len(length(args))
      order_args[seq.int(from = 1, to = max(idx), by = m + 2)] <- lapply(args, function(x){
        as.factor(seqnames(x))
        })
      order_args[seq.int(from = 2, to = max(idx), by = m + 2)] <- lapply(args, function(x){
        as.factor(strand(x))
        })
      order_args[seq.int(from = 3, to = max(idx), by = m + 2)] <- lapply(args, start)
      order_args[seq.int(from = 4, to = max(idx), by = m + 2) + rep(seq(0, m - 3, by = 1))] <- lapply(args, function(x){unlist(lapply(x@extra_pos, function(xx){lapply(xx, c)}))})
      order_args[seq.int(from = 5, to = max(idx), by = m + 2)] <- lapply(args, function(x){end})
      order_args[idx] <- lapply(args, function(x){end(x)})
      do.call(order, c(order_args, list(na.last = na.last, decreasing = decreasing)))
    }
  }
})
  
setMethod(sort, "CoMeth", function(x, decreasing = FALSE, ...){
  x[order(x, decreasing = decreasing), ] 
})

setMethod("sampleNames", "CoMeth", function(object) {
  colnames(object)
})

setReplaceMethod("sampleNames", signature = signature(object = "CoMeth", value = "ANY"), function(object, value) {
  if (length(value) != length(sampleNames(object))){
    stop("Invalid 'sampleNames' length")
  }
  colnames(object) <- value
  object
  })

#### Functions that work on CoMeth objects ####

#' Obtain the \code{pos} from a \code{\link{CoMeth}} object
#' 
#' The \code{pos} of a \code{\link{CoMeth}} object are the genomic co-ordinates of the cytosines that make up each m-tuple.
#' @param x A CoMeth object
#'
#' @return A DataFrame containing the seqnames and positions (columns) of each m-tuple (rows)
#' @export
getPos <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  m <- getM(x)
  if (m == 1){
    pos  <- DataFrame(seqnames(x), start(x))
  } else if (m == 2) {
    pos <- DataFrame(seqnames(x), start(x), end(x))
  } else {
    pos <- DataFrame(seqnames(x), start(x), sapply(x@extra_pos, function(xx){xx}), end(x))
  }
  colnames(pos) <- c('seqnames', paste0('pos', seq_len(m))) 
  return(pos)
}

#' Get \code{m} from a \code{\link{CoMeth}} object
#'
#' Get the size of m-tuples, i.e. the \code{m} in m-tuples, from a \code{\link{CoMeth}} object.
#' @param x A CoMeth object
#'
#' @return An integer
#' @export
getM <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  m <- colData(x)[, 'm']
  if (!zero_range(m)){
    stop("'m' not equal for all samples in 'CoMeth' object")
  } else{
    return(m[1])
  }
}

#' Get \code{methylation_type} from a \code{\link{CoMeth}} object
#'
#' Get the size of methylation type, e.g. 'CG', 'CHH', 'CG/CHG', etc., from a \code{\link{CoMeth}} object.
#' @param CoMeth A CoMeth object
#'
#' @return An string describing the \code{methylation_type}
#' @export
getMethylationType <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  mt <- colData(x)[['methylation_type']][1]
  if (!all(mt == mt[1])){
    stop("'methylation_type' not equal for all samples in 'CoMeth' object")
  } else{
    return(mt[1])
  }
}

#' Get the \code{coverage} from a \code{\link{CoMeth}} object.
#'
#' The \code{coverage} of a \code{\link{CoMeth}} object is the number of reads for each m-tuple, i.e. \code{rowSums(getCounts(CoMeth))}. The current implementation is slow.
#' @param CoMeth A \code{\link{CoMeth}} object
#'
#' @return A numeric matrix. Each column of the matrix corresponds to a sample and each row to an m-tuple.
#' @export
getCoverage <- function(x) {
  stopifnot(is(x, "CoMeth"))
  Reduce("+", assays(x))
}
