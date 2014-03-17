# QUESTION: Should this be a function of a method?
# TODO: Don't export, this is an internal function

#' An internal function
#' @param object A CoMeth object
#'
#' @return A character vector wit the names of the assays
#' @export
#' @note Should be an internal function. Copied from bsseq package
assayNames <- function(object) {
  names(assays(object))
}

# An internal function
#' @param object A CoMeth object
#' @param names The expected assay names
#'
#' @return NULL if successful, otherwise a message
#' @export
#' @note Should be an internal function. Copied from bsseq package
.checkAssayNames <- function(object, names) {
  nms <- assayNames(object)
  if(!all(names %in% nms)){
    return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                   class(object), paste0(names, collapse = ", ")))
  } else{
    NULL
  }
}

# An internal function
#' @param x a numeric vector.
#' 
#' @return TRUE if all elements of the vector are identical (within machine precision). FALSE in all other cases, including if the vector contains any NAs
#' @export
#' @note Should be an internal function. This function is based on Hadley and John's answer to http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) {
    val <- TRUE
  } 
  if (any(is.na(x)) & any(is.na(x))){
    val <- FALSE
  } else{
    val <- (abs(max(x) - min(x)) < tol)
  }
  
  return(val)
}

#' An internal function for constructing/combining CoMeth objects
#' Combine list-wise data into matrix-like data whilst taking care of common and sample-specific m-tuples, i.e. filling in zeros for 'counts' when a sample doesn't have any observations for that m-tuple. 
#' @param sample_names A character vector containing the names of the samples.
#' @param seqnames A list of character vectors, a list of factors or an \code{\link[IRanges]{RleList}} object containing the sequence names. Each element of the list should be named and match one of the \code{sample_names}.
#' @param pos A list of lists. Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the \code{sample_names}. The number of elements of each sub-list should be equal to \code{m}. For a given sample, each element of the inner-list stores the positions of each m-tuple as an integer vector. E.g. \code{pos[[2]][[1]]} contains for \code{sample2}, all \code{pos1} for each m-tuple as an integer vector.
#' @param counts Must be a list of (integer) lists. Each element of the outer-list should be named and match one of the \code{sample_names}. The number of elements of each sub-list should be equal to \eqn{2 ^ \code{m}}. For a given sample, each element of the inner-list stores the number of times each co-methylation pattern is observed at each m-tuple. E.g. \code{pos[[2]][[1]]} contains for \code{sample2}, how many times the first co-methylation pattern was observed for each m-tuple as an integer vector.
#' @param strand An list of character vectors, list of factors or an \code{\link[IRanges]{RleList}} object containing the strand information of each m-tuple. \strong{WARNING}: m-tuples will not be combined across samples if they are on different strands. 
#' @param m An integer storing the size of the m-tuples, i.e. the \code{m} in m-tuple. Only a single value is accepted, and not a list, because \code{m} must be the same for all samples in a \code{CoMeth} object.
#' 
#' @return A list containing the combined seqnames, pos, counts and strand
#' @note Should be an internal function.
.combine <- function(sample_names, seqnames, pos, counts, m, strand){
  if (length(sample_names) > 1){
    combined_seqnames <- seqnames[[sample_names[1]]]
    combined_strand <- strand[[sample_names[[1]]]]
    combined_coordinates <- matrix(unlist(pos[[sample_names[[1]]]], use.names = FALSE), ncol = m)
    combined_counts <- lapply(counts[[sample_names[1]]], as.matrix)
    for ( i in seq(from = 2, to = length(sample_names), by = 1)){
      pair_seqnames <- c(as.character(combined_seqnames), as.character(seqnames[[sample_names[i]]])) # as.character necessary for if seqnames are factors
      pair_strand <- c(combined_strand, strand[[sample_names[[i]]]])
      pair_coordinates <- cbind(as.factor(pair_seqnames), factor(as.vector(pair_strand), levels = c('+', '-', '*')), rbind(combined_coordinates, matrix(unlist(pos[[sample_names[[i]]]], use.names = FALSE), ncol = m))) # Have to add seqnames and strand as a factor
      in_both_idx1 <- duplicated(pair_coordinates, MARGIN = 1) # Indexes *[[sample_names[i]]] (need to offset by NROW(combined_*))
      in_both_idx2 <- duplicated(pair_coordinates, MARGIN = 1, fromLast = TRUE) # Indexes combined_* (for the first nrow(combined_*) elements)
      in_both_combined_idx <- which(in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_both_new_sample_idx <- which(in_both_idx1[seq(from = nrow(combined_coordinates) + 1, to = length(in_both_idx1), by = 1)])
      in_just_combined_idx <- which(!in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_just_new_sample_idx <- which(!in_both_idx1[seq(from = nrow(combined_coordinates) + 1, to = length(in_both_idx1), by = 1)])
      ## Always combine things in this order: (1) In both, (2) In combined_*, (3) In *[[sample_names[i]]]
      combined_seqnames <- c(as.character(combined_seqnames[in_both_combined_idx]), as.character(combined_seqnames[in_just_combined_idx]), as.character(seqnames[[sample_names[i]]][in_just_new_sample_idx]))
      combined_strand <- c(combined_strand[in_both_combined_idx], combined_strand[in_just_combined_idx], strand[[sample_names[i]]][in_just_new_sample_idx])
      combined_coordinates <- rbind(combined_coordinates[in_both_combined_idx, ], combined_coordinates[in_just_combined_idx, ], matrix(unlist(pos[[sample_names[[i]]]], use.names = FALSE), ncol = m)[in_just_new_sample_idx, ])
      combined_counts <- mapply(FUN = function(combined_counts, this_sample_counts, in_both_combined_idx, in_both_new_sample_idx, in_just_combined_idx, in_just_new_sample_idx){
        val <- rbind(cbind(combined_counts[in_both_combined_idx, ], this_sample_counts[in_both_new_sample_idx]), 
                     cbind(combined_counts[in_just_combined_idx, ], NA_integer_), 
                     cbind(matrix(NA_integer_, nrow = length(in_just_new_sample_idx), ncol = ncol(combined_counts)), this_sample_counts[in_just_new_sample_idx]))
        return(val)
      }, 
      combined_counts = combined_counts, this_sample_counts = counts[[sample_names[[i]]]], MoreArgs = list(in_both_combined_idx = in_both_combined_idx, in_both_new_sample_idx = in_both_new_sample_idx, in_just_combined_idx = in_just_combined_idx, in_just_new_sample_idx = in_just_new_sample_idx), SIMPLIFY = FALSE)
    }
    seqnames <- combined_seqnames
    strand <- combined_strand
    cc <- seq_len(ncol(combined_coordinates))
    names(cc) <- paste0('pos', cc)
    pos <- lapply(cc, FUN = function(cc, combined_coordinates){
      combined_coordinates[, cc]
    }, combined_coordinates = combined_coordinates)
    counts <- combined_counts
    counts <- lapply(counts, FUN = function(x, sample_names){
      colnames(x) <- sample_names
      return(x)
    }, sample_names = sample_names)
  } else{
    seqnames <- seqnames[[sample_names]]
    pos <- pos[[sample_names]]
    counts <- lapply(counts[[sample_names]], as.matrix) # assays slots must be matrix-like objects
    strand <- strand[[sample_names]]
  }
  return(list(seqnames = seqnames, pos = pos, counts = counts, strand = strand))
}
  