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
#' @param m An \code{integer} storing the size of the m-tuples, i.e. the \code{m} in m-tuple. Only a single value is accepted, and not a list, because \code{m} must be the same for all samples in a \code{CoMeth} object.
#' @param sample_names A \code{character} vector containing the names of the samples. Sample names must be unique.
#' @param pos A \code{\link[IRanges]{DataFrameList}} \code{pos} must be named and the names must match those given in \code{sample_names}. The columns of each DataFrame must be: \code{seqnames}, \code{pos1}, ..., \code{posm}, where, for example, \code{posm} is \code{pos3} if \code{m} = 3. \code{seqnames} stored the sequence/chromosome name of the m-tuples. Therefore, the number of columns of each DataFrame is \code{m} + 1 and the number of rows is equal to the number of m-tuples for that particular sample.
#' @param counts A \code{\link[IRanges]{DataFrameList}}. \code{counts} must be named and the names must match those given in \code{sample_names}. The entry in each DataFrame corresponds to the number of times that particular co-methylation pattern (columns) was observed for that particular m-tuple (rows). Therefore, each DataFrame must have the same number of rows as its corresponding DataFrame in \code{pos} and have \eqn{2 ^ \code{m}} columns. The column names of each DataFrame must match that given by \code{make_m_tuple_names(m)}. 
#' @param strand An (optional)\code{\link[IRanges]{RleList}} object containing the strand information of each m-tuple. \strong{WARNING}: If \code{strand} is missing, all m-tuples in the resulting \code{\link{CoMeth}} object will have their strand set to \code{*} to signify that the strand is unknown or irrelevant (such as when methylation measurements have been combined across strands). \strong{WARNING}: m-tuples will not be combined across samples if they are on different strands.
#' 
#' @return A list containing the combined seqnames, pos, counts and strand
#' @note Should be an internal function.
.combine <- function(m, sample_names, pos, counts, strand){
  if (length(sample_names) > 1){
    combined_coordinates <- DataFrame(pos[[1]], strand = strand[[1]])
    combined_counts <- lapply(counts[[sample_names[[1]]]], as.matrix)
    
    for (i in seq(from = 2, to = length(sample_names), by = 1)){
      ## Find any identical m-tuples between the current "combined" data and the "new sample"
      pair_coordinates <- rbind(combined_coordinates, DataFrame(pos[[sample_names[i]]], strand = strand[[sample_names[i]]]))
      in_both_idx1 <- duplicated(pair_coordinates) # Indexes *[[sample_names[i]]] (need to offset by NROW(combined_*))
      in_both_idx2 <- duplicated(pair_coordinates, fromLast = TRUE)
      in_both_combined_idx <- which(in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_both_new_sample_idx <- which(in_both_idx1[seq(from = nrow(combined_coordinates) + 1, to = length(in_both_idx1), by = 1)])
      in_just_combined_idx <- which(!in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_just_new_sample_idx <- which(!in_both_idx1[seq(from = nrow(combined_coordinates) + 1, to = length(in_both_idx1), by = 1)])
      
      ## Combine the "combined" data and the "new sample" data to make the new "combined" data
      ## Always combine things in this order: (1) In both, (2) In just combined_*, (3) In just *[[sample_names[i]]] 
      combined_coordinates <- rbind(combined_coordinates[in_both_combined_idx, ], combined_coordinates[in_just_combined_idx, ], DataFrame(pos[[sample_names[i]]], strand = strand[[sample_names[i]]])[in_just_new_sample_idx, ])
      combined_counts <- mapply(FUN = function(combined_counts, this_sample_counts, in_both_combined_idx, in_both_new_sample_idx, in_just_combined_idx, in_just_new_sample_idx){
      val <- rbind(cbind(combined_counts[in_both_combined_idx, ], this_sample_counts[in_both_new_sample_idx]),
                   cbind(combined_counts[in_just_combined_idx, , drop = FALSE], matrix(NA_integer_, nrow = length(in_just_combined_idx), ncol = 1)), 
                   cbind(matrix(NA_integer_, nrow = length(in_just_new_sample_idx), ncol = ncol(combined_counts)), this_sample_counts[in_just_new_sample_idx, , drop = FALSE]))
      return(val)
      }, combined_counts = combined_counts, this_sample_counts = lapply(counts[[sample_names[[i]]]], as.matrix), MoreArgs = list(in_both_combined_idx = in_both_combined_idx, in_both_new_sample_idx = in_both_new_sample_idx, in_just_combined_idx = in_just_combined_idx, in_just_new_sample_idx = in_just_new_sample_idx), SIMPLIFY = FALSE)
    }
    
    coordinates <- combined_coordinates
    counts <- combined_counts
    
  } else{
    coordinates <- DataFrame(pos[[sample_names[1]]], strand[[sample_names[1]]])
    counts <- lapply(counts[[sample_names[[1]]]], as.matrix) # assays slots must be matrix-like objects
  }
  return(list(coordinates = coordinates, counts = counts))
}
  