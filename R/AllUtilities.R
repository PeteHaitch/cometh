### =========================================================================
### Helper functions not exported 
### =========================================================================

## TODO: Rewrite documentation but don't export function.
#' Check whether all elements of a numeric vector are identical (within machine precision)
#' @param x a numeric vector.
# 
#' @return TRUE if all elements of the vector are identical (within machine precision). FALSE in all other cases, including if the vector contains any NAs
#' 
#' @export
#' @keywords internal
#' 
#' @note This function is based on Hadley and John's answer to http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
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

## TODO: Check documentation and re-write as necessary.
#' An internal function used when constructing/combining CoMeth objects
#' 
#' Combine list-wise data into matrix-like data whilst taking care of common and sample-specific m-tuples, i.e. filling in NA for 'counts' when a sample doesn't have any observations for that m-tuple. 
#' 
#' @param m An \code{integer} storing the size of the m-tuples, i.e. the \code{m} in m-tuple. Only a single value is accepted, and not a list, because \code{m} must be the same for all samples in a \code{CoMeth} object.
#' @param sample_names A \code{character} vector containing the names of the samples. Sample names must be unique.
#' @param pos A \code{\link[IRanges]{DataFrameList}} \code{pos} must be named and the names must match those given in \code{sample_names}. The columns of each DataFrame must be: \code{seqnames}, \code{pos1}, ..., \code{posm}, where, for example, \code{posm} is \code{pos3} if \code{m} = 3. \code{seqnames} stored the sequence/chromosome name of the m-tuples. Therefore, the number of columns of each DataFrame is \code{m} + 1 and the number of rows is equal to the number of m-tuples for that particular sample.
#' @param counts A \code{\link[IRanges]{DataFrameList}}. \code{counts} must be named and the names must match those given in \code{sample_names}. The entry in each DataFrame corresponds to the number of times that particular co-methylation pattern (columns) was observed for that particular m-tuple (rows). Therefore, each DataFrame must have the same number of rows as its corresponding DataFrame in \code{pos} and have \eqn{2 ^ \code{m}} columns. The column names of each DataFrame must match that given by \code{make_m_tuple_names(m)}. 
#' @param strand An (optional)\code{\link[IRanges]{RleList}} object containing the strand information of each m-tuple. \strong{WARNING}: If \code{strand} is missing, all m-tuples in the resulting \code{\link{CoMeth}} object will have their strand set to \code{*} to signify that the strand is unknown or irrelevant (such as when methylation measurements have been combined across strands). \strong{WARNING}: m-tuples will not be combined across samples if they are on different strands.
#' 
#' @keywords internal
#' 
#' @return A list containing the combined seqnames, pos, counts and strand
.combine <- function(sample_names, seqnames, pos, strand, counts){
    
  ## Combine the data
  if (length(sample_names) > 1L){
    combined_coordinates <- DataFrame(seqnames = seqnames[[1L]], pos[[1L]], strand = strand[[1L]])
    combined_counts <- lapply(counts[[sample_names[[1L]]]], as.matrix)
    
    for (i in seq(from = 2L, to = length(sample_names), by = 1L)){
      ## Find any identical m-tuples between the current "combined" data and the "new sample"
      pair_coordinates <- rbind(combined_coordinates, DataFrame(seqnames = seqnames[[sample_names[[i]]]], pos[[sample_names[[i]]]], strand = strand[[sample_names[[i]]]]))
      in_both_idx1 <- duplicated(pair_coordinates) # Indexes *[[sample_names[[i]]]] (need to offset by NROW(combined_*))
      in_both_idx2 <- duplicated(pair_coordinates, fromLast = TRUE)
      in_both_combined_idx <- which(in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_both_new_sample_idx <- which(in_both_idx1[seq(from = nrow(combined_coordinates) + 1L, to = length(in_both_idx1), by = 1L)])
      in_just_combined_idx <- which(!in_both_idx2[seq_len(nrow(combined_coordinates))])
      in_just_new_sample_idx <- which(!in_both_idx1[seq(from = nrow(combined_coordinates) + 1L, to = length(in_both_idx1), by = 1L)])
      
      ## Combine the "combined" data and the "new sample" data to make the "new combined" data
      ## Always combine things in this order: (1) In both, (2) In just combined_*, (3) In just *[[sample_names[[i]]]]
      combined_coordinates <- rbind(combined_coordinates[in_both_combined_idx, ], combined_coordinates[in_just_combined_idx, ], DataFrame(seqnames = seqnames[[sample_names[[i]]]], pos[[sample_names[[i]]]], strand = strand[[sample_names[[i]]]])[in_just_new_sample_idx, ])
      combined_counts <- mapply(FUN = function(combined_counts, this_sample_counts, in_both_combined_idx, in_both_new_sample_idx, in_just_combined_idx, in_just_new_sample_idx){
      val <- rbind(cbind(combined_counts[in_both_combined_idx, , drop = FALSE], this_sample_counts[in_both_new_sample_idx]),
                   cbind(combined_counts[in_just_combined_idx, , drop = FALSE], matrix(NA_integer_, nrow = length(in_just_combined_idx), ncol = 1L)), 
                   cbind(matrix(NA_integer_, nrow = length(in_just_new_sample_idx), ncol = ncol(combined_counts)), this_sample_counts[in_just_new_sample_idx, , drop = FALSE]))
      return(val)
      }, combined_counts = combined_counts, this_sample_counts = lapply(counts[[sample_names[[i]]]], as.matrix), MoreArgs = list(in_both_combined_idx = in_both_combined_idx, in_both_new_sample_idx = in_both_new_sample_idx, in_just_combined_idx = in_just_combined_idx, in_just_new_sample_idx = in_just_new_sample_idx), SIMPLIFY = FALSE)
    }
    
    seqnames <- combined_coordinates[, 1L]
    pos <- as.matrix(combined_coordinates[, -c(1L, ncol(combined_coordinates))])
    strand <- combined_coordinates[, ncol(combined_coordinates)]
    counts <- combined_counts # list elements are already matrices  
    # TODO: Add names to counts?
  } else if (length(sample_names) == 1L){
    seqnames <- seqnames[[sample_names[[1L]]]]
    pos <- as.matrix(pos[[sample_names[[1L]]]])
    strand <- strand[[sample_names[[1L]]]]
    counts <- lapply(counts[[sample_names[[1L]]]], function(x){
      dim(x) <- c(length(x), 1L)
      return(x)
      }) # counts will go in the assays slot of a SummarizedExperiment-type object. Therefore counts must be coerced to matrices - this hack is a fast way to turn a vector into a column matrix.
      # TODO: Add names to counts?
  } else{
    seqnames <- Rle()
    pos <- matrix()
    strand <- Rle()
    counts <- list()
  }
  
  ## Return list of objects
  return(list(seqnames = seqnames, pos = pos, strand = strand, counts = counts))
}

## TODO: Write documentation
## NOTE: This function should only be called if m >= 3.
#' @export
#' 
#' @keywords internal
.findIdentical.MTuples <- function(query, subject, select){
  ## Idea1: Create seqnames:strand:pos vectors for both query and subject and then find all matches  
  ## This is ~6x slower than Idea 2 for a query with 1300 elements and a subject with 100000 elements. 
#   q_pos <- getPos(query)
#   q_k <- paste(seqnames(query), strand(query), do.call("paste", c(split(q_pos, c(col(q_pos))), list(sep = ":"))), sep = ":")
#   s_pos <- getPos(subject)
#   s_k <- paste(seqnames(subject), strand(subject), do.call("paste", c(split(s_pos, c(col(s_pos))), list(sep = ":"))), sep = ":")
#   ## Now, Find all matches and return a Hits object. This seems harder than it should be. It's easy to find the *first* match (using match(q_k, s_k)) but hard to find *all* matches.
#   tmp <- lapply(q_k, function(q, s_k){
#     which(s_k %in% q)
#   }, s_k = s_k)
  
  ## Idea 2: Create (m - 1) x 2-tuples (pairs) from the m-tuples as GRanges and run findOverlaps parameters chosen to select only exact matches.
  ## Then, intersect the Hits objects from each pair to create the final Hits object.
  m <- getM(query) # Assumes m is identical for query and subject
  
  ## Create an index variable
  ii <- seq_len(m - 1L)
  
  ## Create all pairs and run findOverlaps
  pair_hits <- lapply(ii, function(i, q_seqnames, s_seqnames, q_strand, s_strand, q_pos, s_pos){
    findOverlaps(query = GRanges(seqnames = q_seqnames, ranges = IRanges(start = q_pos[, i], end = q_pos[, i + 1]), strand = q_strand), subject = GRanges(seqnames = s_seqnames, ranges = IRanges(start = s_pos[, i], end = s_pos[, i + 1]), strand = s_strand), maxgap = 0L, minoverlap = 1L, type = "equal", select = select)
  }, q_seqnames = seqnames(query), s_seqnames = seqnames(subject), q_strand = strand(query), s_strand = strand(subject), q_pos = getPos(query), s_pos = getPos(subject))
  
  ## Intersect the Hits objects
  tuples_hits <- Reduce(intersect, pair_hits)
  
  return(tuples_hits)
}

#' Make m-tuple names
#'
#' This helper function constructs m-tuple names in the correct (alphabetical) order for a given value of m
#' @param m The size of the m-tuple. Must be an int.
#' 
#' @export
#' 
#' @keywords internal
#' 
#' @examples
#' .make_m_tuple_names(1L)
#' .make_m_tuple_names(2L)
#' .make_m_tuple_names(3L)
#' @return A character vector
.make_m_tuple_names <- function(m){
  if (!is.integer(m) || m < 1){
    stop("'m' must be an int. 'm' must be greater than 0.")
  }
  sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')}))))
}