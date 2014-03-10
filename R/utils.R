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
#' @note Should be an internal function. Based on Hadley and John's answer to http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
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