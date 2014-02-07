# QUESTION: Should this be a function of a method?
# TODO: Don't export, this is an internal function

#' An internal function
#' @param object A CoMeth object
#' 
#' @return A character vector wit the names of the assays
#' @export
#' @note Should be an internal function. Copied from bsseq package
assayNames <- function(object) {
  names(object@assays$field("data"))
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
  if(!all(names %in% nms))
    return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                   class(object), paste0(names, collapse = ", ")))
  else
    NULL
}