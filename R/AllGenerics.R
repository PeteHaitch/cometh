### -------------------------------------------------------------------------
### getPos
###

## TODO: Figure out how to document these generics in the MTuples and/or CoMeth 
## help pages, as appropriate.
#' Extract the genomic positions from an \code{MTuples} object.
#' 
#' @param x An \code{\link{MTuples}} or \code{\link{CoMeth}} object.
#' 
#' @export
#' 
#' @return A numeric matrix. Each row of the matrix is an m-tuple and each 
#' column is a position. Note that this function does not return the seqnames
#' of the m-tuples; use the \code{\link{seqnames}} getter to do this.
setGeneric("getPos",  function(x, ...) {
  standardGeneric("getPos")
})

### -------------------------------------------------------------------------
### IPD
###

#' @export
setGeneric("getIPD", function(x, ...) {
  standardGeneric("getIPD")
})

### -------------------------------------------------------------------------
### getCoverage
###

## I accidentally defined getCoverage as a function, and not a method, in 
## R/methods-CoMeth-class.R. This caused the obscure error:
## "cyclic name space dependency detected when loading ‘cometh’, already 
## loading ‘cometh’"
## This is the exact same problem as described 
## here: https://stat.ethz.ch/pipermail/r-devel/2011-March/060100.html
#' @export
setGeneric("getCoverage", function(x, ...) {
  standardGeneric("getCoverage")
})

### -------------------------------------------------------------------------
### getM
###

#' @export
setGeneric("getM", function(x, ...) {
  standardGeneric("getM")
})

### -------------------------------------------------------------------------
### getMethylationType
###

#' @export
setGeneric("getMethylationType",function(x, ...) { 
  standardGeneric("getMethylationType")
})