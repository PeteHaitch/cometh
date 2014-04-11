### =========================================================================
### All classes 
### =========================================================================


### -------------------------------------------------------------------------
### MTuples 
###

.valid.MTuples.pos <- function(object){
  
  # The "empty" MTuples object
  if (isTRUE(all(is.na(start(object))))){
    m <- NA_integer_
  } else if (isTRUE(all(is.na(object@extraPos)))){
    # m = 1 or 2
    if (isTRUE(all(start(object) == end(object)))){
      m <- 1L
      pos <- as.matrix(start(object))
    } else{
      m <- 2L
      pos <- cbind(start(object), end(object))
    }
  } else{
    # m > 2
    m <- ncol(object@extraPos) + 2L
    pos <- cbind(start(object), object@extraPos, end(object))
  }
  
  if (!is.na(m)){
    if (!allRowsSorted(pos)){
      paste0("positions in each tuple must be sorted.")
    }
  }
  
  if (!is.na(m)){
    if (min(pos) < 0){
      # min(x) < 0 is faster than any(x < 0)
      paste0("positions in each tuple must be positive integers.")
    }
  }
}

.valid.MTuples <- function(object){
  c(.valid.MTuples.pos(object)) # Include all .valid.MTuples.* functions in this vector
}

##
setClass('MTuples', 
         representation(extraPos = "matrix"), # What to do when m < 3? Currently store as matrix of NAs with 1 column. Perhaps better to store as NULL/NA?
         contains = "GRanges",
         validity = .valid.MTuples
         )

### -------------------------------------------------------------------------
### CoMeth 
###

.valid.CoMeth <- function(object){
  c() # Include all .valid.CoMeth.* functions in this vector
}

CoMeth <- setClass('CoMeth', 
                    contains="SummarizedExperiment",
                    validity = .valid.CoMeth)