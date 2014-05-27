### =========================================================================
### MTuples: tuples of genomic positions
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

#' @include AllGenerics.R
#' @export
setMethod("getPos", "MTuples", function(x){
  
  m <- getM(x)
  
  if (is.na(m)){
    pos <- matrix()
  } else {
    if (m == 1L){
      pos <- as.matrix(start(x))    
    } else if (m == 2L){
      pos <- cbind(start(x), end(x))
      } else{
        # m > 2
        pos <- cbind(start(x), x@extraPos, end(x))
      }
    colnames(pos) <- paste0('pos', seq_len(m))
  }
  
  return(pos)
})

#' @include AllGenerics.R
#' @export
setMethod("getM", "MTuples", function(x){
  
  # The "empty" MTuples object
  if (isTRUE(all(is.na(start(x))))){
    m <- NA_integer_
  } else if (isTRUE(all(is.na(x@extraPos)))){
    # m = 1 or 2
    if (isTRUE(all(start(x) == end(x)))){
      m <- 1L
    } else{
      m <- 2L
    }
  } else{
    # m > 2
    m <- ncol(x@extraPos) + 2L
  }
  
  return(m)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
MTuples <- function(seqnames = Rle(), 
                    pos = matrix(), 
                    strand = Rle("*", length(seqnames)), 
                    ...,
                    seqlengths = NULL, 
                    seqinfo = NULL){
  ## Check that all required arguments are not missing
  if (missing(seqnames)){
    stop(sQuote("seqnames"), " missing.\nPlease see the help page for MTuples, which can accessed by typing ", 
         sQuote("?MTuples"), 
         " at the R prompt, for further details of this argument.")
  }
  if (missing(pos)){
    stop(sQuote("pos"), " missing.\nPlease see the help page for MTuples, which can accessed by typing ", 
         sQuote("?MTuples"), 
         " at the R prompt, for further details of this argument.")
  }
  
  ## Check that all required arguments are of the correct type
  ## Or try to coerce to the correct type
  if (!is.matrix(pos)){
    stop(sQuote("pos"), " must be a matrix.\nPlease see the help page for MTuples, which can accessed by typing ", sQuote("?MTuples"), " at the R prompt, for further details of this argument.") 
  }
  
  ## Other argument checks are deferred to the GRanges constructor, which
  ## is called further down.
  m <- ncol(pos)

  if (!isTRUE(all(is.na(pos)))){
    if (m == 1L){
      ranges <- IRanges(start = pos[, 1L], width = 1L)
      extraPos <- matrix(NA_integer_, nrow = length(ranges))
    } else if (m == 2L){
      ## When m = 2, need extra check that start != pos, because this isn't caught by the validity methods when a single m-tuple is passed.
      ## Otherwise this supposed 2-tuple is silently converted to a 1-tuple.
      ## See GitHub issue #8 (https://github.com/PeteHaitch/cometh/issues/8)
      if (!.allRowsSortedCpp(pos)){
        stop(paste0("positions in each m-tuple must be sorted in strictly increasing order, i.e. ", sQuote('pos1'), " < ", sQuote('pos2'), " < ", sQuote('...'), " < ", sQuote('posm')))
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
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, seqlengths = seqlengths, seqinfo = seqinfo, ...)
  
  new("MTuples", gr, extraPos = extraPos)
  
} 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

## Use 'c', with the optional ignore.mcols argument, as inherited from GRanges. Don't define rbind, cbind or combine as these aren't defined for GRanges. 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# TODO: None at this point. Coercion to GRanges, DataFrame and matrix might be useful.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

#' @include AllGenerics.R
#' @export
setMethod("getIPD", "MTuples", function(x){
  m <- getM(x)
  if (m == 1L){
    stop("It does not make sense to compute IPD when ", sQuote('m'), " = 1.")
  } else{
    .rowDiffsCpp(getPos(x)) # matrixStats::rowDiffs(getPos(x)) is a (slower) alternative
  }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compare() and related methods.
###
### compare() is based on the method defined for `GRanges` objects but cannot 
### explicity inherit the method due to the `extraPos` slot in `MTuples` 
### objects.
### However, unlike the compare() method defined for `GRanges`, the compare() 
### method for `MTuples` requires that both `x` and `y` have the same length.
### I also define the element wise (aka "parallel") operators '<=' and '=='.
### The other element wise operators (`!=`, `>=`, `<`, `>`) work out-of-the-box 
### on `MTuples` objects via inheritance from `GRanges` -> `Vector`.

# This is adapted from GenomicRanges:::.GenomicRanges.compare.
#' @export
.MTuples.compare <- function(x, y){
  
  if (length(x) != length(y)){
    stop("Cannot ", sQuote('compare'), " ", sQuote('MTuples'), " objects when ",
         sQuote('length(x)'), " != ", sQuote('length(y)'))
  }
  
  ## Check 'm' is identical
  if (getM(x) != getM(y)){
    stop("Cannot ", sQuote('compare'), " ", sQuote('MTuples'), 
         " objects with different ", sQuote('m'), ".")
  }
  
  ## Pre-comparison step (see above for details).
  ## merge() will fail if 'x' and 'y' don't have compatible underlying
  ## sequences.
  seqinfo <- merge(seqinfo(x), seqinfo(y))
  seqlevels <- seqlevels(seqinfo)
  if (any(diff(match(seqlevels(y), seqlevels)) < 0L)){
    stop("the 2 objects to compare have ", sQuote('seqlevels'), 
         " in incompatible orders.") # Error message differs slightly from that provided by .GenomicRanges.compare
  }
  ## This should only insert new seqlevels in the existing ones i.e. it
  ## should NEVER drop or reorder existing levels
  seqlevels(x) <- seqlevels(y) <- seqlevels
  
  ## This is where .MTuples.compare really differs from .GenomicRanges.compare
  a <- as.integer(seqnames(x)) - as.integer(seqnames(y))
  b <- as.integer(strand(x)) - as.integer(strand(y))
  c <- getPos(x) - getPos(y)
  
  ## Loop over cbind(a, b, c) by row and report the first non-zero element or 
  ## the final element.
  ## Do this without actually forming cbind(a, b, c).
  ## Rcpp solution is > 1000x faster than equivalent R solution when 
  ## length(x) = 2,000,000
  val <- .compareMTuplesCpp(a, b, c)
  
  return(val)
}

setMethod("compare", c("MTuples", "MTuples"), function(x, y) {
  .MTuples.compare(x, y)
})

setMethod("<=", c("MTuples", "MTuples"), function(e1, e2) {
  .MTuples.compare(e1, e2) <= 0L
})

setMethod("==", c("MTuples", "MTuples"), function(e1, e2) {
  .MTuples.compare(e1, e2) == 0L
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### duplicated()  
###
### Can't rely on using duplicated() via inheritance from `Vector` because 
### `MTuples` inherits from `GRanges`, rather than from `Vector` directly.
### Furthermore, simply deferring to the duplicated() method for GRanges would 
### only check pos1, posm of 'pos', which is okay provided m < 3 but not if 
### m >= 3.
### While I could defer to the GRanges method when m < 3, I now have a fast 
### method that works for all m.
###
### anyDuplicated() will be very slightly faster than any(duplicated()).
###
### unique() will work out-of-the-box on an `MTuples` object thanks to the 
### method for `GRanges` objects, which inherits from `Vector` objects.

## TODO: Perhaps move rev to mTuplesHashCpp where it will be faster.
.duplicated.MTuples <- function(x, incomparables = FALSE, fromLast = FALSE){
  if (!identical(incomparables, FALSE)){
    stop(sQuote('duplicated'), " method for ", sQuote('MTuples'), 
         " objects only accepts ", sQuote('incomparables = FALSE'))
  }
  ## Uses the "rowSums hash" approach.
  ## It is very fast and produces identical results to 
  ## base::duplicated.array(x, MARGIN = 1).
  ## .candidateDuplicateMTuplesCpp only returns __candidate__ duplicates.
  ## These candidates need to be further checked using the (slower) 
  ## base::duplicated.array method
  
  if (!fromLast){
    a <- as.integer(seqnames(x))
    b <- as.integer(strand(x))
    C <- getPos(x)
  } else{
    a <- rev(as.integer(seqnames(x)))
    b <- rev(as.integer(strand(x)))
    C <- getPos(x)[seq.int(from = length(x), to = 1), , drop = FALSE]
  }
    d <- .candidateDuplicateMTuplesCpp(a = as.integer(seqnames(x)), 
                        b = as.integer(strand(x)), 
                        C = getPos(x))
    d[d] <- duplicated(cbind(a[d], b[d], C[d, , drop = FALSE], 
                             deparse.level = 0), 
                       incomparables = incomparables, 
                       MARGIN = 1, 
                       fromLast = fromLast)
  return(d)  
}

# duplicated() checks seqnames, strand and positions of each m-tuple. 
## TODO: Is there is a need for an S3/S4 combo? If so, why?
### S3/S4 combo for duplicated.MTuples
#' @export
duplicated.MTuples <- function(x, incomparables = FALSE, ...){
  .duplicated.MTuples(x, incomparables = incomparables, ...)
}
#' @export
setMethod("duplicated", "MTuples", .duplicated.MTuples)

## TODO: .anyDuplicated() should return 0 if no duplicates and index 
## of first duplicate __NOT__ TRUE/FALSE
.anyDuplicated.MTuples <- function(x, incomparables = FALSE){
  if (!identical(incomparables, FALSE)){
    stop(sQuote('anyDuplicated'), " method for ", sQuote('MTuples'), 
         " objects only accepts ", sQuote('incomparables = FALSE'))
  }
  
  a <- as.integer(seqnames(x))
  b <- as.integer(strand(x))
  C <- getPos(x)
  d <- .candidateDuplicateMTuplesCpp(a = as.integer(seqnames(x)), 
                                    b = as.integer(strand(x)), 
                                    C = getPos(x))
  if (isTRUE(any(d))){
    # add (if non-zero) is the index of the first duplicated element in d[d], 
    # instead of d (which is what we really want).
    add <- anyDuplicated(cbind(a[d], b[d], C[d, , drop = FALSE], deparse.level = 0))
    if (!identical(add, 0L)){
      ad <- which(d)[add]
    } else{
      ad <- 0L
    }
  } else{
    ad <- 0L
  }
  return(ad)
}

# anyDuplicated() checks seqnames, strand and positions of each m-tuple.
# Will be very slightly faster than any(duplicated()).
## TODO: Is there is a need for an S3/S4 combo? If so, why?
### S3/S4 combo for duplicated.MTuples
#' @export
anyDuplicated.MTuples <- function(x, incomparables = FALSE, ...){
  .anyDuplicated.MTuples(x, incomparables = incomparables, ...)
}
#' @export
setMethod("anyDuplicated", "MTuples", .anyDuplicated.MTuples)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match()
###
### %in%, findMatches(), countMatches() will work out-of-the-box on `MTuples` 
### objects thanks to the method for Vector objects.

## Loosely based on 'match' method for GenomicRanges objects
#' @export
setMethod("match", c("MTuples", "MTuples"), 
          function(x, table, nomatch = NA_integer_, incomparables = NULL, 
                   ignore.strand = FALSE){
  
  if (!isSingleNumberOrNA(nomatch)){
    stop(sQuote('nomatch'), " must be a single number or ", sQuote('NA'))
  }
  if (!is.integer(nomatch)){
    nomatch <- as.integer(nomatch)
  }
  if (!is.null(incomparables)){
    stop(sQuote('match'), " method for ", sQuote('MTuples'), 
         " objects only accepts ", sQuote('incomparables = NULL'))
  }
  if (!isTRUEorFALSE(ignore.strand)){
    stop(sQuote('ignore.strand'), " must be ", sQuote('TRUE'), " or ", 
         sQuote('FALSE'))
  }
  ## Calling merge() is the way to check that 'x' and 'table' are based on the 
  ## same reference genome.
  merge(seqinfo(x), seqinfo(table))
  
  findOverlaps(x, table, type = "equal", select = "first")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### order() and related methods.
###
### The order() and rank() methods for MTuples objects are consistent with the 
### order implied by compare().
### sort() is defined separately, rather than simply inherited from `GRanges`,
### because it does not use the `by` argument.


# Loosely based on 'order' method for GRanges
#' @export
setMethod(order, "MTuples", function(..., na.last = TRUE, decreasing = FALSE){
  
  if (!isTRUEorFALSE(decreasing)){
    stop("'decreasing' must be TRUE or FALSE")
  }
  
  args <- list(...)
  
  if (!.zero_range(sapply(args, getM))){
    stop("All ", sQuote('MTuples'), " objects must have the same ", sQuote('m'),
         " value")
  }
  
  m <- getM(args[[1]])
  
  if (m < 3){
    ## If m < 3 just defer to the order method defined for GRanges
    ## Can't simply use callNextMethod() because (oddly) there is no order 
    ## method defined for GRanges (rather it is defined for GenomicRanges).
    order(do.call("c", lapply(args, function(x){as(x, "GRanges")})), 
          na.last = na.last, decreasing = decreasing)
  } else{
    ## If m >= 3 then need to define an order method that takes note of the 
    ## 'extra positions' in each m-tuple
    order_args <- vector("list", (m + 2L) * length(args))
    idx <- (m + 2L) * seq_len(length(args))
    order_args[seq.int(from = 1, to = max(idx), by = m + 2)] <- lapply(args, function(x){
      as.factor(seqnames(x))
    })
    order_args[seq.int(from = 2, to = max(idx), by = m + 2)] <- lapply(args, function(x){
      as.factor(strand(x))
    })
    order_args[seq.int(from = 3, to = max(idx), by = m + 2)] <- lapply(args, start)
    order_args[seq.int(from = 4, to = max(idx), by = m + 2) + rep(seq(0, m - 3, by = 1))] <- lapply(args, function(x, m){getPos(x)[, -c(1, m)]}, m = m)
    order_args[seq.int(from = 5 + m - 3, to = max(idx), by = m + 2)] <- lapply(args, function(x){end})
    order_args[idx] <- lapply(args, function(x){end(x)})
    do.call(order, c(order_args, list(na.last = na.last, decreasing = decreasing)))
  }
})

# Loosely based on 'sort' method for GRanges
#' @export
setMethod(sort, "MTuples", function(x, decreasing = FALSE, ignore.strand = FALSE, by){
  
  if (!missing(by)){
    stop("Sorry, the ", sQuote('by'), " argument is not currently implemented.")
  }
  
  if (!isTRUEorFALSE(ignore.strand)){
    stop(sQuote('ignore.strand'), " must be ", sQuote('TRUE'), " or ", sQuote('FALSE'))
  }
  if (ignore.strand) {
    x2 <- unstrand(x)
    i <- order(x2, decreasing = decreasing)
  } else{
    i <- order(x, decreasing = decreasing)
  }
  x[i, , drop = FALSE]
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

## Ensure the extraPos column "sticks" during subsetting, etc.
setMethod(GenomicRanges:::extraColumnSlotNames, "MTuples",
          function(x) {
            c("extraPos")
})

## The show method is adapted from that of GRanges
.makeNakedMatFromMTuples <- function(x){
  lx <- length(x)
  nc <- ncol(mcols(x))
#   m <- getM(x)
#   if (isTRUE(all(is.na(x@extraPos)))){
#     # m = 1 or 2
#     if (m == 1L){
#       pos <- as.matrix(start(x))
#     } else if (m == 2L)
#       pos <- cbind(start(x), end(x))
#     }
#   } else {
#     pos <- cbind(start(x), x@extraPos, end(x))
#   }
#   colnames(pos) <- paste0('pos', seq_len(m))
  pos <- getPos(x)
  ans <- cbind(seqnames = as.character(seqnames(x)), pos, strand = as.character(strand(x)))
  
  ## This code is commented out because otherwise it will repeat the extraPos field in the output
  #extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
  #  if (length(extraColumnNames) > 0L) {
  #      ans <- do.call(cbind, c(list(ans), lapply(GenomicRanges:::extraColumnSlots(x), 
  #          showAsCell)))
  #  }
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), IRanges:::showAsCell), 
                                 list(check.names = FALSE)))
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
  }
  return(ans)
}


## This line is commented out because the show method for MTuples cannot currently use the with.classinfo argument.
#showMTuples <- function(x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE){
showMTuples <- function(x, margin = "", print.seqlengths = FALSE){
  lx <- length(x)
  nc <- ncol(mcols(x))
  if (isTRUE(all(is.na(x@extraPos)))){
    # m = 1 or 2
    if (isTRUE(all(start(x) == end(x)))){
      m <- 1
    } else{
      m <- 2
    }
  } else{
    m <- ncol(x@extraPos) + 2
  }
  
  cat(class(x), " with ", lx, " x ", ifelse(lx == 1L, paste0(m, "-tuple"), 
                                            paste0(m, "-tuples")), " and ", nc, " metadata ", ifelse(nc == 1L, 
                                                                                                     "column", "columns"), ":\n", sep = "")
  
  
  out <- IRanges:::makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromMTuples)
  ## These lines commented out because classinfo is more complicated for MTuples objects than GRanges objects. For example, some of the `pos` information is stored in an IRanges object while some is stored in a matrix.
  #if (with.classinfo) {
  #    .COL2CLASS <- c(seqnames = "Rle", ranges = "IRanges", 
  #        strand = "Rle")
  #    extraColumnNames <- extraColumnSlotNames(x)
  #    .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
  #    classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
  #    stopifnot(identical(colnames(classinfo), colnames(out)))
  #    out <- rbind(classinfo, out)
  #}
  
  if (nrow(out) != 0L){ 
    rownames(out) <- paste0(margin, rownames(out))
  }
  print(out, quote = FALSE, right = TRUE)
  if (print.seqlengths) {
    cat(margin, "---\n", sep = "")
    GenomicRanges:::showSeqlengths(x, margin = margin)
  }
}

setMethod("show", "MTuples",function(object){
  showMTuples(object, margin="  ", print.seqlengths = TRUE)
})
