### =========================================================================
### MTuples: tuples of genomic positions
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

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

MTuples <- function(seqnames = Rle(), pos = matrix(), strand = Rle("*", length(seqnames)), ..., seqlengths = NULL, seqinfo = NULL){

  ## TODO: Argument checks,
  m <- ncol(pos)
  
  if (!isTRUE(all(is.na(pos)))){
    if (m == 1L){
      ranges <- IRanges(start = pos[, 1L], width = 1L)
      extraPos <- matrix(NA_integer_, nrow = length(ranges))
    } else if (m == 2L){
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
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, ..., seqlengths = seqlengths, seqinfo = seqinfo)
  
  new("MTuples", gr, extraPos = extraPos)
  
} 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

# Combine 2 MTuples objects.
# To combine more than 2 MTuples object use Reduce(comine, list), where list is a list of MTuples objects.
setMethod("combine", c("MTuples", "MTuples"), function(x, y, ...){
  # Idea 1: Should use the .combine function with the 'sample_names' argument set to something artificial and the 'counts' argument set to NULL
#   sample_names <- seq_len(2)
#   seqnames <- list(seqnames(x), seqnames(y))
#   names(seqnames) <- sample_names
#   pos <- list(getPos(x), getPos(y))
#   names(pos) <- sample_names
#   strand <- list(strand(x), strand(y))
#   names(strand) <- sample_names
#   counts <- NULL
#   ## TODO: What about seqinfo and mcols?
#   
#   combined_data <- .combine(sample_names = sample_names, seqnames = seqnames, pos = pos, strand = strand, counts = counts)
#   MTuples(seqnames = combined_data$seqnames, pos = combined_data$pos, combined_data$strand, ..., )
  
  
  ## Idea 2: Use z <- c(x, y) and then unique(z), for a suitably defined unique method
  ## Can actually just the unique method that MTuples objects inherit from Vector objects. 
  
  ## Check m is the same for both x and y
  if (!identical(getM(x), getM(y))){
    stop("Can only combine ", sQuote("MTuples"), " object with the same ", sQuote("m"))
  }
  
  ## FIXME: "z <- c(x, y)" will not work if !identical(colnames(mcols(x)), colnames(mcols(y))).
  z <- c(x, y)
  ## FIXME "unique(z)" will ignore mcols(x) and mcols(y) when combining in favour of just keeping mcols(x)
  return(unique(z))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

# None at this point.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Filtering
###

## TODO: Filter by IPD?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setMethod("IPD", "MTuples", function(x){
  rowDiffsCpp(getPos(x)) # matrixStats::rowDiffs(getPos(x)) is a (slower) alternative
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