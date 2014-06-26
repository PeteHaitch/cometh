### =========================================================================
### CoMeth: methylation patterns at m-tuples
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

#' @include AllGenerics.R
#' @rdname CoMeth
setMethod("getPos", "CoMeth", function(x){
  getPos(rowData(x))
  })

#' @include AllGenerics.R
#' @rdname CoMeth
setMethod("getIPD", "CoMeth", function(x){
  getIPD(rowData(x))
})

#' @include AllGenerics.R
setMethod("getMethylationType", "CoMeth", function(x){
  x@colData$methylation_type
})

#' @include AllGenerics.R
#' @export
#' @rdname CoMeth
setMethod("sampleNames", "CoMeth", function(object) {
  colnames(object)
})

## TODO: This might not be properly exported; How do I document a replace 
## method?
# #' @include AllGenerics.R
# setReplaceMethod("sampleNames", 
#                  signature = signature(object = "CoMeth", value = "ANY"), 
#                  function(object, value) {
#   if (length(value) != length(sampleNames(object))){
#     stop("Invalid ", sQuote('sampleNames'), " length")
#   }
#   colnames(object) <- value
#   object
# })

## setMethod("f<-") is equivalent to setReplaceMethod("f")
## See http://stackoverflow.com/questions/24253048/whats-the-difference-between-setmethod-and-set-setreplacemethod?lq=1
#' @include AllGenerics.R
#' @export
#' @rdname CoMeth
setMethod("sampleNames<-", 
                 signature = signature(object = "CoMeth", value = "ANY"), 
                 function(object, value) {
                   if (length(value) != length(sampleNames(object))){
                     stop("Invalid ", sQuote('sampleNames'), " length")
                   }
                   colnames(object) <- value
                   object
                   })

#' @include AllGenerics.R
setMethod("getM", "CoMeth", function(x) {
  getM(rowData(x))
})

setMethod("getCoverage", "CoMeth", function(x) {
  m <- getM(x)
  assay_names <- .make_m_tuple_names(m)
  Reduce("+", assays(x)[assay_names])
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

## cbind combines CoMeth objects from different samples and (possibly) different m-tuples
## TODO: Need to completely re-think how cbind works with new CoMeth constructor.
## TODO: mcols() info of each CoMeth object is lost; fix.
##
# setMethod("cbind", "CoMeth", function(..., deparse.level = 1){
#   
#   args <- unname(list(...))
#   
#   ## Check that all objects have the same seqinfo. 
#   ## Don't want to cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being cbind-ed have the same seqinfo. 
#   ## While it might be nice to to cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from hg18). Rather than handling each case in cbindI will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely cbind-ed
#   seqinfo <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
#   if (is(seqinfo, "try-error")){
#     stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('seqinfo'), ".")
#   }
#   
#   ## Check that all CoMeth objects have the same 'm'
#   m <- sapply(args, getM)
#   if (!.zero_range(m)){
#     stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples, that is, different ", sQuote('m'), ".")
#   }
#   m <- m[1]
#   
#   ## Check that there are no duplicate sample_names
#   sample_names <- unlist(lapply(args, sampleNames))
#   if (any(duplicated(sample_names))){
#     stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")
#   }
#   
#   ## Check that all colData can be safely rbind-ed. Note that even though we are cbind-ing the CoMeth objects, we have to rbind their colData.
#   ## TODO: There might be a better way to do this using GenomicRanges:::.cbind.DataFrame, which looks like it might allow for different columns in the colData. However, I would still need to ensure that all colData contain the 'methylation_type' column.
#   if (!isTRUE(all(sapply(args, function(x, x1){identical(sort(colnames(colData(x))), sort(colnames(colData(x1))))}, x1 = args[[1]])))){
#     stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('colData'), ".") 
#   }
#   
#   ## Extract and combine all data
#   sample_names <- as(unlist(lapply(args, sampleNames)), "CharacterList")
#   methylation_type <- as(unlist(lapply(args, getMethylationType)), "CharacterList")
#   names(methylation_type) <- sample_names
#   # Combine assays
#   
#   
#   counts <- DataFrameList(unlist(lapply(args, function(x, m){
#     lapply(sampleNames(x), function(sn, x, m){
#       #val <- sapply(seq_len(length(assays(x))), function(i, sn, x){
#       val <- sapply(.make_m_tuple_names(m), function(an, sn, x){
#         #assay(x, i)[, sn]
#         assay(x, an)[, sn]
#         }, sn = sn, x = x)
#       colnames(val) <- .make_m_tuple_names(m)
#       val <- DataFrame(val)
#       return(val)}, 
#       x = x, m = m)
#   }, m = m)))
#   names(counts) <- sample_names
#   seqnames <- RleList(unlist(lapply(args, function(x){replicate(n = ncol(x), seqnames(x))})))
#   names(seqnames) <- sample_names
#   pos <- DataFrameList(unlist(lapply(args, function(x){replicate(n = ncol(x), DataFrame(getPos(x)))})))
#   names(pos) <- sample_names
#   strand <- RleList(unlist(lapply(args, function(x){replicate(n = ncol(x), strand(x))})))
#   names(strand) <- sample_names
#   colData <- do.call(what = "rbind", args = lapply(args, function(x){
#     colData(x)[-c(which(colnames(colData(x)) == 'methylation_type'))]
#   }))
#   exptData <- do.call(what = "c", args = lapply(args, exptData))
#   
#   ## Construct CoMeth object from combined data
#   CoMeth(sample_names = sample_names, methylation_type = methylation_type, counts = counts, seqnames = seqnames, pos = pos, seqinfo = seqinfo, strand = strand, colData = colData, exptData = exptData)
# })

## rbind combines Cometh objects from the same samples and different m-tuples.
## rbind won't work on CoMeth objects with different mcols
## TODO: Need to completely re-think how rbind works with the new CoMeth constructor.
##
# setMethod("rbind", "CoMeth", function(..., deparse.level = 1){
# 
#   args <- unname(list(...))
#   
#   ## Check that all objects have the same seqinfo. 
#   ## Don't want to rbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed have the same seqinfo. 
#   ## While it might be nice to to rbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from hg18). Rather than handling each case in rbindI will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed
#   seqinfo <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
#   if (is(seqinfo, "try-error")){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('seqinfo'), ".")
#   }
#   
#   ## Check that all CoMeth objects have the same 'm'
#   m <- sapply(args, getM)
#   if (!.zero_range(m)){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples (i.e. different ", sQuote('m'), ").")
#   }
#   m <- m[1L]
#   
#   ## Check that all CoMeth objects have the same sample_names
#   if (!isTRUE(all(sapply(args, function(x, x1){identical(sampleNames(x), sampleNames(x1))}, x1 = args[[1]])))){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('sample_names'), ".")
#   }
#   
#   ## Check that all mcols of rowData are idential
#   if (!isTRUE(all(sapply(args, function(x, x1){identical(colnames(mcols(x)), colnames(mcols(x1)))}, x1 = args[[1]])))){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('mcols'), ".")
#   }
#   
#   ## Check that all colData is identical
#   if (!isTRUE(all(sapply(args, function(x, x1){identical(colData(x), colData(x1))}, x1 = args[[1]])))){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('colData'), ".") 
#   }
#   
#   ## Check that there are no identical m-tuples
#   if(sapply(seq_len(length(args) - 1), function(i, args){any(overlapsAny(args[[i]], args[[i + 1]], type = 'equal'))}, args = args)){
#     stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")
#   }
#   
#   ## Extract and combine all data (NEW WAY)
#   sample_names <- unique(unlist(lapply(args, sampleNames)))
#   methylation_type <- paste0(sort(unique(unlist(lapply(args, getMethylationType)))), collapse = '/')
#   colData <- colData(args[[1]]) # colData has already been checked and found to be identical for all args
#   mtuples <- do.call("c", lapply(args, rowData))
#   assays <- SimpleList(lapply(names(assays(args[[1]])), function(i, args){
#     do.call("rbind", lapply(args, function(y, i){
#       assay(y, i)
#     }, i = i))
#   }, args = args))
#   names(assays) <- names(assays(args[[1]]))
#   exptData <- do.call(what = "c", args = lapply(args, exptData))
#   
#   ## Don't call the CoMeth constructor but because there's no need to "combine"
#   ## data, just need to initialise the class.
#   if (m == 1L){
#     class <- "CoMeth1"
#   } else if (m == 2L){
#     class <- "CoMeth2"
#   } else{
#     class <- "CoMeth3Plus"
#   }
#   
#   new(class, SummarizedExperiment(assays = assays, rowData = mtuples, colData = colData, exptData = exptData))
# })
# 


# TODO: Need to re-think in light of changes to CoMeth constructor.
# ## combine tries to figure out the combination of cbind and rbind that will automatically combine the CoMeth objects into a single CoMeth object.
# ## At this stage, only simple combinations are covered and more complex combinations will require manual handling.
# ## I'll try to add more advanced combining as the need arises.
# setMethod("combine", "CoMeth", function(x, y, ...){
#   
#   if (class(x) != class(y)){
#     stop(paste("objects must be the same class, but are ", class(x), ", ", class(y), sep = ""))
#   }
#   
#   if (!identical(getM(x), getM(y))){
#     stop("Cannot ", sQuote('combine'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples (i.e. different ", sQuote('m'), ").")
#   }
#   
#   ## Try to rbind if objects contain the same sample names
#   if (identical(sampleNames(x), sampleNames(y))){
#     val <- try(rbind(x, y), silent = TRUE)
#     if (is(val, "try-error")){
#       stop(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
#     }
#   ## Try to cbind if objects contain different sample names
#   } else if (!isTRUE(any(duplicated(c(sampleNames(x), sampleNames(y)))))){
#     val <- try(cbind(x, y), silent = TRUE)
#     if (is(val, "try-error")){
#       stop("Sorry, ", sQuote('combine'), " failed when trying to ", sQuote('cbind'), " the intermediate ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
#     }
#   ## Otherwise we can't combine them automatically :(
#   } else if (isTRUE(any(sampleNames(x) %in% sampleNames(y)))){
#     stop("Sorry, cannot automatically ", sQuote('combine'), " this combination of ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
#   }
#   return(val)
# })

## TODO: Add warnings if assays slot contains non-standard assays that these 
## will be removed.
## TODO: Should this be a function or a method?
## TODO: Figure out where bplapply can be effectively (and safely) subbed in 
## for lapply.
#' Combine counts across strand.
#' 
#' @description 
#' Combine the counts (e.g. \code{M} and \code{U} in a \code{CoMeth1} object) 
#' across strands. Only applicable if the methylation type is CG, which is 
#' generally strand-symmetric, because all other methylation types are 
#' strand-asymmetric. 
#' 
#' @details
#' This will only combine the counts and \strong{all other assays will be 
#' removed}; the exceptions are for \code{CoMeth1} and \code{CoMeth2} objects, 
#' which will have \code{beta} and \code{LOR} re-computed, respectively. The 
#' positions of all m-tuples will be with respect to the position of the 
#' cytosine on the '+' strand, however, the strand of all m-tuples will be set 
#' to \code{*}.
#' 
#' @param x A \code{CoMeth} object with \code{methylation_type == 'CG'}.
#' @param sort A \code{logical(1)} indicating whether the resulting object 
#' should be sorted (default: \code{TRUE}).
#' 
#' @return A \code{CoMeth} object of the same class as \code{x}. See details.
#' 
#' @export
#' @rdname CoMeth
combineStrands <- function(x, sort = TRUE) {
  
  if (!inherits(x = x, what = "CoMeth")) {
    stop(sQuote('x'), " must inherit from a ", sQuote('CoMeth'), " object.")
  }
  
  if (!all(getMethylationType(x) == 'CG')) {
    stop("Can only combine counts across strands if the ", 
         sQuote('methylation_type'), " of all samples is ", sQuote('CG'), ".")
  }
  
  if (any(strand(x) == '*')) {
    stop("Some of the counts have already been combined across strands (i.e. ", 
         sQuote('any(strand(x) == '*')'), ' is TRUE).')
  }
  
  # Split the rowData based on strand. 
  # Throw away "*"-strand since there aren't any m-tuples on that strand.
  ## TODO: split produces a GRangesList not an MTuplesList (which isn't yet 
  ## defined). The main problem with this is that y doesn't "look" like a List 
  ## of MTuples objects because it defers to the GRangesList "show" method.
  print("Combining...")
  y <- split(rowData(x), strand(rowData(x)))[1:2]
  ol <- findOverlaps(y[[1]], shift(unstrand(y[[2]]), shift = -1L, 
                                   use.names = FALSE), 
                     type = 'equal')
  plus_ol <- findOverlaps(y[[1]], rowData(x), type = 'equal')
  neg_ol <- findOverlaps(y[[2]], rowData(x), type = 'equal')
  # plus_only and neg_only are with respect to x.
  plus_only <- subjectHits(plus_ol)[countQueryHits(ol) == 0]
  neg_only <- subjectHits(neg_ol)[countSubjectHits(ol) == 0]
  # Order is both, plus_only, neg_only.
  # neg_only m-tuples need to be shifted -1L (wrt '+' strand).
  rowData <- unstrand(c(y[[1]][queryHits(ol)], y[[1]][countQueryHits(ol) == 0], 
          shift(y[[2]][countSubjectHits(ol) == 0], shift = -1L)))
  ci <- grep(pattern = '^[MU]', x = names(assays(x)), value = TRUE)
  ## TODO: Is parallelisation worth it?
  both_assays <- lapply(assays(x)[ci], function(counts, ol) {
    counts[queryHits(ol), , drop = FALSE] + 
      counts[subjectHits(ol), , drop = FALSE]
  }, ol = ol)
  plus_assays <- lapply(assays(x)[ci], function(counts, plus_only) {
    counts[plus_only, , drop = FALSE]
  }, plus_only = plus_only)
  neg_assays <- lapply(assays(x)[ci], function(counts, neg_only) {
    counts[neg_only, , drop = FALSE]
  }, neg_only = neg_only)
  assays <- mapply(function(b, p, n) {
    rbind(b, p, n)
  }, b = both_assays, p = plus_assays, n = neg_assays, SIMPLIFY = FALSE)
  
  m <- getM(rowData)
  if (m == 1L) {
    class <- 'CoMeth1'
    assays <- c(assays, list(beta = .beta(assays)))
  } else if (m == 2L) {
    class <- 'CoMeth2'
    assays <- c(assays, list(LOR = .LOR(assays)))
  } else{
    class <- 'CoMeth3Plus'
    # No need to update assays
  }
  
  # Create the CoMeth object
  print("Creating CoMeth object...")
  z <- new(class, SummarizedExperiment(assays = assays, rowData = rowData, 
                                       colData = colData(x), 
                                       exptData = exptData(x), 
                                       verbose = FALSE))
  if (sort) {
    z <- sort(z)
  }
  return(z)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison
###

## TODO: Check the other comparison methods work (i.e. those defined for MTuples)
## compare defers to the method defined for MTuples.
#' Compare \code{MTuples}.
#' @export
#' @rdname CoMeth
setMethod("compare", c("CoMeth", "CoMeth"), function(x, y) {
  
  if (getM(x) != getM(y)){
    stop("Cannot ", sQuote('compare'), " ", sQuote('CoMeth'), 
         " objects with different ", sQuote('m'))
  }
  compare(rowData(x), rowData(y))
})

#' @export
#' @rdname CoMeth
setMethod("order", "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
  
  args <- list(...)
  
  if (!.zero_range(sapply(args, getM))){
    stop("All ", sQuote('CoMeth'), " objects must have the same ", sQuote('m'), 
         " value")
  }
  
  args <- lapply(list(...), rowData)
  do.call("order", c(args, list(na.last = na.last, decreasing = decreasing)))  
})

## TODO: Add an is.unsorted method.
## TODO: sort(CoMeth, ignore.strand = TRUE) doesn't seem to be supported.
  


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Filtering
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

# TODO: Function collapseStrand(). Should only work with symmetric methylation_type, e.g. CG, CHG. Will be hard to do for more complication methylation_type, e.g. CG/CHG, without using the reference genome (slow).


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

## TODO: Include methylation_type and 'm' when show-ing a CoMeth object.
## TODO: Rename "methylation patterns", as "assays" (or something similar).

## The show method is adapted from that of SummarizedExperiment
#' @export
#' @rdname CoMeth
setMethod("show", "CoMeth", function(object){ 
  selectSome <- BiocGenerics:::selectSome
  scat <- function(fmt, vals = character(), exdent = 2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(BiocGenerics:::selectSome(vals), collapse = " ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
  }
  cat("class:", class(object), "\n")
  cat("dim:", dim(object), "\n")
  expt <- names(exptData(object))
  if (is.null(expt)) 
    expt <- character(length(exptData(object)))
  scat("exptData(%d): %s\n", expt)
  nms <- names(assays(object, withDimnames = FALSE))
  if (is.null(nms)) 
    nms <- character(length(assays(object, withDimnames = FALSE)))
  scat("methylation patterns(%d): %s\n", nms)
  dimnames <- dimnames(object)
  dlen <- sapply(dimnames, length)
  if (dlen[[1]]) 
    scat("rownames(%d): %s\n", dimnames[[1]])
  else scat("rownames: NULL\n")
  scat("rowData metadata column names(%d): %s\n", names(mcols(rowData(object))))
  if (dlen[[2]]) 
    scat("sampleNames(%d): %s\n", dimnames[[2]])
  else cat("colnames: NULL\n")
  scat("colData names(%d): %s\n", names(colData(object)))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Intra range methods
###

## TODO: Check that these can all be safely inherited from SummarizedExperiment.