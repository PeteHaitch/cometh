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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison
###

## TODO: Check the other comparison methods work (i.e. those defined for MTuples)
## compare defers to the method defined for MTuples.
#' Compare \code{MTuples}.
#' @export
#' @rdname CoMeth
setMethod("compare", c("CoMeth", "CoMeth"), function(x, y){
  
  if (getM(x) != getM(y)){
    stop("Cannot ", sQuote('compare'), " ", sQuote('CoMeth'), " objects with different ", sQuote('m'))
  }
  compare(rowData(x), rowData(y))
})

#' @export
#' @rdname CoMeth
setMethod("order", "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
  
  args <- list(...)
  
  if (!.zero_range(sapply(args, getM))){
    stop("All ", sQuote('CoMeth'), " objects must have the same ", sQuote('m'), " value")
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