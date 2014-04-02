### =========================================================================
### Tuple objects
### -------------------------------------------------------------------------

## A Tuple object stores the genomic co-ordinates of tuples, which are sets of single positions on the same chromosome. 
## A tuple consists of a seqname, a strand, and the positions.
## All co-ordinates are 1-based.
## For example, chr1:+:(3, 7, 9) is a text-representation of tuple containing the three points (3, 7, 9) on the forward strand of chromosome 1. 
## Tuples are different to ranges (e.g. GRanges or IRanges) in that any point not explicity listed is not included in the tuple (e.g. chr1:+:4, chr1:+:5 or chr1:-:3 in the above example).

## The class definition is adapted from Martin Morgan's reply to my question on Bioc-Devel (https://stat.ethz.ch/pipermail/bioc-devel/2014-February/005242.html) and Michael Lawerence's reply to my question on BioC-Help (https://stat.ethz.ch/pipermail/bioconductor/2014-March/058490.html)

.Tuple <- setClass("Tuple", contains = "DataFrame", representation(seqinfo = "Seqinfo", m = "integer"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("nrow", "Tuple", function(x){
  callNextMethod()
})

## NOTE: It might make sense to define ncol as "ncol - 2" (i.e. x@m).
setMethod("ncol", "Tuple", function(x){
  callNextMethod()
})

## Doesn't work and isn't needed since it inherits from DataFrame in any case
#setMethod("rownames", "Tuple", function(x, do.null = TRUE, prefix = "row"){
#  callNextMethod()
#})

## Doesn't work and isn't needed since it inherits from DataFrame in any case
#setMethod("colnames", "Tuple", function(x, do.null = TRUE, prefix = "col"){
#  callNextMethod()
#})

setReplaceMethod("rownames", "Tuple", function(x, value){
  callNextMethod()
})

setReplaceMethod("colnames", "Tuple", function(object, value){
  stop("Replacing ", sQuote('colnames'), " of ", sQuote(class(object))," object is not supported.")
})

setMethod("dim", "Tuple", function(x){
  c(nrow(x), ncol(x))
})

setReplaceMethod("dim", "Tuple", function(object, value){
  stop("Replacing ", sQuote('dim'), " of ", sQuote(class(object)), " is not supported.")
})

setMethod("dimnames", "Tuple", function(x){
  callNextMethod()
  })

setReplaceMethod("dimnames", "Tuple", function(object, value){
  stop("Replacing ", sQuote('dimnames'), " of ", sQuote(class(object))," object is not supported because replacing ", sQuote('colnames'), " of ", sQuote(class(object))," object is not supported. If you wish to just replace the ", sQuote('rownames'), " of ", sQuote(class(object)), " please use the ", sQuote('rownames'), " replacement function.")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

validTuple <- function(object){
  
  msg <- NULL
  
  ## Check 'm'
  if (object@m != (ncol(object) - 2)){
    msg <- validMsg(msg, paste0("Invalid ", sQuote('m'), ". Should have that ", sQuote('m'), " + 2 = ", sQuote('ncol(Tuple)'), " but ", sQuote('m'), " = ", object@m, " and ", sQuote('ncol(Tuple)'), " = ", ncol(object)))
  }
  
  ## Check colnames
  if (!identical(colnames(object), c('seqnames', 'strand', paste0('pos', seq_len(object@m))))){
    msg <- validMsg(msg, paste0("Invalid ", sQuote('colnames(Tuple)'), ". Should be: ", paste0(sapply(c('seqnames', 'strand', paste0('pos', seq_len(object@m))), sQuote, USE.NAMES = FALSE), collapse = ", "), "."))
  }
  
  ## Check pos1, ..., posm are positive integers
  pos_positive_integer <- isTRUE(all(sapply(seq(from = 3, to = ncol(object), by = 1), FUN = function(i, object){
    is.integer(object[[i]]) && isTRUE(all(object[[i]] > 0))
  }, object = object)))
  if (!pos_positive_integer){
    msg <- validMsg(msg, paste0(c("All elements of ", paste0(sapply(paste0('pos', seq_len(object@m - 1)), sQuote, USE.NAMES = FALSE), sep = ", "), sQuote(paste0('pos', object@m)), " must be positive integers."), collapse = ""))
  }
  
  ## Check pos1, ..., posm columns are ordered left-to-right
  pos_ordered <- isTRUE(all(sapply(seq(from = 3, to = ncol(object) - 1, by = 1), FUN = function(i, object){
    isTRUE(all(object[[i]] < object[[(i + 1)]]))
  }, object = object)))
  if (!pos_ordered){
    msg <- validMsg(msg, paste0(c("Must have ", paste0(sapply(paste0('pos', seq_len(object@m - 1)), sQuote, USE.NAMES = FALSE), sep = " < " ), sQuote(paste0('pos', object@m)), " for each ", object@m, "-tuple."), collapse = ''))
  }
  
  ## Check all seqnames have an entry in 'seqinfo'
  if (any(!levels(object[['seqnames']]) %in% seqlevels(object@seqinfo))){
    msg <- validMsg(msg, paste0("All values of the ", sQuote('seqnames'), " column must have a corresponding entry in the ", sQuote('seqinfo'), " field."))
  }
  
  ## Check seqlengths of seqinfo are valid, i.e. greater than the biggest position for that chromosome (unless the seqlength is NA)
  max_pos <- sapply(X = levels(object[[1]]), FUN = function(seq_name, object){
    x <- object[object[['seqnames']] == seq_name, ]
    max(sapply(seq(from = 3, to = ncol(object)), FUN = function(i, x){max(x[[i]])}, x = x))}, object = object)
  if (isTRUE(any(sapply(names(max_pos), FUN = function(i, max_pos, seqinfo){
    max_pos[i] > seqlengths(seqinfo)[i]
  }, max_pos = max_pos, seqinfo = object@seqinfo, USE.NAMES = FALSE)))){
    msg <- validMsg(msg, paste0("Invalid ", sQuote('seqinfo'), ". ", sQuote('seqlengths(seqinfo)'), " for each ", sQuote('seqlevel'), " must be greater than the largest position for tuples on that ", sQuote('seqlevel'), "."))
  }
    
  ## Return validity of object
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}
setValidity("Tuple", validTuple)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

Tuple <- function(pos, seqinfo){
  
  ## Check colnames
  if (!identical(colnames(pos), c('seqnames', 'strand', paste0('pos', seq_len(ncol(pos) - 2))))){
    stop("Invalid ", sQuote('colnames(pos)'), ". Should be: ", paste0(sapply(c('seqnames', 'strand', paste0('pos', seq_len(object@m))), sQuote, USE.NAMES = FALSE), collapse = ", "), ".")
  }
  
  ## Check seqinfo is a Seqinfo object
  if (!is(seqinfo, "Seqinfo")){
    stop(sQuote('seqinfo'), " must be a ", sQuote('Seqinfo'), " object")
  }
  
  ## Convert seqnames and strand to Factor-Rle objects if they aren't already
  pos[["seqnames"]] <- Rle(factor(as.character(pos[["seqnames"]]), levels = seqlevels(seqinfo), ordered = TRUE))
  pos[["strand"]] <- Rle(factor(as.character(pos[["strand"]]), levels = c("+", "-", "*"), ordered = TRUE))
  ## note constructor: base class(es) as first and unnamed arg
  .Tuple(pos, seqinfo = seqinfo, m = as.integer(ncol(pos) - 2))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod(show, signature(object = "Tuple"), function(object) {
  # Modified version of the show method for DataTable
  nhead <- get_showHeadLines()
  ntail <- get_showTailLines()
  nr <- nrow(object)
  m <- object@m
  cat(class(object), " with ", nr, 
      ifelse(nr == 1, paste0(" x ", m, "-tuple"), paste0(" x ", m, "-tuples")), 
      sep = "")
  cat('\n---\n')
  if (nr > 0 && m > 0) {
    nms <- rownames(object)
    if (nr < (nhead + ntail + 1L)) {
      out <- as.matrix(format(as.data.frame(lapply(object, 
                                                   showAsCell), optional = TRUE)))
      if (!is.null(nms)) 
        rownames(out) <- nms
    }
    else {
      out <- rbind(as.matrix(format(as.data.frame(lapply(object, function(x) showAsCell(head(x, nhead))), optional = TRUE))), 
                   rbind(rep.int("...", m + 2)), as.matrix(format(as.data.frame(lapply(object, function(x){showAsCell(tail(x, ntail))}), optional = TRUE))))
      rownames(out) <- IRanges:::.rownames(nms, nr, nhead, ntail)
    }
    classinfo <- matrix(unlist(lapply(object, function(x) {
      paste0("<", classNameForDisplay(x)[1], ">")
    }), use.names = FALSE), nrow = 1, dimnames = list("", colnames(out)))
    out <- rbind(classinfo, out)
    print(out, quote = FALSE, right = TRUE)
  }
  cat('---\n')
  print(as.data.frame(object@seqinfo))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "Tuple",
          function(x, i, j, ..., drop = FALSE){
            if (!missing(drop)){
              warning(sQuote('drop'), " ignored when subsetting ", sQuote(class(x)))
            }
            if (!missing(j)){
              stop("Cannot subset by j when subsetting ", sQuote(class(x)))
            }
            if (length(list(...)) > 0L){
              stop("parameters in ", sQuote('...'), " not supported")
            }
            ## Use the "[" method defined for DataFrames
            callNextMethod()
          }
)

## "[[" is inherited from DataFrame. Don't define it via callNextMethod() as it breaks things for reasons I do not understand.

setReplaceMethod("[", "Tuple", 
                 function(x, value){
                   stop("Replacement of ", sQuote(class(object)), " values via ", sQuote('['), " is not supported.")
                 }
)

setReplaceMethod("[[", "Tuple", 
                 function(x, i, j, value){
                   stop("Replacement of ", sQuote(class(x)), " columns via ", sQuote('[['), " is not supported.")
                   }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison.
###

setMethod("duplicated", "Tuple", function(x, incomparables = FALSE, ...){
  callNextMethod()
})

## FIXME: Generates a warning due to how subset method for Tuple class
setMethod("unique", "Tuple", function(x, incomparables = FALSE, ...){
  callNextMethod()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

## as.list, as.data.frame, as(from, "data.frame") all inherited from DataFrame

setMethod("as.matrix", "Tuple", function(x){
  ## Can't go directly to as.matrix because Rle-factor columns (seqnames, strand) aren't properly converted
  as.matrix(as.data.frame(x))
})
