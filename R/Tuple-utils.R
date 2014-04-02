### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ordering and Sorting
###

## The order method for a Tuple object; orders the rows of a Tuple object.
## Each row of a Tuple is column-wise ordered left-to-right.
## 'seqnames' is ordered with respect to 'seqlevels(seqinfo)'.
## 'strand' is ordered + < - < * (i.e. the same as GRanges objects)
## So, for example, chr1:+:(3, 5, 11) < chr1:-:(3, 5, 11) < chr1:+:(3, 5, 11) < chr1:(3, 7, 9) < chr1:(3, 7, 11).
## WARNING: I don't really understand when you would use the case where "..." contains multiple objects but I've included it anyway based on selectMethod("order", "GenomicRanges").
setMethod("order", "Tuple", function(..., na.last = TRUE, decreasing = FALSE){
  
  if (!isTRUEorFALSE(decreasing)){
    stop(sQuote('decreasing'), " must be TRUE or FALSE")
  }
  
  args <- list(...)
  
  if (!zero_range(sapply(args, function(x){x@m}))){
    stop("All ", sQuote('Tuple'), " objects must have the same ", sQuote('m'), ".")
  }
  
  do.call(order, as.data.frame(pos))
})

setMethod(sort, "Tuple", function(x, decreasing = FALSE, ...){
  x[order(x, decreasing = decreasing), ] 
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting and Combining
###

setMethod("cbind", "Tuple", function(..., deparse.level = 1){
  stop("Cannot cbind ", sQuote('Tuple'), " objects")
})

setMethod("rbind", "Tuple", function(..., deparse.level = 1){
  
  args <- list(...)
  
  ## Ensure that all Tuple objects have the same seqinfo
  ## Don't want to rbind Tuple objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same Tuple object provided all Tuple objects that are being rbind-ed have the same seqinfo. 
  ## While it might be nice to to rbind Tuple objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2 from the same genome), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from mm9). Rather than handling each case in rbind I will require to the user to update the seqinfo of the Tuple in an appropriate way to allow for them to be safely rbind-ed.
  seq_info <- try(do.call("merge", lapply(args, function(x){x@seqinfo})), silent = TRUE)
  if (is(seq_info, "try-error")){
    stop("Can only rbind ", sQuote('Tuple')," objects with identical ", sQuote('seqinfo'), ".")
  }
  
  ## Ensure all objects have the same 'm'
  if (!zero_range(sapply(args, function(x){x@m}))){
    stop("Cannot rbind ", sQuote('Tuple'), " objects with the different ", sQuote('m'), ".")
  }
  
  ## Use the rbind method defined for DataFrame objects
  ## I can't simply use callNextMethod() because this ignores the 'm' and 'seqinfo' slots of a Tuple object.
  pos <- do.call(rbind, lapply(args, function(x){as(x, "DataFrame")}))
  
  ## Add back 'm' and 'seqinfo'
  Tuple(pos, seq_info)
})

## split produces a CompressedSplitDataFrameList of Tuple objects by inheritance magic

## TODO: mstack. Not quite sure what it should do.

## TODO: aggregate. Not quite sure what it should do.

