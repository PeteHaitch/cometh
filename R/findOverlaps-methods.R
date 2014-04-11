### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Finding overlaps
###

## AWAITING RESPONSE FROM BIOC-DEVEL (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005549.html) THAT IS RELEVANT TO THESE METHOD DEFINITIONS

## There is a specially defined method for findOverlaps when both the query and the subject are MTuples objects.
## This is to allow for "exact" matching between MTuples.
## If either the subject or the query is not a MTuples object then it defers to the findOverlaps method defined for GRanges objects. 

## Similarly, there is a specially defined method for findOverlaps when either the query or the subject or both are CoMeth objects.
## I cannot simply defer to the findOverlaps method defined for SummarizedExperiment because it does not allow for 'type = equal' nor 'select = last' or 'select = arbitrary' (at least in GenomicRanges_1.14.4).

setMethod("findOverlaps", 
          signature = c("MTuples", "MTuples"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...){
            
            ## Require special handling if type == 'equal' and m > 2, otherwise defer to GRanges method
            .local <- function(query, subject, maxgap = 0L, minoverlaps = 1L, type = c("any", "start", "end", "within", "equal"), select = c("all", "first", "last", "arbitrary"), ignore.strand = FALSE){
              
              ## Check Tuples are compatible (i.e. have the same m)
              q_m <- getM(query)
              s_m <- getM(subject)
              if (q_m != s_m){
                stop("Cannot ", sQuote("findOverlaps"), " between ", sQuote(class(query)), " and ", sQuote(class(subject)), " if they have different ", sQuote("m"))
              }
              
              ## Argument matching
              select <- match.arg(select)
              type <- match.arg(type)
              
              ## If type isn't 'equal' then just use the findOverlaps method defined for GRanges objects
              if (type != 'equal'){
                callNextMethod()
              }
              
              ## The 'maxgap' and 'minoverlap' parameters aren't used when type is 'equal'.
              ## Report a warning if these have been set by the user (actually, if they differ from their defaults)
              if (maxgap != 0L || minoverlap != 1L){
                warning(sQuote("maxgap"), " and ", sQuote("minoverlap"), " ignored when ", sQuote("type = equal"), ".")
                maxmap <- 0L
                minoverlap <- 1L
              }
              
              ## The internal function .findIdentical.MTuples() hasn't been tested with circular chromosomes.
              seqinfo <- merge(seqinfo(query), seqinfo(subject))
              if (isTRUE(any(isCircular(seqinfo)))) {
                stop("Cannot handle circular chromosomes.")
              }
              
              ## If type is 'equal' then use the internal function .findIdentical.MTuples().
              hits <- .findIdentical.MTuples(query, subject, select)
              
              return(hits)
            }
            
            .local(query, subject, maxgap, minoverlaps, type, select, ...)
          }
)

setMethod("findOverlaps", 
          signature = c("CoMeth", "CoMeth"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...) 
          {
            .local <- function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ignore.strand = FALSE) 
            {
              findOverlaps(rowData(query), rowData(subject), maxgap = maxgap, 
                           minoverlap = minoverlap, type = match.arg(type), 
                           select = match.arg(select), ignore.strand = ignore.strand)
            }
            .local(query, subject, maxgap, minoverlap, type, select, ...)
          }
)

setMethod("findOverlaps", 
          signature = c("CoMeth", "MTuples"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...) 
          {
            .local <- function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ignore.strand = FALSE) 
            {
              findOverlaps(rowData(query), subject, maxgap = maxgap, 
                           minoverlap = minoverlap, type = match.arg(type), 
                           select = match.arg(select), ignore.strand = ignore.strand)
            }
            .local(query, subject, maxgap, minoverlap, type, select, ...)
          }
)

setMethod("findOverlaps", 
          signature = c("MTuples", "CoMeth"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...) 
          {
            .local <- function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ignore.strand = FALSE) 
            {
              findOverlaps(query, rowData(subject), maxgap = maxgap, 
                           minoverlap = minoverlap, type = match.arg(type), 
                           select = match.arg(select), ignore.strand = ignore.strand)
            }
            .local(query, subject, maxgap, minoverlap, type, select, ...)
          }
)