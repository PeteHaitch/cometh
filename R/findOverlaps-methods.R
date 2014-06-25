### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Finding overlaps
###

## TODO: Replace .local() within setMethod calls with a separate (better named) function and then call that function from setMethod
## TODO: What happens if findOverlaps (type = 'equal') with m-tuples (m > 2) when positions are identical but strands differ, e.g. '+' vs. '+' (should be a match) or '+' vs. '-' (should not be a match).
## TODO: Define an overlapsAny method for MTuples - the GRanges method won't 
## work as intended because it uses "select = 'first'" when really we want 
## overlapsAny(MTuples, MTuples) to be identical to  
## countOverlaps(MTuples, MTuples, type = 'equal', select = 'all') != 0. Need to
## think about what we want if type != 'equal'.
## TODO: Really want a better method when m > 2; specifically, it should be 
## possible to hash (chr:strand:pos1:...:posm) and do a fast match rather than 
## using multiple calls to findOverlaps method for GRanges, which includes the 
## overhead of constructing multiple GRanges.

## AWAITING RESPONSE FROM BIOC-DEVEL (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005549.html) THAT IS RELEVANT TO THESE METHOD DEFINITIONS
## REPLY FROM BIOC-DEVEL (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005555.html): Can remove this method once I have updated to GenomicRanges >= 1.15.46 (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005555.html).
## TODO: Once this change is made, update DESCRIPTION with requiriment of GenomicRanges >= 1.15.46.
## TODO: Refinement for ciruclar chromosomes

## There is a specially defined method for findOverlaps when both the query and the subject are MTuples objects.
## This is to allow for "exact" matching between MTuples.
## If either the subject or the query is not a MTuples object then it defers to the findOverlaps method defined for GRanges objects. 

## Similarly, there is a specially defined method for findOverlaps when either the query or the subject or both are CoMeth objects.
## I cannot simply defer to the findOverlaps method defined for SummarizedExperiment because it does not allow for 'type = equal' nor 'select = last' or 'select = arbitrary' (at least in GenomicRanges_1.14.4).

setMethod("findOverlaps", 
          signature = c("MTuples", "MTuples"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...){
            
            ## Require special handling if type == 'equal' and m > 2, otherwise defer to GRanges method
            .local <- function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...){
              
              ## Check Tuples are compatible (i.e. have the same m)
              q_m <- getM(query)
              s_m <- getM(subject)
              if (q_m != s_m){
                stop("Cannot ", sQuote("findOverlaps"), " between ", sQuote(class(query)), " and ", sQuote(class(subject)), " if they have different ", sQuote("m"), ".")
              }
              m <- q_m
              
              ## Argument matching
              select <- match.arg(select)
              type <- match.arg(type)
              
              ## If type isn't 'equal' then just use the findOverlaps method defined for GRanges objects
              if (type != 'equal'){
                callNextMethod()
              } else{
              
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
                  warning("No special handling of circular chromosomes.")
                }
              
                ## If type is 'equal' and m >= 3 then use the internal function .findIdentical.MTuples().
                ## Else, if type is 'equal' and m = 1 or m = 2 then just defer to the findOverlaps method for GRanges
                if (m >= 3){
                  .findIdentical.MTuples(query, subject, select)
                } else {
                  callNextMethod()
                }
              }
            }
            .local(query, subject, maxgap, minoverlap, type, select, ...)
          }
)

## TODO: Remove this method once I have updated to GenomicRanges version >= 1.15.46 (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005555.html)

setMethod("findOverlaps", 
          signature = c("CoMeth", "CoMeth"),
          function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ...) 
          {
            .local <- function(query, subject, maxgap = 0L, minoverlap = 1L, type = c("equal", "any", "start", "end", "within"), select = c("all", "first", "last", "arbitrary"), ignore.strand = FALSE) 
            {
              findOverlaps(rowData(query), rowData(subject), maxgap = maxgap, 
                           minoverlap = minoverlap, type = match.arg(type), 
                           select = match.arg(select), ...)
            }
            .local(query, subject, maxgap, minoverlap, type, select, ...)
          }
)

## TODO: Remove this method once I have updated to GenomicRanges version >= 1.15.46 (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005555.html)

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

## TODO: Remove this method once I have updated to GenomicRanges version >= 1.15.46 (https://stat.ethz.ch/pipermail/bioc-devel/2014-April/005555.html)

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