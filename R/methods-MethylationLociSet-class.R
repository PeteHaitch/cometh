### =========================================================================
### MethylationLociSet: Basically, a GRanges object of all methylation loci
### in a sample, normally the reference genome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' The constructor function for MethylationLociSet objects.
#' 
#' A MethylationLociSet contains the positions of all methylation loci in the sample. This is usually the set of all methylation loci in the reference genome.
#' Most users will construct these objects using the \code{\link{makeMLS}} 
#' function.
MethylationLociSet <- function(seqnames = Rle(), ranges = IRanges(), 
                               strand = Rle("*", length(seqnames)), seqinfo,
                               methylation_type, ...){
  
  ## Check that all required arguments are not missing
  if (missing(seqinfo)){
    stop(sQuote('seqinfo'), " missing.\nPlease see the help page for MethylationLociSet, which can accessed by typing ", sQuote("?MethylationLociSet"), " at the R prompt, for further details of this argument.")
  }
  if (missing(methylation_type)){
    stop(sQuote('methylation_type'), " missing.\nPlease see the help page for MethylationLociSet, which can accessed by typing ", sQuote("?MethylationLociSet"), " at the R prompt, for further details of this argument.")
  }
  
  ## Check that all required arguments are of the correct class 
  if (!is(seqinfo, 'Seqinfo')){
    stop(sQuote('seqinfo'), " must be a ", sQuote('Seqinfo'), " object.\nPlease see the help page for MethylationLociSet, which can accessed by typing ", sQuote("?MethylationLociSet"), " at the R prompt, for further details of this argument.")
  }
  if (!is(methylation_type, 'character')){
    stop(sQuote('methylation_type'), " must be a ", sQuote('character'), ".\nPlease see the help page for MethylationLociSet, which can accessed by typing ", sQuote("?MethylationLociSet"), " at the R prompt, for further details of this argument.")
  }
  
  ## Make GRanges
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, 
                seqinfo = seqinfo)
  
  ## Construct MethylationLociSet object
  new("MethylationLociSet", gr, methylation_type = methylation_type)
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

## TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeMLS
###

## TODO

#' Make the MethylationLociSet from a reference genome.