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
#' A MethylationLociSet contains the positions of all methylation loci in the 
#' sample. This is usually the set of all methylation loci in the reference 
#' genome.
#' Most users will construct these objects using the \code{\link{makeMLS}} 
#' function.
MethylationLociSet <- function(seqnames = Rle(), ranges = IRanges(), 
                               strand = Rle("*", length(seqnames)), seqinfo,
                               methylation_type, ...){
  
  # Check that all required arguments are not missing
  if (missing(seqinfo)){
    stop(sQuote('seqinfo'), " missing.\nPlease see the help page for ", 
         "MethylationLociSet, which can accessed by typing ", 
         sQuote("?MethylationLociSet"), 
         " at the R prompt, for further details of this argument.")
  }
  if (missing(methylation_type)){
    stop(sQuote('methylation_type'), " missing.\nPlease see the help page for ",
         "MethylationLociSet, which can accessed by typing ", 
         sQuote("?MethylationLociSet"), 
         " at the R prompt, for further details of this argument.")
  }
  
  # Check that all required arguments are of the correct class 
  if (!is(seqinfo, 'Seqinfo')){
    stop(sQuote('seqinfo'), " must be a ", sQuote('Seqinfo'), 
         " object.\nPlease see the help page for MethylationLociSet, which ", 
         "can accessed by typing ", sQuote("?MethylationLociSet"), " at the R ",
         "prompt, for further details of this argument.")
  }
  if (!is(methylation_type, 'character')){
    stop(sQuote('methylation_type'), " must be a ", sQuote('character'), 
         ".\nPlease see the help page for MethylationLociSet, which can ", 
         "accessed by ", "typing ", sQuote("?MethylationLociSet"), 
         " at the R prompt, for further details of this argument.")
  }
  
  # Check that all ranges are of width 1.
  if (!.zero_range(width(ranges)) || !identical(width(ranges[1]), 1L)){
    stop("All ", sQuote("ranges"), " must have ", sQuote("width"), " = 1.\nThe", 
         " start of each range should reflect the position of the cytosine on", 
         " the relevant strand.")
  }
  # Check that all methylation loci have strand information
  if (isTRUE(any(strand == '*'))){
    stop(sQuote('strand'), " must be ", sQuote('+'), " or ", sQuote('-'), " (", 
         sQuote('*'), " is not permitted).")
  }
  # Check that methylation_type is valid
  if (!all(sapply(X = methylation_type, .valid_methylation_type))){
    stop("Invalid ", sQuote('methylation_type'), 
         ".\nPlease see the help page for MethylationLociSet, which can ", 
         "accessed by ", "typing ", sQuote("?MethylationLociSet"), 
         " at the R prompt, for further details of this argument.")
  }
  
  # Make GRanges
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, 
                seqinfo = seqinfo)
  
  # Construct MethylationLociSet object
  # Sort-unique methylation_type to remove any duplicate elements 
  # (which shouldn't happen).
  new("MethylationLociSet", gr, 
      methylation_type = sort(unique(methylation_type)))
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

## TODO: Define c() such that the union of the methylation types is taken

setMethod("show", "MethylationLociSet", function(object) {
  cat("Methylation type:", paste(sort(object@methylation_type), collapse = "/"),
      "\n")
  callNextMethod()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeMLS
###

#' Make a MethylationLociSet for a methylation type and a reference genome.
#' 
#' Construct an object to store the location of all methylation loci of the 
#' given type in a genome sequence. For example, find all CpGs in the human 
#' reference genome (hg19). The returned object is a 
#' \code{\link{MethylationLociSet}}.
#' 
#' @param methylation_type A character. The possible values are "CG", "CHG", 
#' "CHH" or "CNN" or multiple values specified as a character vector, e.g. 
#' c("CG", "CHG") or c("CHG", "CG") for "CG/CHG" methylation.
#' @param bsgenome A \code{\link[BSgenome]{BSgenome}} object of the reference
#' genome.
#' @param ... Other arguments passed to \code{\link[Biostrings]{matchPattern}}.
#' 
#' @note This is basically a wrapper around the 
#' \code{\link[Biostrings]{matchPattern}} function, where 
#' \code{methylation_type} = \code{pattern} and \code{bsgenome} = 
#' \code{subject}.
#' \emph{WARNING:} This will use a considerable amount of memory and take some 
#' time for non-CG methylation.
#' 
#' @return A \link{MethylationLociSet} object.
#' 
#' @export
makeMLS <- function(methylation_type, bsgenome){
  
  # Argument checks
  if (!all(sapply(X = methylation_type, .valid_methylation_type))){
    stop("Invalid ", sQuote('methylation_type.'), 
         ".\nPlease see the help page for MethylationLociSet, which can ", 
         "accessed by ", "typing ", sQuote("?MethylationLociSet"), 
         " at the R prompt, for further details of this argument.")
  }
  if (!is(bsgenome, "BSgenome")){
    stop(sQuote("bsgenome"), " must be a ", sQuote("BSgenome"), " object.")
  }
  
  # Construct BSParams object
  bsparam <- new("BSParams", X = bsgenome, FUN = function(pdict, ...) {
    unlist(matchPDict(pdict = pdict, ...))
  }, 
  simplify = TRUE)
  
  # Convert the methylation_type to a PDict.
  # Then, search the bsgenome for that PDict.
  # Have to do in 2 steps because matchPDict only accepts patterns of the same 
  # width.
  
  # (1) CG methylation
  if ('CG' %in% methylation_type){
    # CG methylation is strand-symmetric.
    # Therefore only need to run search once.
    pdict_cg <- PDict('CG')
    cg_fwd <- bsapply(bsparam, pdict = pdict_cg)
    cg_rev <- cg_fwd
    gr_cg <- c(.irl2gr(cg_fwd, '+', seqinfo(bsgenome)), 
               .irl2gr(cg_rev, '-', seqinfo(bsgenome)))
  } else{
    gr_cg <- GRanges(seqinfo = seqinfo(bsgenome))
  }
  
  # (2) Non-CG methylation
  # Make PDict
  x <- DNAStringSet()
  if('CHG' %in% methylation_type){
    xx <- DNAStringSet(apply(
      expand.grid('C',
                  unlist(strsplit(IUPAC_CODE_MAP['G'], ''), use.names = FALSE),
                  unlist(strsplit(IUPAC_CODE_MAP['H'], ''), use.names = FALSE)),
      FUN = paste, MARGIN = 1, collapse = ''))
    x <- c(x, xx) 
  } 
  if ('CHH' %in% methylation_type){
    xx <- DNAStringSet(apply(
      expand.grid('C', 
                  unlist(strsplit(IUPAC_CODE_MAP['H'], ''), use.names = FALSE),
                  unlist(strsplit(IUPAC_CODE_MAP['H'], ''), use.names = FALSE)), 
      FUN = paste, MARGIN = 1, collapse = ''))
    x <- c(x, xx)
  }
  if ('CNN' %in% methylation_type){
    stop("Sorry, ", sQuote('CNN'), " methylation not currently supported.")
    # FIXME: Searching for CNN will match 16 combinations (CAA -> CTT).
    # Really want to match against C<masked><masked>, I think.
    xx <- DNAStringSet(apply(
      expand.grid('C', 
                  unlist(strsplit(IUPAC_CODE_MAP['N'], ''), use.names = FALSE),
                  unlist(strsplit(IUPAC_CODE_MAP['N'], ''), use.names = FALSE)),
      FUN = paste, MARGIN = 1, collapse = ''))
    x <- c(x, xx)
  }

  # Search for non-CG methylation (if required).
  if (length(x) != 0L){
    # Non-CG methylation is generally not strand-symmetric.
    # Therefore have to run search twice:
    # Firstly, with a the original dictional
    # Secondly,  with a reverse-complemented dictionary
    # Can't simply combine dictionaries and run together because there is no 
    # strand information on the returned output.
    pdict_non_cg_fwd <- PDict(x)
    pdict_non_cg_rev <- PDict(reverseComplement(x))
    non_cg_fwd <- bsapply(bsparam, pdict = pdict_non_cg_fwd)
    non_cg_rev <- bsapply(bsparam, pdict = pdict_non_cg_rev)
    
    gr_non_cg <- c(.irl2gr(non_cg_fwd, '+', seqinfo(bsgenome)), 
                   .irl2gr(non_cg_rev, '-', seqinfo(bsgenome)))
  } else{
    gr_non_cg <- GRanges(seqinfo = seqinfo(bsgenome))
  }
  
  # Construct MethylationLociSet
  new("MethylationLociSet", c(gr_cg, gr_non_cg), 
      methylation_type = methylation_type)
}