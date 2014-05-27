## TODO: If this function is exported then need to add data.table and R.utils
## to Imports in DESCRIPTION.
#' A faster version of \code{\link{read.comethylation}}
#' 
#' Uses \code{data.table::fread} for faster reading of \code{comethylation} 
#' \code{.tsv} output files. See https://github.com/PeteHaitch/cometh/issues/18 
#' for benchmarking and a discussion of pros and cons of this approach.
#' @keywords internal
fread.comethylation <- function(files, sample_names, methylation_types, seqinfo,
                                verbose = FALSE) {
  
  # Check seqinfo is supplied
  if (missing(seqinfo)){
    stop("Must supply ", sQuote('seqinfo'))
  }
  
  # Read in each file and store in a list
  x <- lapply(files, function(file, verbose){
    if (grepl("\\.gz$", file)){
      file <- gunzip(file, temporary = TRUE)
    } else if (grepl("\\.bz2$", file)){
      stop("Sorry, can't handle ", sQuote('bzip2'), " compressed files.")
    } else {
      # Nothing to do!
    }
    
    tsv <- fread(file, sep = '\t', header = TRUE, verbose = verbose)
    val <- list(counts = DataFrame(tsv[, grep('[MU]', names(tsv)), with = FALSE]), 
                seqnames = Rle(tsv[, chr]), 
                pos = DataFrame(tsv[, grep('^pos', names(tsv)), with = FALSE]))
    
    return(val)
  }, verbose = verbose)
  
  # Convert list of files into objects that can be passed to the CoMeth constructor
  names(x) <- sample_names
  sample_names <- as(sample_names, "CharacterList")
  methylation_type <- as(methylation_types, "CharacterList")
  names(methylation_type) <- sample_names
  
  CoMeth(sample_names = sample_names, methylation_type = methylation_type, counts = DataFrameList(lapply(X = x, function(y){y$counts})), seqnames = RleList(lapply(X = x, function(y){y$seqnames})), pos = DataFrameList(lapply(X = x, function(y){y$pos})), seqinfo = seqinfo)
}