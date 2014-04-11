### =========================================================================
### CoMeth: methylation patterns at m-tuples
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("getPos", "CoMeth", function(x){
  getPos(rowData(x))
  })

setMethod("getMethylationType", "CoMeth", function(x){
  x@colData$methylation_type
})

setMethod("sampleNames", "CoMeth", function(x) {
  colnames(x)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

## TODO: Enforce seqinfo argument?
## TODO: Currently passsing "..." to both the MTuples() and SummarizedExperiment() constructors; this will likely lead to weird behaviour.
#CoMeth <- function(m, sample_names, pos, counts, strand, methylation_type, seqinfo, sort_cometh = TRUE){
CoMeth <- function(counts = DataFrameList(), seqnames = RleList(), pos = DataFrameList(), strand = RleList(), seqinfo = NULL, seqlengths = NULL, sample_names = CharacterList(), methylation_type = CharacterList(), colData = DataFrame(), exptData = SimpleList(), ..., verbose = FALSE){
  
  ## TODO: Argument checks, including checking for missing arguments

  if (!.zero_range(ncol(pos))){
    stop("Cannot combine ", sQuote("CoMeth"), " with different values of ", sQuote("m")) 
  }
  
  ## Create strand info if none supplied
  if (length(strand) == 0L){
    strand <- RleList(lapply(seqnames, function(x){Rle("*", length(x))}))
  }
  
  ## Add the correct levels to strand if they don't already exist
  ## Code based on GenomicRanges:::newGRanges
  strand <- endoapply(X = strand, FUN = function(x){
    if (!is.factor(runValue(x)) || !identical(levels(runValue(x)), levels(strand()))){ 
      runValue(x) <- strand(runValue(x))
    }
    return(x)
  })
  
  ## Join methylation_type elements into a single string.
  methylation_type <- paste0(sort(unique(unlist(strsplit(unlist(methylation_type), split = '/')))), collapse = '/')

  ## Combine the list-wise data into matrix-like data whilst taking care of common and sample-specific m-tuples, i.e. filling in zeros when a sample doesn't have any observations for that m-tuple.
  combined_data <- .combine(sample_names = sample_names, seqnames = seqnames, pos = pos, strand = strand, counts = counts)
    
  ## Construct rowData of CoMeth object, which is a MTuples object.
  mtuples <- MTuples(seqnames = combined_data$seqnames, pos = combined_data$pos, strand = combined_data$strand, ..., seqlengths = seqlengths, seqinfo = seqinfo)
  
  ## Construct colData of CoMeth object.
  ## Include any coLData supplied as an argument to the CoMeth() constructor.
  if (identical(dim(colData), c(0L, 0L))){
    colData <- DataFrame(methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(sample_names)), row.names = unlist(sample_names))
  } else{
    colData <- cbind(DataFrame(m = rep(ncol(pos), length(sample_names)), methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(sample_names)), row.names = unlist(sample_names)), colData)
  }
  
  ## Construct CoMeth object
  new("CoMeth", SummarizedExperiment(assays = combined_data$counts, rowData = mtuples, colData = colData, exptData = exptData, ..., verbose = verbose))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Filtering
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read .tsv file(s) from comethylation and construct a CoMeth object
###

read.comethylation <- function(files, methylation_types, sample_names, m, seqlengths = NULL, seqinfo = NULL, quiet = FALSE, ) {
  
  ## TODO: Argument checks
  ## TODO: Is m necessary?
  ## TODO: Check that seqinfo is supplied if CoMeth makes it mandatory.
  
  x <- lapply(files, function(file){
    if (grepl("\\.gz$", file)){
      con <- gzfile(file)
    } else if (grepl("\\.bz2$", file)){
      con <- bzfile(file)
    } else {
      con <- file(file)
    }
    
    # Construct header and column types. The first column is a character, the rest are int
    column_headers <- strsplit(readLines(con, n = 1), '\t')[[1]]
    what0 <- replicate(length(column_headers), character(0))
    names(what0) <- column_headers
    what0[-1] <- replicate(length(column_headers) - 1, integer(0))
  
    tsv <- scan(con, skip = 1, what = what0, sep = "\t", quote = "", na.strings = "NA", quiet = quiet, comment.char = "")
    close(con)
    val <- list(counts = DataFrame(tsv[grep('[MU]', names(tsv))]), seqnames = Rle(tsv[['chr']]), pos = DataFrame(tsv[grep('^pos', names(tsv))]))
    
    return(val)
  })
  
  names(x) <- sample_names
  methylation_type <- as(methylation_types, "CharacterList")
  names(methylation_type) <- sample_names
  
  CoMeth(counts = DataFrameList(lapply(X = x, function(y){y$counts})), seqnames = RleList(lapply(X = x, function(y){y$seqnames})), pos = DataFrameList(lapply(X = x, function(y){y$pos})), seqinfo = seqinfo, sample_names = sample_names, methylation_type = methylation_type)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

## The show method is adapted from that of SummarizedExperiment
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