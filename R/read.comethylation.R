## TODO: Remove dependency on R.utils
## TODO: Add bzip support.
## TODO: Profile.
## TODO: Add offset (which is a parameter of .LOR)?


#' Read \code{tsv} output files from \code{comethylation} software.
#'
#' @description
#' Read the \code{tsv} output files from the Python software 
#' \code{comethylation} and construct a single \code{\link{CoMeth}} object. All 
#' files should be of the same size m-tuples, e.g. all 2-tuples. Files may have 
#' different \code{methylation_type}, although this is not recommended. All 
#' files should be mapped against the same reference genome, which is supplied 
#' as the \code{seqinfo} argument.
#'
#' @param files The \code{.tsv} files created by \code{comethylation}. These 
#' files may be compressed with gzip or bzip2.
#' @param sample_names The sample names of each file.
#' @param methylation_types A character vector with the type of methylation 
#' type of each file.
#' @param verbose logical Default \code{FALSE}. If \code{TRUE},
#' \code{read.comethylation()} will report some output when reading the files 
#' with \code{\link[data.table]{fread}} and may produce some output when 
#' constructing the \code{CoMeth} object.
#' @param seqinfo A \code{\link[GenomicRanges]{Seqinfo}} object containing 
#' information about the reference genome of the sample.
#' 
#' @note The compression of \code{file} is determined by the file extension: 
#' \code{.gz} for compressed with gzip and \code{.bz2} for compressed with 
#' bzip2. Otherwise the file is assumed to be uncompressed.
#'
#' @seealso \code{\link{CoMeth}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
#' 
#' @export
read.comethylation <- function(files, sample_names, methylation_types, seqinfo,
                                verbose = FALSE) {
  
  # Check there is a unique sample_name for each filename
  if (length(sample_names) != length(files) | 
        length(sample_names) != length(unique(sample_names))) {
    stop("Each file must have a unique sample name.")
  }
  
  # Construct colData
  if (length(methylation_types) != length(sample_names)) {
    stop("Much specify methylation type of each file")
  }
  methylation_type <- paste0(
    sort(unique(unlist(strsplit(methylation_types, '/')))), collapse = '/')
  if(!.valid_methylation_type(methylation_type)) {
    stop("Invalid methylation type.")
  }
  colData <- DataFrame(methylation_type = 
                         rep(methylation_type, length(sample_names)),
                       row.names = sample_names)
  
  # Check seqinfo is supplied
  if (missing(seqinfo)){
    stop("Must supply ", sQuote('seqinfo'), '.')
  }
  
  # Read in files and return as MTuples plus matrix of counts
  x <- bplapply(files, function(file, verbose) {
    if (grepl("\\.gz$", file)){
      file <- gunzip(file, temporary = TRUE, remove = FALSE)
    } else if (grepl("\\.bz2$", file)){
      stop("Sorry, can't handle ", sQuote('bzip2'), " compressed files.")
    } else {
      # Nothing to do!
    }
    
    tsv <- fread(file, sep = '\t', header = TRUE, verbose = verbose)
    val <- list(mtuples = MTuples(seqnames = tsv[, chr], strand = tsv[, strand], 
                                  pos = as.matrix(tsv[, grep('^pos',names(tsv)), 
                                                      with = FALSE])),
                counts = c(tsv[, grep('[MU]', names(tsv)), 
                               with = FALSE])
    )
    return(val)
  }, verbose = verbose)
  
  # Compare all sample m-tuples (smt) against the sample with the most 
  # m-tuples (i). Create a list of any putative sample-specific m-tuples and 
  # unique-ify it (ssmt). Then, concatenate m-tuples from sample i with ssmt 
  # (mt). By construction, mt will not contain any duplicate m-tuples.
  smt <- lapply(x, '[[', 'mtuples')
  m <- sapply(smt, getM)
  if (!.zero_range(m)) {
    stop("All files must contain the same sized m-tuples.")
  }
  m <- m[1L]
  i <- which.max(sapply(smt, length))
  # Compare all other samples against i
  u <- bplapply(seq_along(smt)[-i], function(ii, x, y){
    Rle(!overlapsAny(x[[ii]], y, type = 'equal'))
  }, x = smt, y = smt[[i]])
  ssmt <- unique(do.call("c", mapply('[', smt[-i], u)))
  mt <- sort(c(smt[[i]], ssmt))
  seqinfo(mt) <- seqinfo
  
  # Find overlaps between smt and mt for each sample.
  ol <- bplapply(smt, function(smt, mt) {
    findOverlaps(smt, mt, type = 'equal')
  }, mt = mt)
  
  # Construct assays
  counts <- lapply(x, '[[', 'counts')
  an <- names(counts[[1]])
  nmt <- length(mt)
  x <- bplapply(an, function(an, counts, ol, nmt) {
    mat <- matrix(NA_integer_, nrow = nmt, ncol = length(counts))
    ri <- unlist(lapply(ol, subjectHits), use.names = FALSE)
    ci <- rep(seq_along(counts), times = sapply(ol, queryLength))
    mat[ri + (ci - 1) * nrow(mat)] <- unlist(lapply(counts, '[[', an), 
                                             use.names = FALSE)
    return(mat)
  }, counts = counts, ol = ol, nmt = nmt)
  names(x) <- an
  if (m == 1L) {
    class <- 'CoMeth1'
    assays <- c(x, list(EP = .EP(x), beta = .beta(x)))
  } else if (m == 2L) {
    assays <- c(combined_data$counts, list(EP = .EP(combined_data$counts)), LOR = list(.LOR(combined_data$counts)))
    class <- 'CoMeth2'
  } else{
    assays <- c(combined_data$counts, list(EP = .EP(combined_data$counts)))
    class <- 'CoMethPlus'
  }
  
  new(class, SummarizedExperiment(assays = assays, rowData = mt, 
                                  colData = colData, verbose = verbose))
}