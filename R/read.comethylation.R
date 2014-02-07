# FIXME: This is slowed by list-to-matrix conversions. Perhaps it would be better to directly pass the columsn of the 'counts' and 'pos' matrices as lists to CoMeth. This would require modifications to CoMeth.

#' Parsing \code{tsv} output from \code{comethylation.py}
#' 
#' Read the \code{tsv} output file from \code{comethylation.py} and construct a \code{\link{CoMeth}} object. 
#' 
#' @param file The \code{.tsv} file created by comethylation.py. This file may be compressed with gzip or bzip2
#' @param m The size of the m-tuples
#' @param methylation_type A character vector with the type of methylation
#' @param sample_name The sample name
#' @param quiet logical: if \code{FALSE} (default), read.comethylation() will print a line, saying how many items (m-tuples) have been read.
#' @param seqlengths The sequence lengths of the reference genome of the sample. Must be an integer vector named with the sequence names and containing the lengths (or NA) for each level(seqnames).
#' @param seqinfo An (optional) \code{\link[GenomicRanges]{Seqinfo}} object containing information about the reference genome of the sample
#' @note The compression of \code{file} is determined by the file extension: \code{.gz} for compressed with gzip and \code{.bz2} for compressed with bzip2. Otherwise the file is assumed to be uncompressed.
#' 
#' @export
#' @seealso \code{\link{CoMeth}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
read.comethylation <- function(file, m, methylation_type, sample_name, quiet = FALSE, seqlengths = NULL, seqinfo = NULL) {
  
  # Construct header and column types. The first column is character, the rest are int
  if (grepl("\\.gz$", file)){
    con <- gzfile(file)
  } else if (grepl("\\.bz2$", file)){
    con <- bzfile(file)
  } else {
    con <- file(file)
  }
  column_headers <- strsplit(readLines(con, n = 1), '\t')[[1]]
  close(con)
  what0 <- replicate(length(column_headers), character(0))
  names(what0) <- column_headers
  what0[-1] <- replicate(length(column_headers) - 1, integer(0))
  
  # Read the file
  if (grepl("\\.gz$", file)){
    con <- gzfile(file)
  } else if (grepl("\\.bz2$", file)){
    con <- bzfile(file)
  } else {
    con <- file(file)
  }
  tsv <- scan(con, skip = 1, what = what0, sep = "\t", quote = "", na.strings = "NA", quiet = quiet, comment.char = "")
  close(con)
  
  # Pass to CoMeth constructor
  seqnames <- tsv[['chr']]
  pos <- matrix(unlist(tsv[grep('pos', names(tsv))]), ncol = m) # FIXME: Slow
  colnames(pos) <- names(tsv[grep('pos', names(tsv))])
  counts <- matrix(unlist(tsv[grep('[MU]', names(tsv))]), ncol = 2 ^ m) # FIXME: Slow
  colnames(counts) <- names(tsv[grep('[MU]', names(tsv))])
  cometh <- CoMeth(seqnames = seqnames, pos = pos, counts = counts, m = as.integer(m), methylation_type = methylation_type, sample_name = sample_name, strand = "*", seqlengths = seqlengths, seqinfo = seqinfo)

  # Return cometh
  return(cometh)
}