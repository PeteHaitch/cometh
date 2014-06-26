## TODO: Use a different delimiter to '.' when adding sample names to count 
## names. Otherwise a sample name with '.' in it will break the code. Whatever 
## delimiter I end up using, add a check that this delimiter is not present in 
## any of the sample_names.
## TODO: Remove dependency on R.utils::gunzip and improve gunzipping and temp
## file handling.
## TODO: Add bzip support.
## TODO: Profile.
## TODO: Add offset as a parameter (which is a parameter of .LOR)?
## TODO: Add timing output?
## TODO: Add min_cov option; only those m-tuples with sequencing coverage > 
## min_cov will be retained. This will greatly reduce the object sizes, 
## particularly for large m, and generally we are not interested in m-tuples 
## with low coverage in any case. Perhaps the minimum coverage should be 
## recorded in a slot of the CoMeth object.
## TODO: See if assays can be stored as DataFrames of Rles and whether this is 
## more efficient than matrices of integers.
## TODO: Compute approximate memory usage.
## TODO: Check that all files only contain either strand %in% '*' or strand 
## %in% c('+', '-'); otherwise the resulting object is a bit of a mess.

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
#' files may be compressed with gzip or bzip2 (not yet implemented). All files 
#' must contain the same sized m-tuples.
#' @param sample_names The sample names of each file. Must be unique.
#' @param methylation_types A character vector with the type of methylation 
#' type of each file. It is recommended, although not enforced, that all samples 
#' have the same methylation type. A sample with multiple methylation types 
#' should have these separated by a forward-slash, for example, 'CG/CHG'. 
#' @param seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing 
#' information about the reference genome of the sample.
#' @param colData A \code{\link[GenomicRanegs]{DataFrame}} describing the samples. 
#' Row names, which are required, become the column names of the \code{CoMeth}. 
#' Must not include a column named \code{methylation_type} andrRow names must 
#' be identical to \code{sample_names}.
#' @param exptData An optional \code{\link[IRanges]{SimpleList}} of arbitrary 
#' content describing the overall experiment.
#' @param verbose logical Default \code{FALSE}. If \code{TRUE},
#' \code{read.comethylation()} will report some output when reading the files 
#' with \code{\link[data.table]{fread}} and may produce some output when 
#' constructing the \code{CoMeth} object.
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
                               colData = DataFrame(), exptData = SimpleList(), 
                               verbose = FALSE) {
  
  # Check there is a unique sample name for each filename
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
  if (identical(dim(colData), c(0L, 0L))) {
    colData <- DataFrame(methylation_type = 
                           rep(methylation_type, length(sample_names)), 
                         row.names = sample_names)
  } else{
    if (any(colnames(colData) == 'methylation_type')) {
      stop(paste0(sQuote('colData'), " must not include a column called ", 
                  sQuote('methylation_type'), "."))
    }
    if (!identical(row.names(colData), sample_names)) {
      stop(sQuote('row.names'), " of ", sQuote('colData'), " must match ",
           sQuote('sample_names'), ".")
    }
    cbind(colData, DataFrame(methylation_type = 
                               rep(methylation_type, length(sample_names)), 
                             row.names = sample_names))
  }
  
  # Check seqinfo is supplied
  if (missing(seqinfo)){
    stop("Must supply ", sQuote('seqinfo'), '.')
  }
  
  # Decompress files (if required).
  if (any(grepl("\\.gz$", files) || any(grepl("\\.bz2$", files)))) {
    print("Decompressing files...")
    x <- do.call("rbind", bplapply(files, function(file, verbose) {
      if (grepl("\\.gz$", file)){
        file <- gunzip(file, remove = FALSE)
        file_to_remove <- file
      } else if (grepl("\\.bz2$", file)){
        stop("Sorry, can't handle ", sQuote('bzip2'), " compressed files.")
        # TODO: Uncomment once I have a way to deal with bzip files
        # file <- bzip(file)
        # file_to_remove <- file
      } else {
        # Nothing to do!
        file_to_remove <- NA
      }
      return(data.frame(file = file, file_to_remove = file_to_remove, 
                        stringsAsFactors = FALSE))
    }))
  } else{
    # No files require decompressing
    x <- data.frame(file = files, file_to_remove = rep(NA, length(files)), 
                    stringsAsFactors = FALSE)
  } 
  
  # Read in file(s) serially and, if more than one file, merge these files.
  if (length(files) == 1L) {
    print(paste0("Reading file ", files[1]))
  } else {
    print(paste0("Reading file ", files[1]))
  }
  mtsv <- fread(input = x$file[1], header = TRUE, verbose = verbose)
  keys <- c("chr", "strand", 
            grep(pattern = '^pos', x = names(mtsv), value = TRUE))
  setkeyv(x = mtsv, cols = keys, verbose = verbose)
  # Append the sample name to the names of the counts, e.g. 'M' -> 'M.sample1'
  wn <- which(!(names(mtsv) %in% keys))
  setnames(mtsv, c(names(mtsv)[seq_len(ncol(mtsv))[-wn]], 
                   paste(names(mtsv)[wn], sample_names[1], sep = ".")))
  rm(wn)
  if (length(files) > 1) {
    for (i in seq.int(from = 2, to = length(files), by = 1)) {
      print(paste0("Reading and merging ", files[i]))
      ## TODO: try-catch the merge in case files with different sized m-tuples
      ## are supplied by the user because the error message supplied by 
      ## merge is a bit hard to understand.
      tmp <- fread(input = x$file[i], header = TRUE, verbose = verbose)
      setkeyv(x = tmp, cols = keys, verbose = verbose)
      wn <- which(!(names(tmp) %in% keys))
      setnames(tmp, c(names(tmp)[seq_len(ncol(tmp))[-wn]], 
                      paste(names(tmp)[wn], sample_names[i], sep = ".")))
      mtsv <- merge(mtsv, tmp, by = keys, all = TRUE)
      rm(tmp, wn)
    }
  }
  # Remove the decompressed files (if required)
  unlink(x$file_to_remove)
  
  # Combine the counts into the assays list
  print("Creating assays...")
  assay_names <- unique(sapply(strsplit(grep(pattern = '\\.', x = names(mtsv), 
                                             value = TRUE),
                                        split = '\\.'), '[[', 1))
  assays <- lapply(assay_names, function(an, mtsv) {
    pat <- paste0('^', an, '.')
    i <- grep(pattern = pat, x = names(mtsv), value = TRUE)
    sn <- sapply(strsplit(x = i, split = '\\.'), '[[', 2)
    as.matrix(setnames(mtsv[, i, with = FALSE], sn))
  }, mtsv = mtsv)
  names(assays) <- assay_names
  
  rowData <- MTuples(seqnames = mtsv[, chr], strand = mtsv[, strand], 
                     pos = as.matrix(mtsv[, grep(pattern = '^pos', 
                                                 x = colnames(mtsv)), 
                                          with = FALSE]), 
                     seqinfo = seqinfo)
  rm(mtsv)
  
  m <- getM(rowData)
  if (m == 1L) {
    class <- 'CoMeth1'
    assays <- c(assays, list(beta = .beta(assays)))
  } else if (m == 2L) {
    class <- 'CoMeth2'
    assays <- c(assays, list(LOR = .LOR(assays)))
  } else{
    class <- 'CoMeth3Plus'
    # No need to update assays
  }
  
  # Create the CoMeth object
  print("Creating CoMeth object...")
  new(class, SummarizedExperiment(assays = assays, rowData = rowData, 
                                  colData = colData, exptData = exptData, 
                                  verbose = verbose))
}