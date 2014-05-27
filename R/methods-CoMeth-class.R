### =========================================================================
### CoMeth: methylation patterns at m-tuples
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

#' @include AllGenerics.R
#' @export
setMethod("getPos", "CoMeth", function(x){
  getPos(rowData(x))
  })

#' @include AllGenerics.R
#' @export
setMethod("getIPD", "CoMeth", function(x){
  getIPD(rowData(x))
})

#' @include AllGenerics.R
#' @export
setMethod("getMethylationType", "CoMeth", function(x){
  x@colData$methylation_type
})

#' @include AllGenerics.R
#' @export
setMethod("sampleNames", "CoMeth", function(object) {
  colnames(object)
})

#' @include AllGenerics.R
#' @export
setReplaceMethod("sampleNames", signature = signature(object = "CoMeth", value = "ANY"), function(object, value) {
  if (length(value) != length(sampleNames(object))){
    stop("Invalid ", sQuote('sampleNames'), " length")
  }
  colnames(object) <- value
  object
})

#' @include AllGenerics.R
#' @export
setMethod("getM", "CoMeth", function(x) {
  getM(rowData(x))
})

#' Get the \code{coverage} from a \code{\link{CoMeth}} object.
#'
#' The \code{coverage} of a \code{\link{CoMeth}} object is the number of reads for each m-tuple, i.e. \code{rowSums(getCounts(CoMeth))}. The current implementation is slow.
#' @param x A \code{\link{CoMeth}} object
#'
#' @return A numeric matrix. Each column of the matrix corresponds to a sample and each row to an m-tuple.
#' @include AllGenerics.R
#' @export
setMethod("getCoverage", "CoMeth", function(x) {
  m <- getM(x)
  assay_names <- .make_m_tuple_names(m)
  Reduce("+", assays(x)[assay_names])
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

## TODO: As of GenomicRanges_1.16.0, assays can be DataFrame objects and not just matrices. Decide which is a better type for the CoMeth class.

#' The constructor function for CoMeth objects
#'
#' This is the general function to construct a \code{CoMeth} object. Users will generally use the \code{\link{read.comethylation}} function to construct a \code{CoMeth} object. The constructor allows for multiple samples in the same object by passing most arguments as a *List where each element of the *List corresponds to the arguments for a given sample. These multiple samples will then be appropriately combined into a single \code{CoMeth} object.
#' 
#' @param sample_names A \code{\link[IRanges]{CharacterList}} containing the 
#' names of the samples. Sample names must be unique.
#' @param methylation_type A \code{\link[IRanges]{CharacterList}} storing the 
#' type of methylation loci for these m-tuples. \code{methylation_type} must be 
#' named and the names must match those given in \code{sample_names}. For each 
#' sample the possible values are "CG", "CHG", "CHH" or "CNN" or multiple values 
#' specified as a character vector, e.g. c("CG", "CHG") or c("CHG", "CG") for 
#' "CG/CHG" methylation. The \code{methylation_type} of the resulting 
#' \code{CoMeth} object is the union of this argument, e.g. if 
#' \code{methylation_type = list('sample1' = 'CG', 'sample2' = 'CHH')} then the 
#' \code{methylation_type} of the resulting \code{CoMeth} object is 
#' \code{CG/CHH}; this also gives a warning.
#' @param counts A \code{\link[IRanges]{DataFrameList}}. \code{counts} must be 
#' named and the names must match those given in \code{sample_names}. The 
#' entries in each DataFrame corresponds to the number of times that particular 
#' co-methylation pattern (columns) was observed for that particular m-tuple 
#' (rows). Therefore, each DataFrame must have the same number of rows as its 
#' corresponding DataFrame in \code{pos} and have \eqn{2 ^ \code{m}} columns. 
#' All samples must have the same m, that is, all DataFrames 
#' in \code{counts} must have the same number of columns.
#' The column names of each DataFrame must match those given by 
#' \code{.make_m_tuple_names(m)}. 
#' @param seqnames a \code{\link[IRanges]{RleList}}. \code{seqnames} must be 
#' named and the names must match those given in \code{sample_names}. Each Rle 
#' is the \code{\link[GenomicRanges]{seqnames}} of each m-tuple for that sample.
#' @param pos A \code{\link[IRanges]{DataFrameList}}. \code{pos} must be named 
#' and the names must match those given in \code{sample_names}. The columns of 
#' each DataFrame must be: \code{seqnames}, \code{pos1}, ..., \code{posm}, 
#' where, for example, \code{posm} is \code{pos3} if \code{m} = 3. 
#' Therefore, the number of columns of each DataFrame is m, the size of the 
#' m-tuples, and the number of rows is equal to the number of m-tuples for that 
#' particular sample. All samples must have the same m, that is, all DataFrames 
#' in \code{pos} must have the same number of columns.
#' @param seqinfo A \code{\link[GenomicRanges]{Seqinfo}} object containing 
#' information about the reference genome of the samples. Only a single value is 
#' accepted because all samples must be mapped against the same reference 
#' genome. However, multiple genomes per sample are allowed in order to 
#' accommodate a spiked-in unmethylated genome (normally lambda phage), which is 
#' a common step in a bisulfite-sequencing protocol. \strong{NOTE}: 
#' \code{seqinfo} is a required argument in order to construct a \code{CoMeth} 
#' object, whereas for many Bioconductor objects it is an optional argument.
#' @param strand An optional \code{\link[IRanges]{RleList}} object containing 
#' the strand information of each m-tuple. \strong{WARNING}: If \code{strand} is 
#' not supplied, all m-tuples in the resulting \code{CoMeth} object will 
#' have their strand set to \code{*} to signify that the strand is unknown or 
#' irrelevant. 
#' \strong{WARNING}: m-tuples will not be combined across samples if they are on 
#' different strands. 
#' @param colData An optional \code{\link[IRanges]{DataFrame}} describing the 
#' samples. Row names must match the \code{sample_names}, otherwise an error is
#' returned.
#' \code{sample_names}.
#' @param exptData An optional \code{\link[IRanges]{SimpleList}} of arbitrary 
#' content describing the overall experiment. 
#' @param ... Additional arguments passed to the internal call to the 
#' \code{\link{MTuples}} constructor. See \code{\link{MTuples}} for details on 
#' what these might be.
#' @param verbose A \code{logical(1)} indicating whether messages about data 
#' coercion during construction should be printed.
#' 
#' @export
#' @seealso \code{\link{read.comethylation}} for a function to read in the 
#' \code{.tsv} output file of \code{comethylation} and construct a \code{CoMeth}
#' object. 
#' @seealso \code{\link[GenomicRanges]{SummarizedExperiment}} for the class that
#' \code{CoMeth} extends.
#' 
#' @return A \code{CoMeth1} (if \code{m} \eqn{= 1}), \code{CoMeth2} (if 
#' \code{m} \eqn{= 2}) or \code{CoMeth3Plus} (if \code{m} \eqn{>= 3}) object. 
#' All these are concrete subclasses of the VIRTUAL \code{\link{CoMeth}} class.
#' 
#' 
#' cat("TODO")
CoMeth <- function(sample_names = CharacterList(), methylation_type = CharacterList(), counts = DataFrameList(), seqnames = RleList(), pos = DataFrameList(), seqinfo = Seqinfo(), strand = RleList(), colData = DataFrame(), exptData = SimpleList(), ..., verbose = FALSE){
  
  ## Check that all required arguments are not missing
  if (missing(sample_names)){
    stop(sQuote("sample_names"), " missing.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (missing(methylation_type)){
    stop(sQuote('methylation_type'), " missing. Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (missing(counts)){
    stop(sQuote('counts'), " missing. Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (missing(seqnames)){
    stop(sQuote('seqnames'), " missing. Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (missing(pos)){
    stop(sQuote('pos'), " missing. Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (missing(seqinfo)){
    stop(sQuote('seqinfo'), " missing. Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }

  ## Check that all required arguments are of the correct class 
  if (!is(sample_names, 'CharacterList')){
    stop(sQuote('sample_names'), " must be a ", sQuote('CharacterList'), ". Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!is(methylation_type, 'CharacterList')){
    stop(sQuote('methylation_type'), " must be a ", sQuote('CharacterList'), ". Please see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!is(counts, "DataFrameList")){
    stop(sQuote('counts'), " must be a ", sQuote('DataFrameList'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!is(seqnames, "RleList")){
    stop(sQuote('seqnames'), " must be an ", sQuote("RleList"), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!is(pos, "DataFrameList")){
    stop(sQuote('pos'), " must be a ", sQuote("DataFrameList"), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!is(seqinfo, 'Seqinfo')){
    stop(sQuote('seqinfo'), " must be a ", sQuote('Seqinfo'), " object.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  
  ## Check that 'sample_names' are unique
  if (!identical(unique(unlist(sample_names)), unlist(sample_names))){
    stop("Each element of ", sQuote("sample_names"), " must be unique.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  
  ## Check that all required arguments that are *Lists (e.g. RleList, DataFrameList, etc.) have the same names as the 'sample_names'
  if (!setequal(names(methylation_type), unlist(sample_names))){
    stop("Names of ", sQuote('methylation_type'), " must match those in ", sQuote('sample_names'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!setequal(names(counts), unlist(sample_names))){
    stop("Names of ", sQuote('counts'), " must match those in ", sQuote('sample_names'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!setequal(names(seqnames), unlist(sample_names))){
    stop("Names of ", sQuote('seqnames'), " must match those in ", sQuote('sample_names'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!setequal(names(pos), unlist(sample_names))){
    stop("Names of ", sQuote('pos'), " must match those in ", sQuote('sample_names'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }  
  
  ## Check that all required arguments have the correct dimensions
  if (!.zero_range(ncol(counts)) || !identical(log2(ncol(counts[[1]])), round(log2(ncol(counts[[1]])), 0))){
    stop(sQuote('ncol(counts)'), " must be identical for all elements of ", sQuote('counts'), " and should be a power of 2.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of these arguments.")
  }
  m <- as.integer(log2(ncol(counts[[1]])))
  if (!identical(nrow(counts[unlist(sample_names)]), nrow(pos[unlist(sample_names)]))){
    stop(sQuote('nrow(counts)'), " must be identical to ", sQuote('nrow(pos)'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of these arguments.")
  }
  if (!identical(sapply(seqnames[unlist(sample_names)], length), nrow(pos))){
    stop("Length of each ", sQuote('Rle'), " element of ", sQuote("seqnames"), " must be identical to ", sQuote('nrow(pos)'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of these arguments.")
  }
  if (!isTRUE(all(ncol(pos) == m))){
    stop(sQuote('m'), " is inferred to be ", m, " from ", sQuote("counts"), " data but not all ", sQuote("DataFrame"), " elements of ", sQuote("pos"), " have the same inferred value of ", sQuote("m"), ", namely ", sQuote('ncol(pos)'), " = (", paste0(ncol(pos), collapse = ", "), ")\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of these arguments.")
  }

  ## Check that all elements of *List required arguments are the correct type, correct dimensions and, if required, have the correct names.
  ## NOTE: These checks can't be made until 'm' is computed during the "check correctness of dimensions and names of required arguments" of the *Lists themselves.
  if (!all(unlist(methylation_type) %in% c('CG', 'CHG', 'CHH', 'CNN'))){
    stop(sQuote('methylation_type'), " for each sample must be ", sQuote('CG'), ", ", sQuote('CHG'), ", ", sQuote('CHH'), " or ", sQuote('CNN'), ", or a vector of some combination of these, e.g., ", sQuote("c('CG', 'CHG')"), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!all(sapply(counts, function(x, m){isTRUE(all.equal(sapply(x, class), rep('integer', 2 ^ m), check.attributes = F))}, m = m))){
    stop("Class of columns in each ", sQuote('DataFrame'), " element of ", sQuote("counts"), " must be: ", paste0(rep(sQuote('integer'), 2 ^ m), collapse = ", "), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (any(sapply(counts, function(x, m){!identical(colnames(x), .make_m_tuple_names(m))}, m = m))){
    stop(sQuote('m'), " is set to ", m, " so colnames for all elements of ", sQuote('counts'), " must be: ", paste0(sQuote(.make_m_tuple_names(m)), collapse = ", "), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (!all(sapply(pos, function(x, m){isTRUE(all.equal(sapply(x, class), rep('integer', m), check.attributes = F))}, m = m))){
    stop("Class of columns in each ", sQuote("DataFrame"), " element of ", sQuote('pos'), "  must be: '", paste0(rep('integer', m), collapse = "', '"), "'.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  if (any(sapply(pos, function(x, m){!identical(colnames(x), paste0('pos', seq_len(m)))}, m = m))){
    stop(sQuote('m'), " is set to ", m, " so colnames for all elements of ", sQuote('pos'), " must be: ", paste0(sQuote(paste0('pos', seq_len(m))), collapse = ', '), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  
  ## If supplied, check the (non-required) 'strand' argument.
  ## Otherwise, construct the default, which encodes "no strand information".
  if (!is(strand, "RleList")){
    stop(sQuote('strand'), " must be an ", sQuote('RleList'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  ### Create strand info if none supplied
  if (length(strand) == 0L){
    strand <- RleList(lapply(seqnames, function(x){Rle("*", length(x))}))
  }
  ### Add the correct levels to strand if they don't already exist
  ### Code copied from GenomicRanges:::newGRanges
  strand <- endoapply(X = strand, FUN = function(x){
    if (!is.factor(runValue(x)) || !identical(levels(runValue(x)), levels(strand()))){ 
      runValue(x) <- strand(runValue(x))
    }
    return(x)
  })
  if (!identical(sapply(strand, length), sapply(pos, nrow))){
    stop("Length of each ", sQuote('Rle'), " element in ", sQuote('strand'), " must equal the number of rows of its corresponding element in ", sQuote('pos'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  
  ## If supplied, check the (non-required) 'colData' argument.
  if (!is(colData, "DataFrame")){
    stop(sQuote('colData'), ", if supplied, must be a ", sQuote('DataFrame'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }

  ## If supplied, check the (non-required) 'exptData' argument.
  if (!is(exptData, "SimpleList")){
    stop(sQuote('exptData'), ", if supplied, must be a ", sQuote('SimpleList'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
  }
  
  ## Warnings about particular parameter combinations
  if (isTRUE(any(unlist(methylation_type) != methylation_type[[1]]))){
    warning("The supplied ", sQuote('methylation_type'), " parameter says that samples contain data for different methylation types. The union of these methylation types will be used as the ", sQuote('methylation_type'), " of the returned ", sQuote('CoMeth'), " object.")
  }
  if (c("*", "+", "-") %in% as.vector(unlist(strand, use.names = FALSE)) %*% c(1, 0.5, 0.5) > 1){
    warning("The supplied ", sQuote('strand'), " argument contains ", sQuote('*'), " as well as at least one of ", sQuote('+'), " or ", sQuote('-'), ". m-tuples will not be combined across samples if they are on different strands. This means that m-tuples with the exact same genomic co-ordinates across samples, but with different strands (e.g. ", sQuote('*'), " and ", sQuote('+'), "), will not be combined.\nTo combine all m-tuples based on genomic co-ordinates regardless of strand, do not set the ", sQuote('strand'), " parameter when calling the ", sQuote('CoMeth'), " constructor.")
  }

  ## Combine data from multiple samples (also works if only a single sample is provided).
  ### Join elements of 'methylation_type' into a single string.
  methylation_type <- paste0(sort(unique(unlist(strsplit(unlist(methylation_type), split = '/')))), collapse = '/')
  ### Combine the *List-wise data into matrix-like data, whilst taking care of common and sample-specific m-tuples, i.e. filling in NAs when a sample doesn't have any observations for that m-tuple.
  combined_data <- .combine(sample_names = sample_names, seqnames = seqnames, pos = pos, strand = strand, counts = counts)
  
  ## Construct rowData of CoMeth object, which is a MTuples object.
  mtuples <- MTuples(seqnames = combined_data$seqnames, pos = combined_data$pos, strand = combined_data$strand, seqlengths = NULL, seqinfo = seqinfo, ...)
  
  ## If colData not supplied, then create the default colData, which simply stores the 'sample_names' and 'methylation_type'
  ## If coLData supplied, then combine the coLData with the default.
  if (identical(dim(colData), c(0L, 0L))){
    colData <- DataFrame(methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(sample_names)), row.names = unlist(sample_names))
  } else{
    if (!setequal(row.names(colData), unlist(sample_names))){
      stop(sQuote('row.names(colData)'), " must be identical to ", sQuote('sample_names'), ".\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of this argument.")
    } else{
      colData <- DataFrame(methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(sample_names)), colData, row.names = unlist(sample_names))
    }
  }
  
  ## Construct CoMeth object
  if (m == 1L){
    assays <- c(combined_data$counts, list(EP = .EP(combined_data$counts)), beta = list(.beta(combined_data$counts)))
    class <- "CoMeth1"
  } else if (m == 2L){
    assays <- c(combined_data$counts, list(EP = .EP(combined_data$counts)), LOR = list(.LOR(combined_data$counts, m)))
    class <- "CoMeth2"
  } else{
    assays <- c(combined_data$counts, list(EP = .EP(combined_data$counts)))
    class <- "CoMeth3Plus"
  }
  new(class, SummarizedExperiment(assays = assays, rowData = mtuples, 
                                  colData = colData, exptData = exptData, 
                                  verbose = verbose))
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

## cbind combines CoMeth objects from different samples and (possibly) different m-tuples
## TODO: mcols() info of each CoMeth object is lost; fix.
#' @export
setMethod("cbind", "CoMeth", function(..., deparse.level = 1){
  
  args <- unname(list(...))
  
  ## Check that all objects have the same seqinfo. 
  ## Don't want to cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being cbind-ed have the same seqinfo. 
  ## While it might be nice to to cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from hg18). Rather than handling each case in cbindI will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely cbind-ed
  seqinfo <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
  if (is(seqinfo, "try-error")){
    stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('seqinfo'), ".")
  }
  
  ## Check that all CoMeth objects have the same 'm'
  m <- sapply(args, getM)
  if (!.zero_range(m)){
    stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples, that is, different ", sQuote('m'), ".")
  }
  m <- m[1]
  
  ## Check that there are no duplicate sample_names
  sample_names <- unlist(lapply(args, sampleNames))
  if (any(duplicated(sample_names))){
    stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")
  }
  
  ## Check that all colData can be safely rbind-ed. Note that even though we are cbind-ing the CoMeth objects, we have to rbind their colData.
  ## TODO: There might be a better way to do this using GenomicRanges:::.cbind.DataFrame, which looks like it might allow for different columns in the colData. However, I would still need to ensure that all colData contain the 'methylation_type' column.
  if (!isTRUE(all(sapply(args, function(x, x1){identical(sort(colnames(colData(x))), sort(colnames(colData(x1))))}, x1 = args[[1]])))){
    stop("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('colData'), ".") 
  }
  
  ## Extract and combine all data
  sample_names <- as(unlist(lapply(args, sampleNames)), "CharacterList")
  methylation_type <- as(unlist(lapply(args, getMethylationType)), "CharacterList")
  names(methylation_type) <- sample_names
  counts <- DataFrameList(unlist(lapply(args, function(x, m){
    lapply(sampleNames(x), function(sn, x, m){
      #val <- sapply(seq_len(length(assays(x))), function(i, sn, x){
      val <- sapply(.make_m_tuple_names(m), function(an, sn, x){
        #assay(x, i)[, sn]
        assay(x, an)[, sn]
        }, sn = sn, x = x)
      colnames(val) <- .make_m_tuple_names(m)
      val <- DataFrame(val)
      return(val)}, 
      x = x, m = m)
  }, m = m)))
  names(counts) <- sample_names
  seqnames <- RleList(unlist(lapply(args, function(x){replicate(n = ncol(x), seqnames(x))})))
  names(seqnames) <- sample_names
  pos <- DataFrameList(unlist(lapply(args, function(x){replicate(n = ncol(x), DataFrame(getPos(x)))})))
  names(pos) <- sample_names
  strand <- RleList(unlist(lapply(args, function(x){replicate(n = ncol(x), strand(x))})))
  names(strand) <- sample_names
  colData <- do.call(what = "rbind", args = lapply(args, function(x){
    colData(x)[-c(which(colnames(colData(x)) == 'methylation_type'))]
  }))
  exptData <- do.call(what = "c", args = lapply(args, exptData))
  
  ## Construct CoMeth object from combined data
  CoMeth(sample_names = sample_names, methylation_type = methylation_type, counts = counts, seqnames = seqnames, pos = pos, seqinfo = seqinfo, strand = strand, colData = colData, exptData = exptData)
})

## rbind combines Cometh objects from the same samples and different m-tuples.
## rbind won't work on CoMeth objects with different mcols
#' @export
setMethod("rbind", "CoMeth", function(..., deparse.level = 1){

  args <- unname(list(...))
  
  ## Check that all objects have the same seqinfo. 
  ## Don't want to rbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed have the same seqinfo. 
  ## While it might be nice to to rbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from hg18). Rather than handling each case in rbindI will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed
  seqinfo <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
  if (is(seqinfo, "try-error")){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('seqinfo'), ".")
  }
  
  ## Check that all CoMeth objects have the same 'm'
  m <- sapply(args, getM)
  if (!.zero_range(m)){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples (i.e. different ", sQuote('m'), ").")
  }
  m <- m[1L]
  
  ## Check that all CoMeth objects have the same sample_names
  if (!isTRUE(all(sapply(args, function(x, x1){identical(sampleNames(x), sampleNames(x1))}, x1 = args[[1]])))){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('sample_names'), ".")
  }
  
  ## Check that all mcols of rowData are idential
  if (!isTRUE(all(sapply(args, function(x, x1){identical(colnames(mcols(x)), colnames(mcols(x1)))}, x1 = args[[1]])))){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('mcols'), ".")
  }
  
  ## Check that all colData is identical
  if (!isTRUE(all(sapply(args, function(x, x1){identical(colData(x), colData(x1))}, x1 = args[[1]])))){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " objects with different ", sQuote('colData'), ".") 
  }
  
  ## Check that there are no identical m-tuples
  if(sapply(seq_len(length(args) - 1), function(i, args){any(overlapsAny(args[[i]], args[[i + 1]], type = 'equal'))}, args = args)){
    stop("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")
  }
  
  ## Extract and combine all data (NEW WAY)
  sample_names <- unique(unlist(lapply(args, sampleNames)))
  methylation_type <- paste0(sort(unique(unlist(lapply(args, getMethylationType)))), collapse = '/')
  colData <- colData(args[[1]]) # colData has already been checked and found to be identical for all args
  mtuples <- do.call("c", lapply(args, rowData))
  assays <- SimpleList(lapply(names(assays(args[[1]])), function(i, args){
    do.call("rbind", lapply(args, function(y, i){
      assay(y, i)
    }, i = i))
  }, args = args))
  names(assays) <- names(assays(args[[1]]))
  exptData <- do.call(what = "c", args = lapply(args, exptData))
  
  ## Don't call the CoMeth constructor but because there's no need to "combine"
  ## data, just need to initialise the class.
  if (m == 1L){
    class <- "CoMeth1"
  } else if (m == 2L){
    class <- "CoMeth2"
  } else{
    class <- "CoMeth3Plus"
  }
  
  ## TODO: REMOVE print statement
  print("Trying to creating class...")
  new(class, SummarizedExperiment(assays = assays, rowData = mtuples, colData = colData, exptData = exptData))
})

## combine tries to figure out the combination of cbind and rbind that will automatically combine the CoMeth objects into a single CoMeth object.
## At this stage, only simple combinations are covered and more complex combinations will require manual handling.
## I'll try to add more advanced combining as the need arises.
#' @export
setMethod("combine", "CoMeth", function(x, y, ...){
  
  if (class(x) != class(y)){
    stop(paste("objects must be the same class, but are ", class(x), ", ", class(y), sep = ""))
  }
  
  if (!identical(getM(x), getM(y))){
    stop("Cannot ", sQuote('combine'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples (i.e. different ", sQuote('m'), ").")
  }
  
  ## Try to rbind if objects contain the same sample names
  if (identical(sampleNames(x), sampleNames(y))){
    val <- try(rbind(x, y), silent = TRUE)
    if (is(val, "try-error")){
      stop(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
    }
  ## Try to cbind if objects contain different sample names
  } else if (!isTRUE(any(duplicated(c(sampleNames(x), sampleNames(y)))))){
    val <- try(cbind(x, y), silent = TRUE)
    if (is(val, "try-error")){
      stop("Sorry, ", sQuote('combine'), " failed when trying to ", sQuote('cbind'), " the intermediate ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
    }
  ## Otherwise we can't combine them automatically :(
  } else if (isTRUE(any(sampleNames(x) %in% sampleNames(y)))){
    stop("Sorry, cannot automatically ", sQuote('combine'), " this combination of ", sQuote('CoMeth'), " objects. You might try to manually combine these ", sQuote('CoMeth'), " objects one-by-one using the ", sQuote('rbind'), " and ",  sQuote('cbind'), " methods.\nPlease see the help page for CoMeth, which can accessed by typing ", sQuote("?CoMeth"), " at the R prompt, for further details of combining multiple ", sQuote('CoMeth'), " objects.")
  }
  return(val)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison
###

## TODO: Check the other comparison methods work (i.e. those defined for MTuples)
## compare defers to the method defined for MTuples.
#' @export
setMethod("compare", c("CoMeth", "CoMeth"), function(x, y){
  
  if (getM(x) != getM(y)){
    stop("Cannot ", sQuote('compare'), " ", sQuote('CoMeth'), " objects with different ", sQuote('m'))
  }
  compare(rowData(x), rowData(y))
})

#' @export
setMethod(order, "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
  
  args <- list(...)
  
  if (!.zero_range(sapply(args, getM))){
    stop("All ", sQuote('CoMeth'), " objects must have the same ", sQuote('m'), " value")
  }
  
  args <- lapply(list(...), rowData)
  do.call("order", c(args, list(na.last = na.last, decreasing = decreasing)))  
})
  


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Filtering
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read .tsv file(s) from comethylation and construct a CoMeth object
###

#' Parse \code{tsv} output files from the Python software \code{comethylation}.
#'
#' Read the \code{tsv} output files from the Python software \code{comethylation} and construct a \code{\link{CoMeth}} object.
#'
#' @param files The \code{.tsv} files created by \code{comethylation}. These files may be compressed with gzip or bzip2.
#' @param sample_names The sample names of each file.
#' @param methylation_types A character vector with the type of methylation type of each file.
#' @param quiet logical: if \code{FALSE} (default), \code{read.comethylation()} will print a line, saying how many items (m-tuples) have been read per file.
#' @param seqinfo A \code{\link[GenomicRanges]{Seqinfo}} object containing information about the reference genome of the sample.
#' @note The compression of \code{file} is determined by the file extension: \code{.gz} for compressed with gzip and \code{.bz2} for compressed with bzip2. Otherwise the file is assumed to be uncompressed.
#'
#' @details All files should be of the same size m-tuples, e.g. all 2-tuples.
#'
#' @export
#' @seealso \code{\link{CoMeth}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
read.comethylation <- function(files, sample_names, methylation_types, seqinfo, quiet = FALSE) {
  
  ## Check seqinfo is supplied
  if (missing(seqinfo)){
    stop("Must supply ", sQuote('seqinfo'))
  }
  
  # Read in each file and store in a list
  x <- lapply(files, function(file, quiet){
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
  }, quiet = quiet)
  
  # Convert list of files into objects that can be passed to the CoMeth constructor
  names(x) <- sample_names
  sample_names <- as(sample_names, "CharacterList")
  methylation_type <- as(methylation_types, "CharacterList")
  names(methylation_type) <- sample_names
  
  CoMeth(sample_names = sample_names, methylation_type = methylation_type, counts = DataFrameList(lapply(X = x, function(y){y$counts})), seqnames = RleList(lapply(X = x, function(y){y$seqnames})), pos = DataFrameList(lapply(X = x, function(y){y$pos})), seqinfo = seqinfo)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

# TODO: Function collapseStrand(). Should only work with symmetric methylation_type, e.g. CG, CHG. Will be hard to do for more complication methylation_type, e.g. CG/CHG, without using the reference genome (slow).


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

## TODO: Include methylation_type and 'm' when show-ing a CoMeth object

## The show method is adapted from that of SummarizedExperiment
#' @export
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