#### TODOs ####
# TODO: Make documentation more like GRanges, i.e. ?CoMeth describes the class and there is a subsection titled "Construction" that described the constructor, etc.
# TODO: cbind, rbind, combine, comparse, %<=% (and similar). Require that 'm', 'seqinfo', 'methylation_type' are identical when rbind-ing CoMeth objects. Allow for multiple samples in a CoMeth object
# TODO: Function collapseStrand(). Should only work with symmetric methylation_type, e.g. CG, CHG. Will be hard to do for more complication methylation_type, e.g. CG/CHG, without using the reference genome (slow). Add the 'collapse_strand' option to CoMeth constructor.
# TODO: Perhaps simplify the description of the 'pos' and 'counts' parameters to note that these can be list of data.frame or data.frame-like objects (DataFrame) since a data.frame is just a list of lists.
# TODO: Check for any orphan TODOs

#### CoMeth class definition, constructor and validity function ####
#' The CoMeth class
.CoMeth <- setClass("CoMeth", contains = "SummarizedExperiment")

# TODO: Think about whether 'methylation_type' should really be combined or should an error be thrown if not all samples in a CoMeth object have the same 'methylation_type'?
#' The constructor function for CoMeth objects
#'
#' Users will generally use the \code{\link{read.comethylation}} to initialise a \code{CoMeth} object, however, this is the general method to initialise a \code{CoMeth} object. The initialisation allows for multiple samples in the same object by passing most arguments as a list where each element of the list corresponds to the arguments for a given sample. These multiple samples will then be appropriately combined into a single \code{CoMeth} object.
#' @param m An \code{integer} storing the size of the m-tuples, i.e. the \code{m} in m-tuple. Only a single value is accepted, and not a list, because \code{m} must be the same for all samples in a \code{CoMeth} object.
#' @param sample_names A \code{character} vector containing the names of the samples. Sample names must be unique.
#' @param pos A \code{\link[IRanges]{DataFrameList}} \code{pos} must be named and the names must match those given in \code{sample_names}. The columns of each DataFrame must be: \code{seqnames}, \code{pos1}, ..., \code{posm}, where, for example, \code{posm} is \code{pos3} if \code{m} = 3. \code{seqnames} stored the sequence/chromosome name of the m-tuples. Therefore, the number of columns of each DataFrame is \code{m} + 1 and the number of rows is equal to the number of m-tuples for that particular sample.
#' @param counts A \code{\link[IRanges]{DataFrameList}}. \code{counts} must be named and the names must match those given in \code{sample_names}. The entry in each DataFrame corresponds to the number of times that particular co-methylation pattern (columns) was observed for that particular m-tuple (rows). Therefore, each DataFrame must have the same number of rows as its corresponding DataFrame in \code{pos} and have \eqn{2 ^ \code{m}} columns. The column names of each DataFrame must match that given by \code{make_m_tuple_names(m)}. 
#' @param strand An (optional)\code{\link[IRanges]{RleList}} object containing the strand information of each m-tuple. \strong{WARNING}: If \code{strand} is missing, all m-tuples in the resulting \code{\link{CoMeth}} object will have their strand set to \code{*} to signify that the strand is unknown or irrelevant (such as when methylation measurements have been combined across strands). \strong{WARNING}: m-tuples will not be combined across samples if they are on different strands. 
#' @param methylation_type A \code{\link[IRanges]{CharacterList}} storing the type of methylation loci for these m-tuples. \code{methylation_type} must be named and the names must match those given in \code{sample_names}. For each sample the possible values are "CG", "CHG", "CHH" or "CNN" or multiple values specified as a character vector, e.g. c("CG", "CHG") or c("CHG", "CG") for "CG/CHG" methylation. The \code{methylation_type} of the resulting \code{CoMeth} object is the union of this argument, e.g. if \code{methylation_type = list('sample1' = 'CG', 'sample2' = 'CHH')} then the \code{methylation_type} of the resulting \code{CoMeth} object is \code{CG/CHH}.
#' @param seqinfo A \code{\link[GenomicRanges]{Seqinfo}} object containing information about the reference genome of the samples. Only a single value is accepted, and not a list, because all samples must be mapped against the same reference genome. However, multiple genomes per sample are allowed to accommodate a spiked-in unmethylated genome (normally lambda phage), which is a common step in a bisulfite-sequencing protocol.
#' @param sort_cometh \code{logical}. Should the resulting \code{CoMeth} object be sorted by genomic co-ordinates? \strong{WARNING}: CoMeth objects are sorted by seqnames (sort order in accordance with \code{seqinfo}), strand (sort order: \code{+} > \code{-} > \code{*}) and finally by m-tuples. m-tuples are ordered element-wise left-to-right, so that \code{chr1:*:(3, 5, 11)} < \code{chr1:*:(3, 7, 9)} < \code{chr1:*:(3, 7, 11)}.
#'
#' @export
#' @seealso \code{\link{read.comethylation}}
#' @return A \code{\link{CoMeth}} object
#' @examples
#' cat("TODO")
CoMeth <- function(m, sample_names, pos, counts, strand, methylation_type, seqinfo, sort_cometh = TRUE){
  
  ## Check 'sample_names'
  if (missing(sample_names) || !identical(unique(sample_names), sample_names)){
    stop("Need 'sample_names'. Must be a character vector where each element is a unique sample name.")
  }
  
  ## Check 'm'
  ## Must be checked before 'counts' is checked because it uses make_m_tuple_names(m)
  if (missing(m) || !is.integer(m) || length(m) > 1 || m < 1){
    stop("'m' must be specified and must be a single, positive integer.")
  }
  
  ## Check 'pos'
  if (missing(pos)){
    stop("'pos' missing. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!is(pos, "DataFrameList")){
    stop("'pos' must be a DataFrameList. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (is.null(names(pos)) || !all(names(pos) %in% sample_names) || length(names(pos)) != length(unique(names(pos)))){
    stop("Names of 'pos' must match those in 'sample_names'. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!identical(colnames(pos)[[1]], c('seqnames', paste0('pos', seq_len(m)))) || any(sapply(colnames(pos), function(x, y){x != y}, y = colnames(pos)[[1]]))){
    stop("'m' is set to ", m, " so colnames for all elemnts of 'pos' must be: 'seqnames', ", paste0('pos', seq_len(m), collapse = "', '"), "'.\nPlease see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!all(sapply(pos, function(x, m){isTRUE(all.equal(sapply(x, class), c('Rle', rep('integer', m)), check.attributes = F))}, m = m))){
    stop("Class of columns in each DataFrame element of 'pos' must be: '", paste0(c('Rle', rep('integer', m)), collapse = "', '"), "'.\nPlease see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  
  ## Check 'counts'
  if (missing(counts)){
    stop("'counts' missing. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!is(counts, "DataFrameList")){
    stop("'counts' must be a DataFrameList. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (is.null(names(counts)) || !all(names(counts) %in% sample_names) || length(names(counts)) != length(unique(names(counts)))){
    stop("Names of 'counts' must match those in 'sample_names'. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (any(nrow(counts) != nrow(pos))){
    stop("nrow(counts) must be equal to nrow(pos).  Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of these arguments.")
  }
  if (!identical(colnames(counts)[[1]],  make_m_tuple_names(m)) || any(sapply(colnames(counts), function(x, y){x != y}, y = colnames(counts)[[1]]))){
    stop("'m' is set to ", m, " so colnames for all elements of 'counts' must be: '", paste0(make_m_tuple_names(m), collapse = "', '"), "'.\nPlease see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!all(sapply(counts, function(x, m){isTRUE(all.equal(sapply(x, class), rep('integer', 2 ^ m), check.attributes = F))}, m = m))){
    stop("Class of columns in each DataFrame element of 'counts' must be: '", paste0(rep('integer', 2 ^ m), collapse = "', '"), "'.\nPlease see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  
  ## Check 'strand'
  if (missing(strand)){
    warning("'strand' missing so setting to '*' for all m-tuples")
    strand <- RleList(lapply(sample_names, function(i, pos){Rle('*', nrow(pos[[i]]))}, pos = pos))
    names(strand) <- sample_names
  } 
  if (!is(strand, "RleList")){
    stop("'strand' must be an RleList. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!identical(sapply(strand, length), sapply(pos, nrow))){
    stop("Length of each element in 'strand' must equal the number of rows of its corresponding element in 'pos'.")
  }
  if (c("*", "+", "-") %in% as.vector(unlist(strand, use.names = FALSE)) %*% c(1, 0.5, 0.5) > 1){
    warning("'strand' contains '*' as well as at least one of '+' or '-'. m-tuples will not be combined across samples if they are on different strands. This means that m-tuples with the exact same genomic co-ordinates across samples, but with different strands (e.g. '*' and '+'), will not be combined.\nTo combine all m-tuples based on genomic co-ordinates regardless of strand, do not use the 'strand' parameter when calling the CoMeth constructor.")
  }
  
  ## Check 'methylation_type'
  if (missing(methylation_type)){
    stop("'methylation_type' missing. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!is(methylation_type, 'CharacterList')){
    stop("'methylation_type' must be a CharacterList. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (is.null(names(methylation_type)) || !all(names(methylation_type) %in% sample_names) || length(names(methylation_type)) != length(unique(names(methylation_type)))){
    stop("Names of 'methylation_type' must match those in 'sample_names'. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!all(unlist(methylation_type) %in% c('CG', 'CHG', 'CHH', 'CNN'))){
    stop("'methylation_type' for each sample must be 'CG', 'CHG', 'CHH' and 'CNN' or a vector of some combination of these, e.g. 'c('CG', 'CHG')'.\nPlease see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  } 

  ## Check 'seqinfo'
  if (missing(seqinfo)){
    stop("'seqinfo' missing. Please see the help page for CoMeth, which can accessed by typing '?CoMeth' at the R prompt, for further details of this argument.")
  }
  if (!is(seqinfo, 'Seqinfo')){
    stop("'seqinfo' must be a Seqinfo object.")
  }
  
  ## Check 'sort_cometh'
  if (!isTRUEorFALSE(sort_cometh)){
    stop("'sort_cometh' must be TRUE or FALSE")
  }
  
  ## Combine the list-wise data into matrix-like data whilst taking care of common and sample-specific m-tuples, i.e. filling in zeros when a sample doesn't have any observations for that m-tuple.
  methylation_type <- paste0(sort(unique(unlist(methylation_type))), collapse = '/')
  combined_data <- .combine(m = m, sample_names = sample_names, pos = pos, counts = counts, strand =  strand)
  
  ## Construct rowData of CoMeth object
  if (m > 2){
    # Need to store other positions if m > 2
    gr <- GRanges(seqnames = combined_data[['coordinates']][, 1], ranges = IRanges(start = combined_data[['coordinates']][, 2], end = combined_data[['coordinates']][, m + 1]), strand = combined_data[['coordinates']][, m + 2], seqinfo = seqinfo, combined_data[['coordinates']][, -c(1, 2, m + 1, m + 2), drop = FALSE])
    } else{
      gr <- GRanges(seqnames = combined_data[['coordinates']][, 1], ranges = IRanges(start = combined_data[['coordinates']][, 2], end = combined_data[['coordinates']][, m + 1]), strand = combined_data[['coordinates']][, m + 2], seqinfo = seqinfo)
    }
  
  ## Construct assays of CoMeth object
  assays <- combined_data[['counts']]
  
  ## Construct colData of CoMeth object
  colData <- DataFrame(m = rep(m, length(sample_names)), methylation_type = rep(paste0(sort(methylation_type), collapse = '/'), length(sample_names)), row.names = sample_names)
  
  ## Construct CoMeth object
  cometh <- SummarizedExperiment(assays = assays, rowData = gr, colData = colData)
  cometh <- .CoMeth(cometh)

  ## Sort CoMeth object if required
  if (sort_cometh){
    cometh <- sort(cometh)
  }
  
  ## Return CoMeth object
  return(cometh)
}

## CoMeth validity function
validCoMeth <- function(object){
  ## Check that m is identical for all samples
  ## Can't use getM(object) because getM checks whether object is a valid CoMeth object which it isn't when it gets called by the CoMeth constructor
  
  msg <- NULL
  
  ## Check m is correctly set and makes sense
  m <- colData(object)[['m']][1] 
  if (!zero_range(colData(object)[['m']])){
    msg <- validMsg(msg, "'m' must be identical for all samples")
  }
  if (!is.integer(m) || m < 1){
    msg <- validMsg(msg, "'m' must be a positive integer")
  }
  
  ## Check that there are 2^m assays
  if (log2(length(assays(object))) != m){
    msg <- validMsg(msg, "length(assays) does not equal 2 ^ m")
  }
  
  ## Check that there are m positions. Have to avoid call to getPos to avoid recursive call to validObject. It's not easy to define a quick way of doing this, particularly when m = 1 or 2 because then there are no "extra positions" from which to infer the "true m".
  if (m == 1){
    if (any(grepl('pos', colnames(mcols(object))))){
      msg <- validMsg(msg, paste0("'m' is ", m, " but 'inferred m' = ",  2 + sum(grepl('pos', colnames(mcols(object)))), " based on positions of m-tuples")) 
    } else if (!isTRUE(all.equal(start(object), end(object)))){
      msg <- validMsg(msg, paste0("'m' is ", m, " but 'inferred m' = 2 based on positions of m-tuples"))
    }
  } else if (m == 2){
      if (any(start(object) == end(object))){
        msg <- validMsg(msg, paste0("'m' is ", m, " but 'inferred m' = 1 based on positions of m-tuples"))
      } else if (any(grepl('pos', colnames(mcols(object))))){
        msg <- validMsg(msg, paste0("'m' is ", m, " but 'inferred m' = ",  2 + sum(grepl('pos', colnames(mcols(object)))), " based on positions of m-tuples")) 
      }
  } else{
    if (sum(grepl('pos', colnames(mcols(object)))) != (m - 2)){
      msg <- validMsg(msg, paste0("'m' is ", m, " but 'inferred m' = ",  2 + sum(grepl('pos', colnames(mcols(object)))), " based on positions of m-tuples")) 
    }
  }

  
  ## Check assay names
  msg <- validMsg(msg, .checkAssayNames(object, make_m_tuple_names(m)))
  
  ## Check colData contains 'm' and 'methylation_type' once and only once. Other colData is allowed
  if (!all(sapply(c('^m$', '^methylation_type$'), function(x, y){sum(grepl(x, y)) == 1}, y = colnames(colData(object))))){
    stop("colData of 'CoMeth' must contain columns 'm' and 'methylation_type' once each, and once each only.")
  }
  
  ## Check rowData is GRanges
  if (class(rowData(object)) != "GRanges"){
    msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
  }
  
  ## Check that all 'counts' are non-negative
  ## Note from bsseq: benchmarking shows that min(assay()) < 0 is faster than any(assay() < 0) if it is false
  if (min(sapply(assays(object), min, na.rm = TRUE), na.rm = TRUE) < 0) {
    msg <- validMsg(msg, "'counts' has negative entries")
  }
  
  ## Check that all 'pos' are non-negative
  if (m == 1){
    if (min(start(object) < 0)){
      msg <- validMsg(msg, "'pos' has negative entries")
    }
  } else if (m == 2){
    if (min(start(object) < 0) || min(end(object) < 0)){
      msg <- validMsg(msg, "'pos' has negative entries")
    }
  } else{
    if ((min(start(object)) < 0) || (min(end(object)) < 0) || (min(mcols(object)[, grepl('pos', colnames(mcols(object)))]) < 0)){
      msg <- validMsg(msg, "'pos' has negative entries")
    }
  }
 
  ## Return validity of object
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}
setValidity("CoMeth", validCoMeth)

#### Methods for CoMeth objects ####

setMethod(show, "CoMeth", function(object){
  cat("An object of type 'CoMeth' with:\n")
  cat(paste0(" ", nrow(object), " ", getMethylationType(object), " ", getM(object), "-tuples\n"))
  cat(paste0(" ", ncol(object), " samples\n"))
  #callNextMethod()
})

setMethod(length, "CoMeth", function(x){
  nrow(x)
})

## "[" method for a CoMeth object is inherited from SummarizedExperiment

## The cbind method for a CoMeth object differs to that for a SummarizedExperiment. cbind allows for the objects to have different 'pos' data (getPos()) provided all objects have the same seqinfo. Samples without 'counts' data for a particular m-tuple will have the corresponding value filled with NA. The 'm' of all objects must be identical. 'methylation_type' may differ between object but the CoMeth object will have the combined 'methylation_type', e.g. cbind-ing a "CG" with "CHG" will result in a new object with a "CG/CHG" methylation_type.  sample_names must be unique across all objects being combined. The way cbind works is effectively to extract all the relevant information from the objects being combined and then create a new CoMeth object with the CoMeth constructor.
setMethod("cbind", "CoMeth", function(..., deparse.level = 1){
  args <- unname(list(...))
  
  ## Ensure that all objects have the same seqinfo. 
  ## Don't want to rbind/cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed/cbind-ed have the same seqinfo. 
  ## While it might be nice to to rbind/cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from mm9). Rather than handling each case in rbind/cbind I will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed/cbind-ed
  seq_info <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
  if (is(seq_info, "try-error")){
    stop("Can only cbind 'CoMeth' objects with identical 'seqinfo'.")
  }
  
  ## Ensure all sampleNames are unique across the objects being combined
  sample_names <- unlist(lapply(args, sampleNames))
  if (length(sample_names) != length(unique(sample_names))){
    stop("Cannot cbind 'CoMeth' objects containing duplicate 'sample_names'.")
  }
  
  ## Ensure all objects have the same 'm'
  m <- sapply(args, getM)
  if (!zero_range(m)){
    stop("Cannot cbind 'CoMeth' objects with the different 'm'.")
  }
  m <- m[1]
  
  ## Construct 'methylation_type' (combining if necessary)
  methylation_type <- sapply(args, getMethylationType)
  if (!all(methylation_type == methylation_type[[1]])){
    warning("Combining 'CoMeth' objects with different 'methylation_type'. The 'methylation_type' of the new CoMeth object is: ", paste0(sort(unique(methylation_type)), collapse = '/'))
  }
  methylation_type <- CharacterList(lapply(sample_names, function(x){sort(unique(unlist(methylation_type)))}))
  names(methylation_type) <- sample_names
  
  ## Construct pos, counts and strand to be passed to CoMeth constructor
  pos <- DataFrameList(unlist(lapply(args, function(x){replicate(length(sampleNames(x)), getPos(x))}), recursive = FALSE))
  names(pos) <- sample_names
  counts <- DataFrameList(unlist(lapply(args, function(x){getCounts(x)}), recursive = FALSE))
  names(counts) <- sample_names  
  strand <- RleList(unlist(lapply(args, function(x){replicate(length(sampleNames(x)), strand(x))}), recursive = FALSE))
  names(strand) <- sample_names
  
  ## Use CoMeth constructor
  CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info)
})


## rbind combines objects with different m-tuples and the same subjects (sampleNames). The sampleNames must match or an error is thrown. Duplicate columns of colData must contain the same data.
## Can't simply use the rbind method defined for a SummarizedExperiment because we also need to check  the "extra positions" before concluding that an m-tuple is unique
setMethod("rbind", "CoMeth", function(..., deparse.level = 1){
  args <- unname(list(...))
  
  ## Ensure that all CoMeth objects have the same seqinfo. Don't want to rbind/cbind CoMeth objects from different genomes. This "merge" of the seqinfo will still allow for multiple genomes (e.g. human and lambda_phage) within the same CoMeth object provided all CoMeth objects that are being rbind-ed/cbind-ed have the same seqinfo. While it might be nice to to rbind/cbind CoMeth objects with different seqinfo in some cases (e.g. one object contains data on chr1 and another on chr2), it is a bad idea in other cases (e.g. one object contains data from hg19 and one from mm9). Rather than handling each case in rbind/cbind I will require to the user to update the seqinfo of the CoMeth in an appropriate way to allow for them to be safely rbind-ed/cbind-ed
  seq_info <- try(do.call("merge", lapply(args, seqinfo)), silent = TRUE)
  if (is(seq_info, "try-error")){
    stop("Can only rbind 'CoMeth' objects with identical 'seqinfo'.")
  }
  
  ## Ensure all sampleNames are identical across the objects being combined
  sample_names <- lapply(args, sampleNames)
  if (!all(sapply(sample_names, function(x, y){identical(x, y)}, y = sample_names[[1]]))){
    stop("Cannot rbind CoMeth objects with different sampleNames")
  }
  sample_names <- sample_names[[1]]
  
  ## Ensure all objects have the same 'm'
  m <- sapply(args, getM)
  if (!zero_range(m)){
    stop("Cannot rbind 'CoMeth' objects with the different 'm'.")
  }
  m <- m[1]
  
  ## Construct 'methylation_type' (combining if necessary)
  methylation_type <- sapply(args, getMethylationType)
  if (!all(methylation_type == methylation_type[1])){
    warning("Combining 'CoMeth' objects with different 'methylation_type'. The 'methylation_type' of the new CoMeth object is: ", paste0(sort(unique(methylation_type)), collapse = '/'))
  }
  methylation_type <- CharacterList(lapply(sample_names, function(x){sort(unique(methylation_type))}))
  names(methylation_type) <- sample_names

  ## Check that all 'pos' information is unique across samples
  pos <- do.call('rbind', lapply(args, getPos))
  strand <- do.call('c', lapply(args, strand))
  if (any(duplicated(cbind(pos, strand)))){
    stop("Cannot rbind 'CoMeth' objects with duplicate m-tuples")
  }
  
  ## Construct pos, counts and strand data to pass to CoMeth constructor
  pos <- DataFrameList(replicate(length(sample_names), do.call(rbind, lapply(args, getPos))))
  names(pos) <- sample_names
  counts <- DataFrameList(lapply(sample_names, function(sample_name, counts){do.call('rbind', lapply(counts, function(x, sample_name){x[[sample_name]]}, sample_name = sample_name))}, counts = lapply(args, getCounts)))
  names(counts) <- sample_names
  strand <- RleList(replicate(length(sample_names), do.call(c, lapply(args, strand))))
  names(strand) <- sample_names
  
  ## Use CoMeth constructor
  CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info)  
})

## TODO: Simplify. Just use getPos() on each object (assuming getPos returns seqnames, strand, pos1, ..., posm; where strand is a properly-ordered factor) and then use the order method defined for DataFrame
## The order method for a CoMeth
## m-tuples are ordered element-wise left-to-right, so that chr1:(3, 5, 11) < chr1:(3, 7, 9) < chr1:(3, 7, 11)
## WARNING: I don't really understand when you would use the case where "..." contains multiple objects but I've included it anyway based on selectMethod("order", "GenomicRanges")
setMethod(order, "CoMeth", function(..., na.last = TRUE, decreasing = FALSE){
# Based on selectMethod("order", "GenomicRanges")
  if (!isTRUEorFALSE(decreasing)){
    stop("'decreasing' must be TRUE or FALSE")
  }
  
  args <- list(...)
  
  if (!zero_range(sapply(args, getM))){
    stop("All 'CoMeth' objects must have the same 'm' value")
  }
  
  m <- getM(args[[1]])
  
  if (m < 3){
    ## If m < 3 just defer to the order method defined for GRanges
    do.call("order", sapply(args, rowData))
  } else{
    ## If m >= 3 then need to define an order method that takes note of the 'extra positions' in each m-tuple
    order_args <- vector("list", (m + 2L) * length(args))
    idx <- (m + 2L) * seq_len(length(args))
    order_args[seq.int(from = 1, to = max(idx), by = m + 2)] <- lapply(args, function(x){
      as.factor(seqnames(x))
      })
    order_args[seq.int(from = 2, to = max(idx), by = m + 2)] <- lapply(args, function(x){
      as.factor(strand(x))
      })
    order_args[seq.int(from = 3, to = max(idx), by = m + 2)] <- lapply(args, start)
    order_args[seq.int(from = 4, to = max(idx), by = m + 2) + rep(seq(0, m - 3, by = 1))] <- lapply(args, function(x, m){getPos(x)[, -c(1, 2, m)]}, m = m)
    order_args[seq.int(from = 5 + m - 3, to = max(idx), by = m + 2)] <- lapply(args, function(x){end})
    order_args[idx] <- lapply(args, function(x){end(x)})
    do.call(order, c(order_args, list(na.last = na.last, decreasing = decreasing)))
  }
})
  
## The sort method for a CoMeth object
setMethod(sort, "CoMeth", function(x, decreasing = FALSE, ...){
  x[order(x, decreasing = decreasing), ] 
})

## The sampleNames getter and setter for a CoMeth object
setMethod("sampleNames", "CoMeth", function(object) {
  colnames(object)
})

setReplaceMethod("sampleNames", signature = signature(object = "CoMeth", value = "ANY"), function(object, value) {
  if (length(value) != length(sampleNames(object))){
    stop("Invalid 'sampleNames' length")
  }
  colnames(object) <- value
  object
})

#### Functions that work on CoMeth objects ####
## TODO: Should I return the strand along with the seqnames and pos? If so, should the function really be called getCoordinates>
#' Obtain the \code{pos} from a \code{\link{CoMeth}} object
#' 
#' The \code{pos} of a \code{\link{CoMeth}} object are the genomic coordinates (seqname and positions) of the cytosines that make up each m-tuple.
#' @param x A CoMeth object
#'
#' @return A DataFrame containing the seqnames, strand and positions (columns) of each m-tuple (rows)
#' @note \code{getPos} does not return the \code{strand} of the \code{\link{CoMeth}} object. The \code{strand} can be obtained by \code{strand(x)} where \code{x} is a \code{\link{CoMeth}} object.
#' @export
getPos <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  m <- getM(x)
  if (m == 1){
    pos  <- DataFrame(seqnames(x), start(x))
  } else if (m == 2) {
    pos <- DataFrame(seqnames(x), start(x), end(x))
  } else {
    pos <- DataFrame(seqnames(x), start(x), mcols(rowData(x))[, paste0('pos', seq(from = 2, to = m - 1, by = 1))], end(x))
  }
  colnames(pos) <- c('seqnames', paste0('pos', seq_len(m))) 
  return(pos)
}

#' Obtain the \code{counts} from a \code{\link{CoMeth}} object
#' 
#' The \code{counts} of a \code{\link{CoMeth}} object are the number of each co-methylation pattern is observed at each m-tuple (rows) for each sample (columns)
#' @param x A CoMeth object
#' @param sample_names An optional character vector listing for which samples the counts data should be returned. If not specified the counts data for all samples in the CoMeth object are returned.
#'
#' @return A list of matrices. Each element of the list is the counts data for a single sample. Each matrix contains the number of times each co-methylation pattern is observed at each m-tuple for that sample.
#' @export
getCounts <- function(x, sample_names) {
  stopifnot(is(x, "CoMeth"))
  
  if (missing(sample_names)){
    sample_names <- sampleNames(x)
  }
  
  if (any(!sample_names %in% sampleNames(x))){
    stop("Could not find 'counts' data for sample_names '", paste0(sample_names[which(!sample_names %in% sampleNames(x))], collapse = ', '), "' in 'x'")
  }
  
  counts <- lapply(sample_names, function(y, x){sapply(assayNames(x), function(z, y, x){assay(x, z)[, y, drop = FALSE]}, y = y, x = x)}, x = x)
  names(counts) <- sample_names
  
  return(counts)
}

#' Get \code{m} from a \code{\link{CoMeth}} object
#'
#' Get the size of m-tuples, i.e. the \code{m} in m-tuples, from a \code{\link{CoMeth}} object.
#' @param x A CoMeth object
#'
#' @return An integer
#' @export
getM <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  m <- colData(x)[, 'm']
  if (!zero_range(m)){
    stop("'m' not equal for all samples in 'CoMeth' object")
  } else{
    return(m[1])
  }
}

#' Get \code{methylation_type} from a \code{\link{CoMeth}} object
#'
#' Get the size of methylation type, e.g. 'CG', 'CHH', 'CG/CHG', etc., from a \code{\link{CoMeth}} object.
#' @param CoMeth A CoMeth object
#'
#' @return An string describing the \code{methylation_type}
#' @export
getMethylationType <- function(x) {
  stopifnot(is(x, "CoMeth"))
  
  mt <- colData(x)[['methylation_type']][1]
  if (!all(mt == mt[1])){
    stop("'methylation_type' not equal for all samples in 'CoMeth' object")
  } else{
    return(mt[1])
  }
}

#' Get the \code{coverage} from a \code{\link{CoMeth}} object.
#'
#' The \code{coverage} of a \code{\link{CoMeth}} object is the number of reads for each m-tuple, i.e. \code{rowSums(getCounts(CoMeth))}. The current implementation is slow.
#' @param CoMeth A \code{\link{CoMeth}} object
#'
#' @return A numeric matrix. Each column of the matrix corresponds to a sample and each row to an m-tuple.
#' @export
getCoverage <- function(x) {
  stopifnot(is(x, "CoMeth"))
  Reduce("+", assays(x))
}
