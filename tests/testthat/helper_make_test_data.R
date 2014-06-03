### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to makes test data for the MTuples or CoMeth constructor. 
###
#' Function to makes test data for the MTuples or CoMeth constructor.
#' 
#' @param m The size of the m-tuples.
#' @param n The number of m-tuples.
#' @param s The number of samples.
#' @param sim_counts Whether to simulate counts (i.e. if wanting to simulate 
#' CoMeth constructor) or not (i.e. if wanting to simulate for MTuples 
#' constructor).
#' 
#' @return A list of arguments for the MTuples constructor (if sim_counts = 
#' FALSE) or a list of arguments for the CoMeth constructor (if sim_counts = 
#' TRUE).
make_MTuples_or_CoMeth_data <- function(m, n, s, sim_counts = TRUE) {
  
  ## Only simulate a single sample if using sim_counts = FALSE otherwise the number of returned m-tuples is not n
  if (!sim_counts){
    s <- 1
  }
  
  val <- lapply(seq_len(s), FUN = function(ss, sim_counts = sim_counts){
    p <- c(round(0.6 * n, 0), round(0.3 * n, 0), round(0.1 * n, 0))
    seqnames <- Rle(rep(paste0('chr', c(1, 2, 'X')), times = p))
    pos <- DataFrame(matrix(sort(sample(1:(n * m * 10), m * n, replace = FALSE)), ncol = m, byrow = TRUE, dimnames = list(NULL, paste0('pos', 1:m))))
    if (sim_counts){
      counts <- DataFrame(matrix(rpois(2^m * n, 4), ncol = 2^m, dimnames = list(NULL, sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')})))))))
      val <- list(seqnames = seqnames, pos = pos, counts = counts)
    } else{
      val <- list(seqnames = seqnames, pos = pos)
    }
    return(val)
  }, sim_counts = sim_counts)
  names(val) <- paste0('sample', seq_len(s))
  
  seqinfo <- Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), seqlengths = c(249250621, 243199373, 155270560), genome = 'hg19')
  if (sim_counts){
    methylation_type <- as(rep('CG', length(val)), "CharacterList")
    names(methylation_type) <- names(val)
    seqnames <- RleList(lapply(val, FUN = function(x){x$seqnames}))
    pos <- DataFrameList(lapply(val, FUN = function(x){x$pos}))
    counts <- DataFrameList(lapply(val, FUN = function(x){x$counts}))
    sample_names <- as(names(val), "CharacterList")
    val <- list(sample_names = sample_names, methylation_type = methylation_type, counts = counts, seqnames = seqnames, pos = pos, seqinfo = seqinfo)
    #val <- CoMeth(sample_names = sample_names, methylation_type = methylation_type, counts = counts, seqnames = seqnames, pos = pos, seqinfo = seqinfo, verbose = TRUE)
  } else{
    seqnames <- val[[1]]$seqnames
    pos <- as.matrix(val[[1]]$pos)
    val <- list(seqnames = seqnames, pos = pos, seqinfo = seqinfo)
    #val <- MTuples(seqnames = seqnames, pos = pos, seqinfo = seqinfo)
  }
  
  return(val)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to makes test data for the MethylationLociSet constructor.
###
#' Function to makes test data for the MethylationLociSet constructor.
#' 
#' @param n The number of methylation loci.
#' 
#' @return A list of arguments for the MethylationLociSet constructor.
make_MethylationLociSet_data <- function(n) {
  
  p <- c(round(0.6 * n, 0), round(0.3 * n, 0), round(0.1 * n, 0))
  seqnames <- Rle(rep(paste0('chr', c(1, 2, 'X')), times = p))
  ranges <- IRanges(start = c(sample(249250621, p[1]), 
                              sample(243199373, p[2]), 
                              sample(155270560, p[3])), 
                    width = 1)
  strand <- sample(c('+', '-'), n, replace = TRUE)
  seqinfo <- Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), 
                     seqlengths = c(249250621, 243199373, 155270560), 
                     genome = 'hg19')
  
  val <- list(seqnames = seqnames, ranges = ranges, strand = strand, 
              seqinfo = seqinfo)  
  
  return(val)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Function to makes test data for the betaCor function.
###

#' Function to makes test data for the betaCor function.
#' 
#' @param mls A MethylationLociSet.
#' @param f The fraction of loci to sample from mls.
#' @param s The number of samples.
#'
#' @return A CoMeth1 object
#' 
make_CoMeth1_data <- function(mls, f, s) {

  n <- floor(length(mls) * f)
  i <- sort(sample(length(mls), n))
  rd <- MTuples(seqnames = seqnames(mls)[i], 
                pos = matrix(start(mls)[i], ncol = 1), 
                strand = strand(mls)[i],
                seqinfo = seqinfo(mls))
  val <- lapply(seq_len(s), function(i, n) {
    M <- matrix(rpois(n = n, lambda = 10), ncol = 1)
    U <- matrix(rpois(n = n, lambda = 3), ncol = 1)
    list(M = M, U = U)
  }, n = n)
  M <- cbind(sapply(val, function(x) x$M))
  colnames(M) <- paste0('s', seq_len(s))
  U <- cbind(sapply(val, function(x) x$U))
  colnames(U) <- paste0('s', seq_len(s))
  beta <- .beta(list(M = M, U = U))
  EP <- .EP(list(M = M, U = U))
  assays = SimpleList(M = M, U = U, beta = beta, EP = EP)
  cd <- DataFrame(methylation_type = rep('CG', s))
  rownames(cd) <- paste0('s', seq_len(s))
  
  new("CoMeth1", SummarizedExperiment(rowData = rd, assays = assays, 
                                      colData = cd))
}