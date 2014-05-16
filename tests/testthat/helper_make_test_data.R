## Function to create test data. Returns a list of arguments for the MTuples constructor (if sim_counts = FALSE) or a list of arguments for the CoMeth constructor (if sim_counts = TRUE)
make_test_data <- function(m, n, s, sim_counts = TRUE){
  
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