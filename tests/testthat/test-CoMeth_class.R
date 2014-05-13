## TODO: It's possible to have two SummarizedExperiments (well, CoMeth) objects that have identical slots, all.equal returns TRUE but identical returns FALSE. 
## E.g. a <- make_test_data(3L, 10L, 3L, sim_counts = TRUE), identical(combine(a[, 1], a[, 2], a[, 3]), a)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Define test data.
###
context('Define test data')

## TODO: Move make_test_data() to its own R file and call it in both test-MTuples_class.R and test-CoMeth_class.R. Relatedly, figure out the best way to call this function from these test scripts.
expect_true(FALSE)

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth constructor.
###
context('CoMeth constructor')

test_that("CoMeth constructor returns an CoMeth object when m = 1 for a single sample", {
  m1_s1 <- make_test_data(m = 1L, n = 10L, s = 1L, sim_counts = TRUE)
  m1_s1 <- CoMeth(sample_names = m1_s1$sample_names, methylation_type = m1_s1$methylation_type, counts = m1_s1$counts, seqnames = m1_s1$seqnames, pos = m1_s1$pos, seqinfo = m1_s1$seqinfo)
  expect_that(m1_s1, is_a("CoMeth"))
  rm(m1_s1)
})

test_that("CoMeth constructor returns an CoMeth object when m = 1 for multiple samples", {
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(m1_s3, is_a("CoMeth"))
  rm(m1_s3)
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for a single sample", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(m2_s1, is_a("CoMeth"))
  rm(m2_s1)
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for multiple samples", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(m2_s3, is_a("CoMeth"))
  rm(m2_s3)
})

test_that("CoMeth constructor returns an CoMeth object when m >= 2 for a single sample", {
  m3_s1 <- make_test_data(m = 3L, n = 10L, s = 1L, sim_counts = TRUE)
  m3_s1 <- CoMeth(sample_names = m3_s1$sample_names, methylation_type = m3_s1$methylation_type, counts = m3_s1$counts, seqnames = m3_s1$seqnames, pos = m3_s1$pos, seqinfo = m3_s1$seqinfo)
  expect_that(m3_s1, is_a("CoMeth"))
  rm(m3_s1)
  m4_s1 <- make_test_data(m = 4L, n = 10L, s = 1L, sim_counts = TRUE)
  m4_s1 <- CoMeth(sample_names = m4_s1$sample_names, methylation_type = m4_s1$methylation_type, counts = m4_s1$counts, seqnames = m4_s1$seqnames, pos = m4_s1$pos, seqinfo = m4_s1$seqinfo)
  expect_that(m4_s1, is_a("CoMeth"))
  rm(m4_s1)
})

test_that("CoMeth constructor returns an CoMeth object when m >= 2 for multiple samples", {
  m3_s3 <- make_test_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(sample_names = m3_s3$sample_names, methylation_type = m3_s3$methylation_type, counts = m3_s3$counts, seqnames = m3_s3$seqnames, pos = m3_s3$pos, seqinfo = m3_s3$seqinfo)
  expect_that(m3_s3, is_a("CoMeth"))
  rm(m3_s3)
  m4_s3 <- make_test_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(sample_names = m4_s3$sample_names, methylation_type = m4_s3$methylation_type, counts = m4_s3$counts, seqnames = m4_s3$seqnames, pos = m4_s3$pos, seqinfo = m4_s3$seqinfo)
  expect_that(m4_s3, is_a("CoMeth"))
  rm(m4_s3)  
})

test_that("CoMeth parameter checking works: 'sample_names'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘sample_names’ missing."))
  expect_that(CoMeth(sample_names = unlist(m2_s3$sample_names), methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘sample_names’ must be a ‘CharacterList’."))
  expect_that(CoMeth(sample_names = as(rep('sample_1', 3), "CharacterList"), methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Each element of ‘sample_names’ must be unique."))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'methylation_type'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘methylation_type’ missing."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = as.vector(m2_s3$methylation_type), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘methylation_type’ must be a ‘CharacterList’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = unlist(m2_s3$methylation_type), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘methylation_type’ must be a ‘CharacterList’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type[1:2], counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Names of ‘methylation_type’ must match those in ‘sample_names’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = 'CpG', sample3 = "CHG"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘methylation_type’ for each sample must be ‘CG’, ‘CHG’, ‘CHH’ or ‘CNN’, or a vector of some combination of these, e.g., ‘c\\('CG', 'CHG'\\)’"))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = c('CG', 'CHH'), sample3 = "CHG"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), gives_warning("The supplied ‘methylation_type’ parameter says that samples contain data for different methylation types. The union of these methylation types will be used as the ‘methylation_type’ of the returned ‘CoMeth’ object."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = c('CG', 'CHH'), sample3 = "CHH"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), gives_warning("The supplied ‘methylation_type’ parameter says that samples contain data for different methylation types. The union of these methylation types will be used as the ‘methylation_type’ of the returned ‘CoMeth’ object."))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'counts'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘counts’ missing."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = as.list(m2_s3$counts), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘counts’ must be a ‘DataFrameList’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts[1:2], seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Names of ‘counts’ must match those in ‘sample_names’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){cbind(x[, -3])}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘ncol\\(counts\\)’ must be identical for all elements of ‘counts’ and should be a power of 2."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = c(m2_s3$counts[1:2], DataFrameList(sample3 = m2_s3$counts[[3]][, -3])), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘ncol\\(counts\\)’ must be identical for all elements of ‘counts’ and should be a power of 2."))  
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){DataFrame(lapply(x, as.numeric))}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Class of columns in each ‘DataFrame’ element of ‘counts’ must be: ‘integer’, ‘integer’, ‘integer’, ‘integer’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){colnames(x) <- rev(colnames(x)); x}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘m’ is set to 2 so colnames for all elements of ‘counts’ must be: ‘MM’, ‘MU’, ‘UM’, ‘UU’."))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'seqnames'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘seqnames’ missing."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = as.list(m2_s3$seqnames), pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("‘seqnames’ must be an ‘RleList’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames[1:2], pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Names of ‘seqnames’ must match those in ‘sample_names’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = endoapply(m2_s3$seqnames, function(x){x[-1]}), pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error("Length of each ‘Rle’ element of ‘seqnames’ must be identical to ‘nrow\\(pos\\)’"))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'pos'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, seqinfo = m2_s3$seqinfo), throws_error("‘pos’ missing."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = as.list(m2_s3$pos), seqinfo = m2_s3$seqinfo), throws_error("‘pos’ must be a ‘DataFrameList’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos[1:2], seqinfo = m2_s3$seqinfo), throws_error("Names of ‘pos’ must match those in ‘sample_names’."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = endoapply(m2_s3$pos, function(x){x[-1, ]}), seqinfo = m2_s3$seqinfo), throws_error("‘nrow\\(counts\\)’ must be identical to ‘nrow\\(pos\\)’"))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = endoapply(m2_s3$pos, function(x){colnames(x) <- rev(colnames(x)); x}), seqinfo = m2_s3$seqinfo), throws_error("‘m’ is set to 2 so colnames for all elements of ‘pos’ must be: ‘pos1’, ‘pos2’."))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'seqinfo'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos), throws_error("‘seqinfo’ missing."))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = as.data.frame(m2_s3$seqinfo)), throws_error("‘seqinfo’ must be a ‘Seqinfo’ object."))
  rm(m2_s3)
})

test_that("CoMeth parameter checking works: 'strand'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = Rle()), throws_error("‘strand’ must be an ‘RleList’."))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3 = Rle(factor('*'), length(m2_s3$seqnames[[3]])))
  strand <- endoapply(strand, function(x){levels(x) <- c('*', '+', '-', 'AAA'); x})
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = strand), throws_error("invalid strand levels in 'x': AAA"))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3 = Rle(factor('*'), length(m2_s3$seqnames[[3]])))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = endoapply(strand, function(x){x[-1]})), throws_error("Length of each ‘Rle’ element in ‘strand’ must equal the number of rows of its corresponding element in ‘pos’."))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3= Rle(factor('+'), length(m2_s3$seqnames[[1]])))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = strand), gives_warning("The supplied ‘strand’ argument contains ‘\\*’ as well as at least one of ‘\\+’ or ‘\\-’."))
  rm(m2_s3, strand)
})

test_that("CoMeth parameter checking works: 'colData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = data.frame()), throws_error("‘colData’, if supplied, must be a ‘DataFrame’."))
  colData <- DataFrame(a = c(1:3), row.names = paste0('sample', letters[1:3]))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = colData), throws_error("‘row.names\\(colData\\)’ must be identical to ‘sample_names’."))
  colData <- DataFrame(cancer = c(TRUE, TRUE, FALSE), row.names = paste0('sample', 1:3))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = colData), is_a("CoMeth"))
  rm(m2_s3, colData)
})

test_that("CoMeth parameter checking works: 'exptData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, exptData = list()), throws_error("‘exptData’, if supplied, must be a ‘SimpleList’."))
  exptData <- SimpleList(assay_type = "methylC-seq")
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, exptData = exptData), is_a("CoMeth"))
  rm(m2_s3, exptData)
})

test_that("CoMeth parameter checking works: '...'", {
  m2_s1 <-  make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  cgi <- rep(TRUE, 10)
  expect_that(CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo, cgi = cgi), is_a("CoMeth"))
  expect_that(mcols(CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo, cgi = cgi)), is_identical_to(DataFrame(cgi = cgi)))
  rm(m2_s1, cgi)
})

test_that("CoMeth parameter checking works: 'verbose'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, verbose = TRUE), is_a("CoMeth"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth validity
###
context("CoMeth validity")

test_that("CoMeth validity checking works on good object", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(validObject(m2_s3), is_true())
  rm(m2_s3)
})

test_that("CoMeth validity checking works: '.valid.CoMeth.counts'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  assay(m2_s3, 'MM') <- matrix(-1, ncol = ncol(m2_s3), nrow = nrow(m2_s3))
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: ‘counts’ cannot have negative entries."))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  names(assays(m2_s3)) <- rev(names(assays(m2_s3)))
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: assay names must be: assay names must be: ‘MM’, ‘MU’, ‘UM’, ‘UU’"))
  rm(m2_s3)
})

test_that("CoMeth validity checking works: '.valid.CoMeth.methylation_type'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  colData(m2_s3) <- colData(m2_s3)[, -1]
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: ‘colData’ of ‘CoMeth’ must contain column ‘methylation_type’ once and only once."))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  colData(m2_s3)$methylation_type <- rep('CpG', 3)
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: ‘methylation_type’ for each sample must be ‘CG’, ‘CHG’, ‘CHH’ or ‘CNN’, or some combination of these, e.g., ‘CG/CHG’.\nCombinations must sorted alphabetically and be separated by a forward slash \\(‘/’\\)."))
  rm(m2_s3)
})

test_that("CoMeth validity checking works: '.valid.CoMeth.rowData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  rowData(m2_s3) <- as(rowData(m2_s3), "GRanges")
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: ‘rowData\\(CoMeth\\)’ must be an ‘MTuples’ object."))
  rm(m2_s3)
})

test_that("CoMeth validity checking works: '.valid.CoMeth.noDuplicates'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  rowData(m2_s3)[2] <- rowData(m2_s3)[1]
  expect_that(validObject(m2_s3), throws_error("invalid class “CoMeth” object: ‘CoMeth’ object cannot contain duplicate m-tuples."))
  rm(m2_s3)
})

context("CoMeth methods")

test_that("'getPos' works", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  pos <- unlist(m2_s3$pos)
  pos <- unique(pos)
  row.names(pos) <- NULL
  pos <- as.matrix(pos)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(getPos(m2_s3), is_identical_to(pos))
  ## NOTE: The above test will fail if positions are not unique across samples
  ## The below test, which is a bit convoluted, should work.
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3$pos[[2]] <- rbind(m2_s3$pos[[1]][4:10, ], m2_s3$pos[[2]][1:3, ])
  pos <- cbind(as.vector(unlist(m2_s3$seqnames, use.names = FALSE)), as.matrix(unlist(m2_s3$pos)))
  pos <- unique(pos)
  row.names(pos) <- NULL
  pos <- pos[, -1]
  pos <- matrix(as.integer(pos), ncol = 2, dimnames = list(NULL, paste0('pos', 1:2)))
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(getPos(m2_s3)[order(getPos(m2_s3)[, 'pos1'], getPos(m2_s3)[, 'pos2']), ], is_identical_to(pos[order(pos[, 'pos1'], pos[, 'pos2']), ]))
})

test_that("'getM' works", {
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(getM(m1_s3), is_identical_to(1L))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(getM(m2_s3), is_identical_to(2L))
  m3_s3 <- make_test_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(sample_names = m3_s3$sample_names, methylation_type = m3_s3$methylation_type, counts = m3_s3$counts, seqnames = m3_s3$seqnames, pos = m3_s3$pos, seqinfo = m3_s3$seqinfo)
  expect_that(getM(m3_s3), is_identical_to(3L))
  m4_s3 <- make_test_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(sample_names = m4_s3$sample_names, methylation_type = m4_s3$methylation_type, counts = m4_s3$counts, seqnames = m4_s3$seqnames, pos = m4_s3$pos, seqinfo = m4_s3$seqinfo)
  expect_that(getM(m4_s3), is_identical_to(4L))
  rm(m1_s3, m2_s3, m3_s3, m4_s3)
})

test_that("'sampleNames' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', 1)))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', 1:2)))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', 1:3)))
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', 1:4)))
  rm(m2_s1, m2_s2, m2_s3, m2_s4)
})

test_that("'sampleNames<-' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  sampleNames(m2_s1) <- paste0('sample', letters[1])
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', letters[1])))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  sampleNames(m2_s2) <- paste0('sample', letters[1:2])
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', letters[1:2])))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  sampleNames(m2_s3) <- paste0('sample', letters[1:3])
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', letters[1:3])))
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  sampleNames(m2_s4) <- paste0('sample', letters[1:4])
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', letters[1:4])))
  rm(m2_s1, m2_s2, m2_s3, m2_s4)
})

test_that("'length' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(length(m2_s1), is_identical_to(1L))
  rm(m2_s1)
})

test_that("'nrow' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(nrow(m2_s1), is_identical_to(10L))
  rm(m2_s1)
})

test_that("'ncol' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(ncol(m2_s1), is_identical_to(1L))
  rm(m2_s1)
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(ncol(m2_s2), is_identical_to(2L))
  rm(m2_s2)
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(ncol(m2_s3), is_identical_to(3L))
  rm(m2_s3)
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  expect_that(ncol(m2_s4), is_identical_to(4L))
  rm(m2_s4)
})

test_that("'granges' is just an alias of 'rowData'", {
  ## granges(CoMeth) doesn't really make sense because the rowData is MTuples rather than GRanges.
  ## As of GenomicRanges_1.14.4, granges(SummarizedExperiment) is simply an alias of rowData(SummarizedExperiment), hence, by inheritance, granges(CoMeth) is just an alias of rowData(CoMeth).
  ## This test is in case this behaviour changes in future versions of GenomicRanges.
  m5_s4 <- make_test_data(m = 5L, n = 10L, s = 4L, sim_counts = TRUE)
  m5_s4 <- CoMeth(sample_names = m5_s4$sample_names, methylation_type = m5_s4$methylation_type, counts = m5_s4$counts, seqnames = m5_s4$seqnames, pos = m5_s4$pos, seqinfo = m5_s4$seqinfo)
  expect_that(granges(m5_s4), is_identical_to(rowData(m5_s4)))
  rm(m5_s4)
})

test_that("'cbind' works", {
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(cbind(m1_s2[, 1], m1_s2[, 2]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m1_s2, m1_s2), throws_error("Cannot ‘cbind’ ‘CoMeth’ objects with duplicate ‘sample_names’."))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(cbind(m2_s2[, 1], m2_s2[, 2]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m2_s2, m2_s2), throws_error("Cannot ‘cbind’ ‘CoMeth’ objects with duplicate ‘sample_names’."))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(cbind(m3_s2[, 1], m3_s2[, 2]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m3_s2, m3_s2), throws_error("Cannot ‘cbind’ ‘CoMeth’ objects with duplicate ‘sample_names’."))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(cbind(m4_s2[, 1], m4_s2[, 2]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m4_s2, m4_s2), throws_error("Cannot ‘cbind’ ‘CoMeth’ objects with duplicate ‘sample_names’."))
  expect_that(cbind(m2_s2[, 1], m4_s2[, 2]), throws_error("Cannot ‘cbind’ ‘CoMeth’ objects with the different sized m-tuples, that is, different ‘m’."))
  rm(m1_s2, m2_s2, m3_s2, m4_s2)
})

test_that("'rbind' works", {
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(rbind(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(rbind(m1_s2, m1_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(rbind(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m2_s2, m2_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(rbind(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m3_s2, m3_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(rbind(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(rbind(m4_s2, m4_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  rm(m1_s2, m2_s2, m3_s2, m4_s2)
})

test_that("'combine' works", {
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2))
  expect_that(combine(m1_s2, m1_s2), throws_error("‘combine’ failed when trying to ‘rbind’ the intermediate ‘CoMeth’ objects."))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2))
  expect_that(combine(m2_s2, m2_s2), throws_error("‘combine’ failed when trying to ‘rbind’ the intermediate ‘CoMeth’ objects."))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2))
  expect_that(combine(m3_s2, m3_s2), throws_error("‘combine’ failed when trying to ‘rbind’ the intermediate ‘CoMeth’ objects."))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2))
  expect_that(combine(m4_s2, m4_s2), throws_error("‘combine’ failed when trying to ‘rbind’ the intermediate ‘CoMeth’ objects."))
  rm(m1_s2, m2_s2, m3_s2, m4_s2)
})

test_that("'getCoverage' works", {
  m1_s1 <- make_test_data(m = 1L, n = 10L, s = 1L, sim_counts = TRUE)
  m1_s1 <- CoMeth(sample_names = m1_s1$sample_names, methylation_type = m1_s1$methylation_type, counts = m1_s1$counts, seqnames = m1_s1$seqnames, pos = m1_s1$pos, seqinfo = m1_s1$seqinfo)
  expect_that(getCoverage(m1_s1), is_identical_to(matrix(assay(m1_s1, 'M') + assay(m1_s1, 'U'), ncol = 1, dimnames = list(NULL, 'sample1'))))
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(getCoverage(m1_s2), is_identical_to(matrix(assay(m1_s2, 'M') + assay(m1_s2, 'U'), ncol = 2, dimnames = list(NULL, c('sample1', 'sample2')))))
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(getCoverage(m1_s3), is_identical_to(matrix(assay(m1_s3, 'M') + assay(m1_s3, 'U'), ncol = 3, dimnames = list(NULL, c('sample1', 'sample2', 'sample3')))))
  m1_s4 <- make_test_data(m = 1L, n = 10L, s = 4L, sim_counts = TRUE)
  m1_s4 <- CoMeth(sample_names = m1_s4$sample_names, methylation_type = m1_s4$methylation_type, counts = m1_s4$counts, seqnames = m1_s4$seqnames, pos = m1_s4$pos, seqinfo = m1_s4$seqinfo)
  expect_that(getCoverage(m1_s4), is_identical_to(matrix(assay(m1_s4, 'M') + assay(m1_s4, 'U'), ncol = 4, dimnames = list(NULL, c('sample1', 'sample2', 'sample3', 'sample4')))))
  rm(m1_s1, m1_s2, m1_s3, m1_s4)
})

## TODO: Test that duplicated works when m > = 3
test_that("'duplicated' works", {
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  start(m4_s2) <- 1
  end(m4_s2) <- 10
  m4_s2@rowData@extraPos <- matrix(c(rep(3, nrow(m4_s2)), rep(7, nrow(m4_s2))), ncol = 2)
  expect_that(any(duplicated(m4_s2)), is_true())
})