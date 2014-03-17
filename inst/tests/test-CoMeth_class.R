#### Define test data ####
context("Define test data")

make_test_data <- function(m, n){
  test_data <- lapply(cbind(data.frame(chr = c(rep('chr1', 0.6 * n ), rep('chr2', 0.3 * n), rep('chrX', 0.1 * n)), stringsAsFactors = FALSE), as.data.frame(matrix(sort(sample(1:(n * m * 2), m * n, replace = FALSE)), ncol = m, byrow = T, dimnames = list(NULL, paste0('pos', 1:m)))), as.data.frame(matrix(rpois(2^m * n, 4), ncol = 2^m, dimnames = list(NULL, sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')})))))))), function(x){x})
  seqnames <- as.character(test_data[['chr']])
  pos <- lapply(test_data[grep('pos', names(test_data))], as.vector)
  #names(pos) <- names(test_data[grep('pos', names(test_data))])
  counts <- lapply(test_data[grep('[MU]', names(test_data))], as.vector)
  #names(counts) <- names(test_data[grep('[MU]', names(test_data))])
  return(list(seqnames = seqnames, pos = pos, counts = counts))
}

m <- 3L
a <- make_test_data(m, 200)
b <- make_test_data(m, 100)
d <- make_test_data(m, 1000)
e <- list(seqnames = list(a = a$seqnames, b = b$seqnames, d = d$seqnames), pos = list(a = a$pos, b = b$pos, d= d$pos), counts = list(a = a$counts, b = b$counts, d = d$counts))
seqnames <- e$seqnames
pos <- e$pos
counts <- e$counts
strand <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('+', 1000))
sample_names <- c('a', 'b', 'd')
methylation_type <- 'CG'
seq_info <- Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), seqlengths = c(249250621, 243199373, 155270560), genome = 'hg19')
good_cometh <- CoMeth(sample_names = c('a', 'b', 'd'), seqnames = e[['seqnames']], pos = e[['pos']], counts = e[['counts']], strand = e[['strand']], m = m, methylation_type = 'CG', seqinfo =  Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), seqlengths = rep(1000 * m * 2, 3), genome = 'hg19'))

#### Test CoMeth constructor ####
context("CoMeth constructor")

test_that("CoMeth works on good input: 1 sample",{
  expect_that(CoMeth(sample_names = 'a', seqnames = seqnames['a'], pos = pos['a'], counts = counts['a'], strand = strand['a'], m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  })

test_that("CoMeth works on good input: multiple sample",{
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
})

test_that("CoMeth parameter checking works: 'seqnames'",{
  expect_that(CoMeth(sample_names = sample_names, seqnames = RleList(seqnames), pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  expect_that(CoMeth(sample_names = sample_names, seqnames = lapply(seqnames, as.factor), pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  rle_seqnames <- lapply(seqnames, rle)
  expect_that(CoMeth(sample_names = sample_names, seqnames = rle_seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'seqnames'."))  
  null_seqnames <- NULL
  expect_that(CoMeth(sample_names = sample_names, seqnames = null_seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'seqnames'"))
  list_of_null_seqnames <- list(NULL)
  expect_that(CoMeth(sample_names = sample_names, seqnames = list_of_null_seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'seqnames'"))
  expect_that(CoMeth(sample_names = sample_names, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'seqnames'"))
})

test_that("CoMeth parameter checking works: 'pos'", {
  pos_missing_element <- pos[1:2]
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos_missing_element, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'pos'"))
  pos_as_df <- lapply(pos, function(x){lapply(x, as.data.frame)})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos_as_df, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'pos'"))
  no_name_pos <- lapply(pos, function(x){names(x) <- NULL; x})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = no_name_pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Names of each sub-list of 'pos' must be:"))  
  wrong_name_pos <- lapply(pos, function(x){names(x) <- paste0('pos_', 1:3); x})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = wrong_name_pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Names of each sub-list of 'pos' must be:"))  
})

test_that("CoMeth parameter checking works: 'counts'", {
  rev_name_counts <- lapply(counts, function(x){names(x) <- rev(paste0(make_m_tuple_names(3L))); x})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = rev_name_counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  no_name_counts <- lapply(counts, function(x){names(x) <- NULL; x})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = no_name_counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error(" Names of each sub-list of 'counts' must be:"))
  bad_name_counts <- lapply(counts, function(x){names(x) <- rev(paste0(make_m_tuple_names(2L))); x})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = bad_name_counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error(" Names of each sub-list of 'counts' must be:"))
})

test_that("CoMeth parameter checking works: 'm'", {
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = 3, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("'m' must be specified and must be a single, positive integer."))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = list(3L, 3L, 3L), methylation_type = methylation_type, seqinfo = seqinfo), throws_error("'m' must be specified and must be a single, positive integer."))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = c(3L, 3L, 3L), methylation_type = methylation_type, seqinfo = seqinfo), throws_error("'m' must be specified and must be a single, positive integer."))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = 2.3, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("'m' must be specified and must be a single, positive integer."))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = -1L, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("'m' must be specified and must be a single, positive integer."))
  wrong_m <- 2L
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = wrong_m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("Need 'pos'."))
})

test_that("CoMeth parameter checking works: 'methylation_type'", {
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = c('CG', 'CHG'), seqinfo = seqinfo), is_a("CoMeth"))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = c('CGG', 'CHG'), seqinfo = seqinfo), throws_error("'methylation_type' must be specified."))
})

test_that("CoMeth parameter checking works: 'strand'", {
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = lapply(strand, as.vector), m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = lapply(strand, as.factor), m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  null_strand <- NULL
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = null_strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), is_a("CoMeth"))
  list_of_null_strand <- vector('list', m)
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = list_of_null_strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("If 'strand' is not NULL then it must be a list of character vectors, a list of factors or an RleList object containing the strand of each m-tuple."))
  rle_strand <- lapply(strand, function(x){rle(as.vector(x))})
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = rle_strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo), throws_error("If 'strand' is not NULL then it must be a list of character vectors, a list of factors or an RleList object containing the strand of each m-tuple."))  
})

test_that("CoMeth parameter checking works: 'seqinfo'", {
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type), throws_error("'seqinfo' must be specified."))  
})

test_that("CoMeth parameter checking works: 'sort_cometh'", {
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo, sort_cometh = TRUE), is_a("CoMeth"))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo, sort_cometh = FALSE), is_a("CoMeth"))
  expect_that(CoMeth(sample_names = sample_names, seqnames = seqnames, pos = pos, counts = counts, strand = strand, m = m, methylation_type = methylation_type, seqinfo = seqinfo, sort_cometh = t), throws_error("'sort_cometh' must"))
})


#### UP TO HERE ####

#### Test CoMeth methods and functions ####
context("CoMeth accessors")

test_that("sampleNames works",{
  expect_that(sampleNames(good_cometh), is_identical_to(c('a', 'b', 'd')))
})

test_that("granges works",{
  expect_that(granges(good_cometh), is_a("GRanges"))
})

test_that("ncol works",{
  expect_that(ncol(good_cometh), equals(m))
})
test_that("getPos works",{
  expect_true(FALSE) # Need to figure out how to test this for a non-trivial example
})
test_that("getM works",{
  expect_that(getM(good_cometh), is_identical_to(m))
})

test_that("getCoverage works",{
  expect_true(FALSE) # Need to figure out how to test this for a non-trivial example
  #expect_that(getCoverage(tsv_cometh), is_identical_to(rowSums(tsv$counts)))
})

context("CoMeth validity")

test_that("validObject returns TRUE on good input", {
  expect_true(validObject(good_cometh))
  })

test_that("validObject returns msg if 'counts' contains negative values", {
  bad_cometh <- good_cometh
  assay(bad_cometh, 'MMM') <- matrix(-1, ncol = ncol(good_cometh), nrow = nrow(good_cometh))
  expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'counts' has negative entries"))
  })

test_that("validObject returns msg if 'extra_pos' contains negative values", {
  bad_cometh <- good_cometh
  bad_cometh@extra_pos <- list('pos2' = rep(-1, nrow(good_cometh)))
  expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'pos' has negative entries"))
  })

test_that("validObject returns msg if 'extra_pos' contains the wrong number of positions", {
  bad_cometh <- good_cometh
  bad_cometh@extra_pos <- list('pos2' = rep(-1, nrow(good_cometh)), 'pos2' = rep(-1, nrow(good_cometh)))
  expect_that(validObject(bad_good_cometh), throws_error("invalid class “CoMeth” object: 2: 'extra_pos' must be a list of length"))
  })

test_that("validObject returns msg if 'm' = 2 and extra_pos slot doesn't contain an empty list", {
  bad_test_data <- make_test_data(2, 100)
  bad_cometh <- CoMeth(sample_names = 'bad', seqnames = list('bad' = bad_test_data$seqnames), pos = list('bad' = bad_test_data$pos), counts = list('bad' = bad_test_data$counts), m = 2L, strand = NULL, methylation_type = methylation_type, seqinfo = seqinfo)
  bad_cometh@extra_pos <- list(1:100)
  expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'extra_pos' must be an empty list if 'm' = 2."))
  })
    
test_that("validObject returns msg if element of 'extra_pos' is of the wrong length", {
    bad_cometh <- good_cometh
    bad_cometh@extra_pos <- list('pos2' = rep(-1, 100))
    expect_true(FALSE)
    #expect_that(validObject(bad_cometh), throws_error()) # Throws error but not a very helpful one
  })

  test_that("validObject returns msg if 'pos' are not sorted",{
    #bad_cometh <- good_cometh
    #assay(bad_cometh, 'pos2') <- matrix(1, ncol = 1, nrow = 100)
    expect_true(FALSE)
    #expect_that(validObject(bad_cometh), throws_error("Unsort_comethed row in 'pos'"))
  })