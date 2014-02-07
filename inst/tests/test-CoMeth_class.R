context("CoMeth constructor")

make_test_data <- function(m, n){
  test_data <- lapply(cbind(data.frame(chr = c(rep('chr1', 0.6 * n ), rep('chr2', 0.3 * n), rep('chrX', 0.1 * n)), stringsAsFactors = FALSE), as.data.frame(matrix(sort(sample(1:(n * m * 10), m * n, replace = FALSE)), ncol = m, byrow = T, dimnames = list(NULL, paste0('pos', 1:m)))), as.data.frame(matrix(rpois(2^m * n, 4), ncol = 2^m, dimnames = list(NULL, sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')})))))))), I)
  seqnames <- as.character(test_data[['chr']])
  pos <- matrix(unlist(test_data[grep('pos', names(test_data))]), ncol = m)
  colnames(pos) <- names(test_data[grep('pos', names(test_data))])
  counts <- matrix(unlist(test_data[grep('[MU]', names(test_data))]), ncol = 2 ^ m)
  colnames(counts) <- names(test_data[grep('[MU]', names(test_data))])
  return(list(seqnames = seqnames, pos = pos, counts = counts))
}

tsv <- make_test_data(2L, 100)

test_that("CoMeth works on good input",{
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), is_a("CoMeth"))
  })

test_that("CoMeth parameter checking works: 'seqnames'",{
  expect_that(CoMeth(seqnames = Rle(tsv$seqnames), tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), is_a("CoMeth"))
  expect_that(CoMeth(seqnames = factor(tsv$seqnames), tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), is_a("CoMeth"))
  expect_that(CoMeth(seqnames = rle(tsv$seqnames), tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), throws_error("Need 'seqnames'. Must be an Rle object, a character vector or a factor containing the sequence names"))
  expect_that(CoMeth(seqnames = NULL, tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), throws_error("Need 'seqnames'. Must be an Rle object, a character vector or a factor containing the sequence names"))
  expect_that(CoMeth(tsv$pos, tsv$counts, 2L, 'CG', 'tsv'), throws_error("Need 'seqnames'. Must be an Rle object, a character vector or a factor containing the sequence names"))
})

test_that("Cometh parameter checking works: 'pos'", {
  expect_that(CoMeth(tsv$seqnames, as.data.frame(tsv$pos), tsv$counts, 2L, 'CG', 'tsv'), throws_error("Need 'pos'. Must be a numeric matrix."))
  expect_that(CoMeth(tsv$seqnames, as.vector(tsv$pos), tsv$counts, 2L, 'CG', 'tsv'), throws_error("Need 'pos'. Must be a numeric matrix."))
  bad_pos <- tsv$pos
  dimnames(bad_pos) <- NULL
  expect_that(CoMeth(seqnames = tsv$seqnames, bad_pos, tsv$counts, 2L, 'CG', 'tsv'), throws_error("Column names of 'pos' must be: pos1, pos2, ..., posm"))
  colnames(bad_pos) <- paste0('pos_', 1:2)
  expect_that(CoMeth(seqnames = tsv$seqnames, bad_pos, tsv$counts, 2L, 'CG', 'tsv'), throws_error("Column names of 'pos' must be: pos1, pos2, ..., posm"))
})

test_that("Cometh parameter checking works: 'counts'", {
  expect_that(CoMeth(tsv$seqnames, tsv$pos, as.data.frame(tsv$counts), 2L, 'CG', 'tsv'), throws_error("Need 'counts'. Must be a numeric matrix."))
  expect_that(CoMeth(tsv$seqnames, tsv$pos, as.data.frame(tsv$counts), 2L, 'CG', 'tsv'), throws_error("Need 'counts'. Must be a numeric matrix."))
  bad_counts <- tsv$counts
  dimnames(bad_counts) <- NULL
  expect_that(CoMeth(seqnames = tsv$seqnames, tsv$pos, bad_counts, 2L, 'CG', 'tsv'), throws_error(paste0("Column names of 'counts' must be: ", paste0(make_m_tuple_names(2L), collapse = ', '))))
  colnames(bad_counts) <- rev(paste0(make_m_tuple_names(2L)))
  expect_that(CoMeth(seqnames = tsv$seqnames, tsv$pos, bad_counts, 2L, 'CG', 'tsv'), throws_error(paste0("Column names of 'counts' must be: ", paste0(make_m_tuple_names(2L), collapse = ', '))))
})

test_that("Cometh parameter checking works: 'm'", {
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, 2, 'CG', 'tsv'), throws_error("'m' must be specified. Must be a single int and must be greater than 0"))
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, c(2L, 2L), 'CG', 'tsv'), throws_error("'m' must be specified. Must be a single int and must be greater than 0"))
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, c(2L, 3L), 'CG', 'tsv'), throws_error("'m' must be specified. Must be a single int and must be greater than 0"))
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, 2.3, 'CG', 'tsv'), throws_error("'m' must be specified. Must be a single int and must be greater than 0"))
  expect_that(CoMeth(tsv$seqnames, tsv$pos, tsv$counts, -1L, 'CG', 'tsv'), throws_error("'m' must be specified. Must be a single int and must be greater than 0"))
  expect_that(CoMeth(seqnames = tsv$seqnames, tsv$pos, tsv$counts, 3L, 'CG', 'tsv'), throws_error("ncol\\(pos) should be equal to m")) # Note the annoying escaping of special characters
  expect_that(CoMeth(seqnames = tsv$seqnames, cbind(tsv$pos, tsv$pos), tsv$counts, 2L, 'CG', 'tsv'), throws_error("ncol\\(pos) should be equal to m")) # Note the annoying escaping of special characters
  expect_that(CoMeth(seqnames = tsv$seqnames, tsv$pos, cbind(tsv$counts, tsv$counts), 2L, 'CG', 'tsv'), throws_error("ncol\\(counts) should be equal to 2\\^m")) # Note the annoying escaping of special characters

})

test_that("Cometh parameter checking works: 'methylation_type'", {
  expect_that(CoMeth(seqnames = Rle(tsv$seqnames), tsv$pos, tsv$counts, 2L, c('CG', 'CHG'), 'tsv'), is_a("CoMeth"))
  expect_that(CoMeth(seqnames = Rle(tsv$seqnames), tsv$pos, tsv$counts, 2L, c('CGG', 'CHG'), 'tsv'), throws_error("'methylation_type' must be specified and must be a character vector. 'CG', 'CHG', 'CHH' and 'CNN', or a vector of some combination of these, are the only valid 'methylation_type'")) 

})

context("CoMeth accessors")
tsv <- make_test_data(3L, 100)
tsv_cometh <- CoMeth(tsv$seqnames, tsv$pos, tsv$counts, 3L, 'CG', 'tsv')
test_that("sampleNames works",{
  expect_that(sampleNames(tsv_cometh), is_identical_to('tsv'))
})

test_that("granges works",{
  expect_that(granges(tsv_cometh), is_a("GRanges"))
})

test_that("ncol works",{
  expect_that(ncol(tsv_cometh), equals(1))
})
test_that("getpos works",{
   expect_that(getPos(tsv_cometh), is_identical_to(tsv$pos))
})
test_that("getM works",{
  expect_that(getM(tsv_cometh), is_identical_to(3L))
})

test_that("getCoverage works",{
  expect_that(getCoverage(tsv_cometh), is_identical_to(rowSums(tsv$counts)))
})

context("CoMeth validity")
  tsv <- make_test_data(3L, 100)
  good_cometh <- CoMeth(tsv$seqnames, tsv$pos, tsv$counts, 3L, 'CG', 'tsv')

  test_that("validObject returns TRUE on good input", {
    expect_that(validObject(good_cometh), is_true())
  })

  test_that("validObject returns msg if 'counts' contains negative values", {
    bad_cometh <- good_cometh
    assay(bad_cometh, 'MMM') <- matrix(-1, ncol = 1, nrow = 100)
    expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'counts' has negative entries"))
  })

  test_that("validObject returns msg if 'pos' contains negative values", {
    bad_cometh <- good_cometh
    assay(bad_cometh, 'pos2') <- matrix(-1, ncol = 1, nrow = 100)
    expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'pos' has negative entries"))
  })

  test_that("validObject returns msg if 'pos' are not sorted",{
    bad_cometh <- good_cometh
    assay(bad_cometh, 'pos2') <- matrix(1, ncol = 1, nrow = 100)
    expect_that(validObject(bad_cometh), throws_error("Unsorted row in 'pos'"))
  })