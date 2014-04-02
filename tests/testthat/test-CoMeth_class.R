#### Define test data ####
context("Define test data")

make_test_data <- function(m, n){
  test_data <- lapply(cbind(data.frame(chr = c(rep('chr1', 0.6 * n ), rep('chr2', 0.3 * n), rep('chrX', 0.1 * n)), stringsAsFactors = FALSE), as.data.frame(matrix(sort(sample(1:(n * m * 2), m * n, replace = FALSE)), ncol = m, byrow = T, dimnames = list(NULL, paste0('pos', 1:m)))), as.data.frame(matrix(rpois(2^m * n, 4), ncol = 2^m, dimnames = list(NULL, sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')})))))))), function(x){x})
  pos <- DataFrame(seqnames = Rle(as.character(test_data[['chr']])), lapply(test_data[grep('pos', names(test_data))], as.vector))
  counts <- DataFrame(lapply(test_data[grep('[MU]', names(test_data))], as.vector))
  return(list(pos = pos, counts = counts))
}

set.seed(666)
m <- 3L
a <- make_test_data(m, 200)
b <- make_test_data(m, 100)
d <- make_test_data(m, 1000)
e <- list(pos = DataFrameList(a = a$pos, b = b$pos, d= d$pos), counts = DataFrameList(a = a$counts, b = b$counts, d = d$counts))
pos <- DataFrameList(a = a$pos, b = b$pos, d= d$pos)
counts <-  DataFrameList(a = a$counts, b = b$counts, d = d$counts)
strand <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('*', 1000))
sample_names <- c('a', 'b', 'd')
methylation_type <- CharacterList(a = 'CG', b = 'CG', d = 'CG')
seq_info <- Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), seqlengths = c(249250621, 243199373, 155270560), genome = 'hg19')
good_cometh <- CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info)

a_cometh <- CoMeth(m = m, sample_names = sample_names[1], pos = pos[1], counts = counts[1], strand = strand[1], methylation_type = methylation_type[1], seqinfo = seq_info)
a1_cometh <- CoMeth(m = m, sample_names = sample_names[1], pos = DataFrameList(lapply(pos[1], function(x){x[1:100,]})), counts = DataFrameList(lapply(counts[1], function(x){x[1:100,]})), strand = RleList(lapply(strand[1], function(x){x[1:100]})), methylation_type = methylation_type[1], seqinfo = seq_info)
a2_cometh <- CoMeth(m = m, sample_names = sample_names[1], pos = DataFrameList(lapply(pos[1], function(x){x[101:200,]})), counts = DataFrameList(lapply(counts[1], function(x){x[101:200,]})), strand = RleList(lapply(strand[1], function(x){x[101:200]})), methylation_type = methylation_type[1], seqinfo = seq_info)

b_cometh <- CoMeth(m = m, sample_names = sample_names[2], pos = pos[2], counts = counts[2], strand = strand[2], methylation_type = methylation_type[2], seqinfo = seq_info)
d_cometh <- CoMeth(m = m, sample_names = sample_names[3], pos = pos[3], counts = counts[3], strand = strand[3], methylation_type = methylation_type[3], seqinfo = seq_info)

m1 <- make_test_data(1, 50)
m1_cometh <- CoMeth(m = 1L, sample_names = 'm1', pos = DataFrameList(m1 = m1$pos), counts = DataFrameList(m1 = m1$count), strand = RleList(m1 = Rle('*', 50)), methylation_type = CharacterList(m1 = 'CHH'), seqinfo = seq_info)

#### Test CoMeth constructor ####
context("CoMeth constructor")

test_that("CoMeth works on good input: 1 sample",{
  expect_that(CoMeth(m = m, sample_names = 'a', pos = pos['a'], counts = counts['a'], strand = strand['a'], methylation_type = methylation_type['a'], seqinfo = seq_info), is_a("CoMeth")) # Same as a1_cometh
  })

test_that("CoMeth works on good input: multiple sample",{
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), is_a("CoMeth")) # Same as good_cometh
})

test_that("CoMeth parameter checking works: 'seqnames'",{
  bad_seqnames_in_pos <- DataFrameList(lapply(pos, function(x){x$seqnames <- 1; x}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = bad_seqnames_in_pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("Class of columns in each DataFrame element of 'pos' must be:"))
  no_seqnames_in_pos <- DataFrameList(lapply(pos, function(x){x[, -1]}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = no_seqnames_in_pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' is set to 3 so colnames for all elemnts of 'pos' must be:"))
})

test_that("CoMeth parameter checking works: 'pos'", {
  pos_missing_element <- DataFrameList(lapply(pos, function(x){x[, - 3]}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos_missing_element, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("colnames for all elemnts of 'pos' must be:"))
  pos_short_of_rows <- DataFrameList(lapply(pos, function(x){x[-1, ]}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos_short_of_rows, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("nrow\\(counts\\) must be equal to nrow\\(pos\\)."))
  pos_not_integer <- DataFrameList(lapply(pos, function(x){x[, 2] <- x[, 2] + 0.5; x}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos_not_integer, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("Class of columns in each DataFrame element of 'pos' must be:"))
  expect_that(CoMeth(m = m, sample_names = sample_names, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'pos' missing"))
  pos_bad_names <- DataFrameList(lapply(pos, function(x){colnames(x) <- c('seqnames', paste0('pos_', 1:3)); x}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos_bad_names, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' is set to 3 so colnames for all elemnts of 'pos' must be:"))
})

test_that("CoMeth parameter checking works: 'counts'", {
  counts_missing_element <- DataFrameList(lapply(counts, function(x){x[, - 3]}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts_missing_element, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("colnames for all elements of 'counts' must be:"))
  counts_short_of_rows <- DataFrameList(lapply(counts, function(x){x[-1, ]}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts_short_of_rows, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("nrow\\(counts\\) must be equal to nrow\\(pos\\)."))
  counts_not_integer <- DataFrameList(lapply(counts, function(x){x[, 2] <- x[, 2] + 0.5; x}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts_not_integer, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("Class of columns in each DataFrame element of 'counts' must be:"))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'counts' missing"))
  counts_bad_names <- DataFrameList(lapply(counts, function(x){names(x) <- rev(make_m_tuple_names(3L)); x}))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts_bad_names, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("colnames for all elements of 'counts' must be:"))
})

test_that("CoMeth parameter checking works: 'm'", {
  m_not_integer <- 3
  expect_that(CoMeth(m = m_not_integer, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' must be specified and must be a single, positive integer."))
  m_list <- list(3L, 3L, 3L)
  expect_that(CoMeth(m = m_list, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' must be specified and must be a single, positive integer."))
  m_decimal <- 2.3
  expect_that(CoMeth(m = m_decimal, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' must be specified and must be a single, positive integer."))
  m_negative <- -1L
  expect_that(CoMeth(m = m_negative, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' must be specified and must be a single, positive integer."))
  wrong_m <- 2L
  expect_that(CoMeth(m = wrong_m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info), throws_error("'m' is set to 2 so colnames for all elemnts of 'pos' must be:"))
})

test_that("CoMeth parameter checking works: 'methylation_type'", {
  multiple_methylation_types <- CharacterList(a = 'CG', b = 'CG', d = c('CHG', 'CG'))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = multiple_methylation_types, seqinfo = seq_info), is_a("CoMeth"))
  bad_methylation_types <- CharacterList(a = 'CA', b = 'CG', d = 'CG')
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = bad_methylation_types, seqinfo = seq_info), throws_error("'methylation_type' for each sample must be 'CG', 'CHG', 'CHH' and 'CNN' or a vector of some combination of these, e.g. "))  
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, , seqinfo = seq_info), throws_error("'methylation_type' missing"))
})

test_that("CoMeth parameter checking works: 'strand'", {
  strand_with_plus_and_star <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('+', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_plus_and_star, methylation_type = methylation_type, seqinfo = seq_info), gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'."))
  strand_with_minus_and_star <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('-', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_minus_and_star, methylation_type = methylation_type, seqinfo = seq_info), gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'."))
  strand_with_plus_and_minus_and_star <- RleList(a = Rle('*', 200), b = Rle('-', 100), d = Rle('+', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_plus_and_minus_and_star, methylation_type = methylation_type, seqinfo = seq_info), gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'."))
  strand_with_plus_and_minus <- RleList(a = Rle('+', 200), b = Rle('-', 100), d = Rle('+', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_plus_and_minus, methylation_type = methylation_type, seqinfo = seq_info), not(gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'.")))
  strand_with_star <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('*', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_star, methylation_type = methylation_type, seqinfo = seq_info), not(gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'.")))
  strand_with_plus <- RleList(a = Rle('+', 200), b = Rle('+', 100), d = Rle('+', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_plus, methylation_type = methylation_type, seqinfo = seq_info), not(gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'.")))  
  strand_with_minus <- RleList(a = Rle('-', 200), b = Rle('-', 100), d = Rle('-', 1000))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_with_minus, methylation_type = methylation_type, seqinfo = seq_info), not(gives_warning("'strand' contains '\\*' as well as at least one of '\\+' or '-'.")))  
  strand_wrong_length <- RleList(a = Rle('-', 200), b = Rle('-', 100), d = Rle('-', 100))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand_wrong_length, methylation_type = methylation_type, seqinfo = seq_info), throws_error("Length of each element in 'strand' must equal the number of rows of its corresponding element in 'pos'."))    
})

test_that("CoMeth parameter checking works: 'seqinfo'", {
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type), throws_error("'seqinfo' missing"))
  
})

test_that("CoMeth parameter checking works: 'sort_cometh'", {
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info, sort_cometh = TRUE), not(throws_error("'sort_cometh' must be TRUE or FALSE")))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info, sort_cometh = FALSE), not(throws_error("'sort_cometh' must be TRUE or FALSE")))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info, sort_cometh = 'TRUE'), throws_error("'sort_cometh' must be TRUE or FALSE"))
  expect_that(CoMeth(m = m, sample_names = sample_names, pos = pos, counts = counts, strand = strand, methylation_type = methylation_type, seqinfo = seq_info, sort_cometh = 1), throws_error("'sort_cometh' must be TRUE or FALSE"))
})

#### Test CoMeth validity ####
context("CoMeth validity")
## TODO: Add more tests of validity checking

test_that("validObject returns TRUE on good input", {
  expect_true(validObject(good_cometh))
})

test_that("validObject returns msg if 'counts' contains negative values", {
  bad_cometh <- good_cometh
  assay(bad_cometh, 'MMM') <- matrix(-1, ncol = ncol(good_cometh), nrow = nrow(good_cometh))
  expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'counts' has negative entries"))
})

test_that("validObject returns msg if 'pos' contains negative values", {
  bad_cometh <- good_cometh
  values(rowData(bad_cometh)) <- DataFrame('pos2' = rep(-1, nrow(good_cometh)))
  expect_that(validObject(bad_cometh), throws_error("invalid class “CoMeth” object: 'pos' has negative entries"))
})

test_that("validObject returns msg if 'pos' are not sorted",{
  #bad_cometh <- good_cometh
  #assay(bad_cometh, 'pos2') <- matrix(1, ncol = 1, nrow = 100)
  expect_true(FALSE)
  #expect_that(validObject(bad_cometh), throws_error("Unsort_comethed row in 'pos'"))
})

test_that("validObject returns msg if 'm' or 'methylation_type' incorrectly specified", {
  bad_cometh <- good_cometh
  colData(bad_cometh) <- colData(bad_cometh)[, 1, drop = FALSE]
  expect_that(validObject(bad_cometh), throws_error("colData of 'CoMeth' must contain columns 'm' and 'methylation_type' once each, and once each only."))
  bad_cometh <- good_cometh
  colData(bad_cometh) <- DataFrame(m = rep(4L, 3), methylation_type = rep('CG', 3))
  expect_that(validObject(bad_cometh), throws_error("1: length\\(assays\\) does not equal"))
})

#### Functions that work on CoMeth objects ####
context("Functions defined for CoMeth")

test_that("getPos works",{
  expect_true(FALSE) # Need to figure out how to test this for a non-trivial example
})

test_that("getM works",{
  expect_that(getM(good_cometh), is_identical_to(m))
})

#### CoMeth methods ####
context("Test CoMeth methods")

test_that("length", {
  expect_that(length(good_cometh), equals(1294))
})

test_that("sampleNames works", {
  expect_that(sampleNames(good_cometh), is_identical_to(c('a', 'b', 'd')))
  good_cometh_with_new_sample_names <- good_cometh
})

test_that("sampleNames<- works", {
  sampleNames(good_cometh_with_new_sample_names) <- c('s1', 's2', 's3')
  expect_that(validObject(good_cometh_with_new_sample_names), is_true())
  expect_that(sampleNames(good_cometh_with_new_sample_names), equals(c('s1', 's2', 's3')))
})

test_that("granges works", {
  expect_that(granges(good_cometh), is_a("GRanges"))
})

test_that("ncol/NCOL works", {
  expect_that(ncol(good_cometh), equals(m))
  expect_that(NCOL(good_cometh), equals(m))
})

test_that("nrow/NROW works", {
  expect_that(nrow(good_cometh), equals(1294))
  expect_that(nrow(good_cometh), equals(1294))
})

test_that("colData works", {
  expect_that(colData(good_cometh), is_identical_to(DataFrame(m = rep(3L, 3), methylation_type = rep('CG', 3), row.names = c('a', 'b', 'd'))))
})

test_that("colData<- works", {
  good_cometh_new_colData <- good_cometh
  colData(good_cometh_new_colData) <- DataFrame(m = rep(3L, 3), methylation_type = rep('CG/CHG', 3), row.names = c('a1', 'b1', 'd1'))
  expect_that(validObject(good_cometh_new_colData), is_true())
})

test_that("cbind works", {
  cbind_cometh <- cbind(a_cometh, b_cometh, d_cometh)
  ## Can't just check objects are identical but can check they are equal and that all slots are identical
  expect_that(cbind_cometh, equals(good_cometh))
  expect_that(rowData(cbind_cometh), is_identical_to(rowData(good_cometh)))
  expect_that(colData(cbind_cometh), is_identical_to(colData(good_cometh)))
  expect_that(exptData(cbind_cometh), is_identical_to(exptData(good_cometh)))
  expect_that(assays(cbind_cometh), is_identical_to(assays(good_cometh)))  
  ## Cases where cbind should break
  b_cometh_bad_seqinfo <- b_cometh
  genome(b_cometh_bad_seqinfo) <- c(rep('mm10', 3))
  expect_that(cbind(a_cometh, b_cometh_bad_seqinfo), throws_error("Can only cbind 'CoMeth' objects with identical 'seqinfo'."))
  expect_that(cbind(good_cometh, good_cometh), throws_error("Cannot cbind 'CoMeth' objects containing duplicate 'sample_names'"))  
  expect_that(cbind(a_cometh, m1_cometh), throws_error("Cannot cbind 'CoMeth' objects with the different 'm'."))
  a_cometh_CHG <- a_cometh
  colData(a_cometh_CHG)$methylation_type <- 'CHG'
  expect_that(cbind(a_cometh_CHG, b_cometh), gives_warning("Combining 'CoMeth' objects with different 'methylation_type'."))
})

test_that("rbind works", {
  ## Can't just check objects are identical but can check they are equal and that all slots are identical
  rbind_cometh <- rbind(a1_cometh, a2_cometh)
  expect_that(rbind_cometh, equals(a_cometh))
  expect_that(rowData(rbind_cometh), is_identical_to(rowData(a_cometh)))
  expect_that(colData(rbind_cometh), is_identical_to(colData(a_cometh)))
  expect_that(exptData(rbind_cometh), is_identical_to(exptData(a_cometh)))
  expect_that(assays(rbind_cometh), is_identical_to(assays(a_cometh)))
  ## Cases where bbind should break
  a2_cometh_bad_seqinfo <- a2_cometh
  genome(a2_cometh_bad_seqinfo) <- rep('mm10', 3)
  expect_that(rbind(a1_cometh, a2_cometh_bad_seqinfo), throws_error("Can only rbind 'CoMeth' objects with identical 'seqinfo'"))
  a2_cometh_rename <- a2_cometh
  sampleNames(a2_cometh_rename) <- 'b'
  expect_that(rbind(a1_cometh, a2_cometh_rename), throws_error("Cannot rbind CoMeth objects with different sampleNames"))
  a_cometh_rename <- a_cometh
  sampleNames(a_cometh_rename) <- 'm1'
  expect_that(rbind(m1_cometh, a_cometh_rename), throws_error("Cannot rbind 'CoMeth' objects with the different 'm'."))
  a2_cometh_CHG <- a2_cometh
  colData(a2_cometh_CHG)$methylation_type <- 'CHG'
  expect_that(rbind(a1_cometh, a2_cometh_CHG), gives_warning("Combining 'CoMeth' objects with different 'methylation_type'."))
})

test_that("compare works", {
  expect_that(zero_range(compare(good_cometh, good_cometh)), is_true())
  expect_that(compare(m1_cometh, rev(m1_cometh)), is_identical_to(compare(rowData(m1_cometh), rev(rowData(m1_cometh)))))
  expect_that(compare(good_cometh, m1_cometh), throws_error("Cannot compare CoMeth object with different 'm'"))
  good_cometh_diff_seqinfo <- good_cometh
  genome(good_cometh_diff_seqinfo) <- 'mm10'
  expect_that(compare(good_cometh, good_cometh_diff_seqinfo), throws_error())
  ## Test ==, <=, !=, >=, <, >
  ## Test with 3-tuples
  expect_that(compare(good_cometh, good_cometh) == 0, is_identical_to(good_cometh == good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) == 0, is_identical_to(good_cometh == rev(good_cometh)))
  expect_that(compare(good_cometh, good_cometh) <= 0, is_identical_to(good_cometh <= good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) <= 0, is_identical_to(good_cometh <= rev(good_cometh)))
  expect_that(compare(good_cometh, good_cometh) != 0, is_identical_to(good_cometh != good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) != 0, is_identical_to(good_cometh != rev(good_cometh)))
  expect_that(compare(good_cometh, good_cometh) >= 0, is_identical_to(good_cometh >= good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) >= 0, is_identical_to(good_cometh >= rev(good_cometh)))
  expect_that(compare(good_cometh, good_cometh) < 0, is_identical_to(good_cometh < good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) < 0, is_identical_to(good_cometh < rev(good_cometh)))
  expect_that(compare(good_cometh, good_cometh) > 0, is_identical_to(good_cometh > good_cometh))
  expect_that(compare(good_cometh, rev(good_cometh)) > 0, is_identical_to(good_cometh > rev(good_cometh)))
  ## Test with 1-tuples
  expect_that(compare(m1_cometh, m1_cometh) == 0, is_identical_to(m1_cometh == m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) == 0, is_identical_to(m1_cometh == rev(m1_cometh)))
  expect_that(compare(m1_cometh, m1_cometh) <= 0, is_identical_to(m1_cometh <= m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) <= 0, is_identical_to(m1_cometh <= rev(m1_cometh)))
  expect_that(compare(m1_cometh, m1_cometh) != 0, is_identical_to(m1_cometh != m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) != 0, is_identical_to(m1_cometh != rev(m1_cometh)))
  expect_that(compare(m1_cometh, m1_cometh) >= 0, is_identical_to(m1_cometh >= m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) >= 0, is_identical_to(m1_cometh >= rev(m1_cometh)))
  expect_that(compare(m1_cometh, m1_cometh) < 0, is_identical_to(m1_cometh < m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) < 0, is_identical_to(m1_cometh < rev(m1_cometh)))
  expect_that(compare(m1_cometh, m1_cometh) > 0, is_identical_to(m1_cometh > m1_cometh))
  expect_that(compare(m1_cometh, rev(m1_cometh)) > 0, is_identical_to(m1_cometh > rev(m1_cometh)))
})

test_that("dim works", {
  expect_that(dim(good_cometh), is_identical_to(c(nrow(good_cometh), ncol(good_cometh))))
})

test_that("dimnames works", {
  expect_that(dimnames(good_cometh), is_identical_to(list(names(rowData(good_cometh)), rownames(colData(x)))))
})

test_that("end works", {
  expect_that(end(good_cometh), is_identical_to(end(rowData(good_cometh))))
})

test_that("duplicated works", {
  expect_that(any(duplicated(good_cometh)), is_false())
  expect_that(duplicated(good_cometh), not(is_identical_to(duplicated(rowData(good_cometh)))))
})

test_that("mcols works", {
  expect_that(mcols(good_cometh), is_identical_to(DataFrame(pos2 = getPos(good_cometh)[, 3])))
})

test_that("mcols<- works", {
  ## TODO: Write test
  expect_that(FALSE, is_true())
})



#### OLD CODE BELOW THIS LINE ####


test_that("getCoverage works",{
  expect_true(FALSE) # Need to figure out how to test this for a non-trivial example
  #expect_that(getCoverage(tsv_cometh), is_identical_to(rowSums(tsv$counts)))
})

