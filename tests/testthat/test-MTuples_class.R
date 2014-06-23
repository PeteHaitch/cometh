### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test MTuples constructor.
###
context("MTuples constructor")

test_that("MTuples constructor returns an MTuples object when m = 0", {
  expect_that(MTuples(seqnames = Rle(), pos = matrix()), is_a("MTuples"))
})

test_that("MTuples constructor returns an MTuples object when m = 1", {
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(1:10, ncol = 1)), is_a("MTuples"))
})

test_that("MTuples constructor returns an MTuples object when m = 2", {
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(1:20, ncol = 2)), is_a("MTuples"))
})

test_that("MTuples constructor returns an MTuples object when m >= 3", {
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(1:30, ncol = 3)), is_a("MTuples"))
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(1:40, ncol = 4)), is_a("MTuples"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test MTuples validity.
###
context("MTuples validity")

test_that("Ensure m-tuple positions are checked to be sorted", {
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(c(10:20, 15:5), ncol = 2)), throws_error(paste0("positions within each 2-tuple must be sorted in strictly increasing order, i.e. ", sQuote('pos1'), " < ", sQuote('pos2'), '.')))
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(c(1:10, 21:30, 11:20), ncol = 3)), throws_error("positions within each 3-tuple must be sorted in strictly increasing order"))
  expect_that(MTuples('chr1', pos = matrix(c(1, 1), ncol = 2)), throws_error("positions within each 2-tuple must be sorted in strictly increasing order")) # Test for GitHub issue #8 (https://github.com/PeteHaitch/cometh/issues/8)
})

test_that("Ensure m-tuples positions are checked to be all positive", {
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(c(-1:10), ncol = 1)), throws_error(paste0("Some values in ", sQuote('pos'), " are negative.")))
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(c(1:15, -8), ncol = 2)), throws_error(paste0("Some values in ", sQuote('pos'), " are negative.")))
  expect_that(MTuples(seqnames = 'chr1', pos = matrix(c(1:15, -8, 17:30), ncol = 3)), throws_error(paste0("Some values in ", sQuote('pos'), " are negative.")))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test MTuples methods.
###
context("getPos method for MTuples objects")

test_that("getPos works when m = 1", {
  expect_that(getPos(MTuples(seqnames = 'chr1', pos = matrix(1:10, ncol = 1))), is_equivalent_to(matrix(1:10, ncol = 1))) # The latter argument will not have column names, hence the need to use is_equivalent_to rather than is_equal_to
})

test_that("getPos works when m = 2", {
  expect_that(getPos(MTuples(seqnames = 'chr1', pos = matrix(1:20, ncol = 2))), is_equivalent_to(matrix(1:20, ncol = 2))) # The latter argument will not have column names, hence the need to use is_equivalent_to rather than is_equal_to
})

test_that("getPos works when m >= 3", {
  expect_that(getPos(MTuples(seqnames = 'chr1', pos = matrix(1:30, ncol = 3))), is_equivalent_to(matrix(1:30, ncol = 3)))
  expect_that(getPos(MTuples(seqnames = 'chr1', pos = matrix(1:40, ncol = 4))), is_equivalent_to(matrix(1:40, ncol = 4)))
})

context("getM method for MTuples objects")

test_that("getM works when m = 1", {
  expect_that(getM(MTuples(seqnames = 'chr1', pos = matrix(1:10, ncol = 1))), is_identical_to(1L))
})

test_that("getM works when m = 2", {
  expect_that(getM(MTuples(seqnames = 'chr1', pos = matrix(1:20, ncol = 2))), is_identical_to(2L))  
})

test_that("getM works when m >= 3", {
  expect_that(getM(MTuples(seqnames = 'chr1', pos = matrix(1:30, ncol = 3))), is_identical_to(3L))
  expect_that(getM(MTuples(seqnames = 'chr1', pos = matrix(1:40, ncol = 4))), is_identical_to(4L))
})

context("c method for MTuples objects")

test_that("c works when there are no common m-tuples and m = 1", {
  m1_x <- MTuples("chr1", pos = matrix(1:10, ncol = 1))
  m1_y <- MTuples('chr1', pos = matrix(11:20, ncol = 1))
  expect_that(c(m1_x, m1_y), is_identical_to(MTuples("chr1", pos = matrix(1:20, ncol = 1))))
})

test_that("c works when all m-tuples are common m-tuples and m = 1", {
  m1_x <- MTuples("chr1", pos = matrix(1:10, ncol = 1))
  expect_that(c(m1_x, m1_x), is_identical_to( MTuples("chr1", pos = matrix(c(1:10, 1:10), ncol = 1))))
})

test_that("c works when there are some common m-tuples and m = 1", {
  m1_x <- MTuples("chr1", pos = matrix(1:10, ncol = 1))
  m1_y <- MTuples('chr1', pos = matrix(6:15, ncol = 1))
  expect_that(c(m1_x, m1_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 6:15), ncol = 1))))
})

test_that("c works with multiple MTuples object and m = 1", {
  m1_x <- MTuples("chr1", pos = matrix(1:10, ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m1_z <- MTuples("chr1", pos = matrix(21:30, ncol = 1))
  expect_that(c(m1_x, m1_y, m1_z), is_identical_to(MTuples("chr1", pos = matrix(1:30, ncol = 1))))
})

test_that("c works when there are no common m-tuples and m = 2", {
  m2_x <- MTuples("chr1", pos = matrix(1:20, ncol = 2))
  m2_y <- MTuples('chr1', pos = matrix(21:40, ncol = 2))
  expect_that(c(m2_x, m2_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 21:30, 11:20, 31:40), ncol = 2))))
})

test_that("c works when all m-tuples are common m-tuples and m = 2", {
  m2_x <- MTuples("chr1", pos = matrix(1:20, ncol = 2))
  expect_that(c(m2_x, m2_x), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 1:10, 11:20, 11:20), ncol = 2))))
})

test_that("c works when there are some common m-tuples and m = 2", {
  m2_x <- MTuples("chr1", pos = matrix(1:20, ncol = 2))
  m2_y <- MTuples('chr1', pos = matrix(11:30, ncol = 2))
  expect_that(c(m2_x, m2_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:20, 11:30), ncol = 2))))
})

test_that("c works with multiple MTuples object and m = 2", {
  m2_x <- MTuples("chr1", pos = matrix(1:20, ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m2_z <- MTuples("chr1", pos = matrix(21:40, ncol = 2))
  expect_that(c(m2_x, m2_y, m2_z), is_identical_to(MTuples('chr1', pos = matrix(c(1:30, 11:40), ncol = 2))))
})

test_that("c works when there are no common m-tuples and m >= 3", {
  m3_x <- MTuples("chr1", pos = matrix(1:30, ncol = 3))
  m3_y <- MTuples('chr1', pos = matrix(31:60, ncol = 3))
  expect_that(c(m3_x, m3_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 31:40, 11:20, 41:50, 21:30, 51:60), ncol = 3))))
  rm(m3_x, m3_y)
  m4_x <- MTuples("chr1", pos = matrix(1:40, ncol = 4))
  m4_y <- MTuples('chr1', pos = matrix(41:80, ncol = 4))
  expect_that(c(m4_x, m4_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 41:50, 11:20, 51:60, 21:30, 61:70, 31:40, 71:80), ncol = 4))))
})

test_that("c works when all m-tuples are common m-tuples and m >= 3", {
  m3_x <- MTuples("chr1", pos = matrix(1:30, ncol = 3))
  expect_that(c(m3_x, m3_x), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 1:10, 11:20, 11:20, 21:30, 21:30), ncol = 3))))
  rm(m3_x)
  m4_x <- MTuples("chr1", pos = matrix(1:40, ncol = 4))
  expect_that(c(m4_x, m4_x), is_identical_to(MTuples("chr1", pos = matrix(c(1:10, 1:10, 11:20, 11:20, 21:30, 21:30, 31:40, 31:40), ncol = 4))))
})

test_that("c works when there are some common m-tuples and m >= 3", {
  m3_x <- MTuples("chr1", pos = matrix(1:30, ncol = 3))
  m3_y <- MTuples('chr1', pos = matrix(11:40, ncol = 3))
  expect_that(c(m3_x, m3_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:20, 11:30, 21:40), ncol = 3))))
  rm(m3_x, m3_y)
  m4_x <- MTuples("chr1", pos = matrix(1:40, ncol = 4))
  m4_y <- MTuples('chr1', pos = matrix(11:50, ncol = 4))
  expect_that(c(m4_x, m4_y), is_identical_to(MTuples("chr1", pos = matrix(c(1:20, 11:30, 21:40, 31:50), ncol = 4))))
  rm(m4_x, m4_y)
})

test_that("c works with multiple MTuples object and m >= 3", {
  m3_x <- MTuples("chr1", pos = matrix(1:30, ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  m3_z <- MTuples("chr1", pos = matrix(21:50, ncol = 3))
  expect_that(c(m3_x, m3_y, m3_z), is_identical_to(MTuples("chr1", pos = matrix(c(1:30, 11:40, 21:50), ncol = 3))))
  rm(m3_x, m3_y, m3_z)
  m4_x <- MTuples("chr1", pos = matrix(1:40, ncol = 4))
  m4_y <- MTuples("chr1", pos = matrix(11:50, ncol = 4))
  m4_z <- MTuples("chr1", pos = matrix(21:60, ncol = 4))
  expect_that(c(m4_x, m4_y, m4_z), is_identical_to(MTuples("chr1", pos = matrix(c(1:30, 11:40, 21:50, 31:60), ncol = 4))))
})

context("IPD method for MTuples objects")

test_that("IPD works when m = 1", {
  m1_x <- MTuples("chr1", pos = matrix(1:10, ncol = 1))
  expect_that(getIPD(m1_x), throws_error(paste0("It does not make sense to compute IPD when ", sQuote('m'), " = 1.")))
})

test_that("IPD works when m = 2", {
  pos <- matrix(sort(sample(10000, size = 20)), ncol = 2)
  m2_x <- MTuples("chr1", pos = pos)
  expect_that(getIPD(m2_x), is_identical_to(matrix(pos[, 2L] - pos[, 1L], ncol = 1)))
})

test_that("IPD works when m >= 3", {
  pos <- matrix(sort(sample(10000, size = 30)), ncol = 3)
  m3_x <- MTuples("chr1", pos = pos)
  expect_that(getIPD(m3_x), is_identical_to(matrix(c(pos[, 2L] - pos[, 1L], pos[, 3L] - pos[, 2L]), ncol = 2)))
  rm(pos, m3_x)
  pos <- matrix(sort(sample(10000, size = 40)), ncol = 4)
  m4_x <- MTuples("chr1", pos = pos)
  expect_that(getIPD(m4_x), is_identical_to(matrix(c(pos[, 2L] - pos[, 1L], pos[, 3L] - pos[, 2L], pos[, 4L] - pos[, 3L]), ncol = 3)))
})

context("'compare'-based methods for MTuples")

test_that("'compare', works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(compare(m1_x, m1_y), is_identical_to(c(-1L, rep(0L, 9))))
  expect_that(compare(m2_x, m2_y), is_identical_to(c(-2L, -2L, rep(0L, 8))))
  expect_that(compare(m3_x, m3_y), is_identical_to(c(-2L, -2L, -2L, rep(0L, 7))))
})

test_that("'<=' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x <= m1_y, is_identical_to(rep(TRUE, 10)))
  expect_that(m2_x <= m2_y, is_identical_to(rep(TRUE, 10)))
  expect_that(m3_x <= m3_y, is_identical_to(rep(TRUE, 10)))
})

test_that("'==' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x == m1_y, is_identical_to(c(FALSE, rep(TRUE, 9))))
  expect_that(m2_x == m2_y, is_identical_to(c(rep(FALSE, 2), rep(TRUE, 8))))
  expect_that(m3_x == m3_y, is_identical_to(c(rep(FALSE, 3), rep(TRUE, 7))))
})

test_that("'>=' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x >= m1_y, is_identical_to(c(FALSE, rep(TRUE, 9))))
  expect_that(m2_x >= m2_y, is_identical_to(c(rep(FALSE, 2), rep(TRUE, 8))))
  expect_that(m3_x >= m3_y, is_identical_to(c(rep(FALSE, 3), rep(TRUE, 7))))
})

test_that("'<' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x < m1_y, is_identical_to(c(TRUE, rep(FALSE, 9))))
  expect_that(m2_x < m2_y, is_identical_to(c(rep(TRUE, 2), rep(FALSE, 8))))
  expect_that(m3_x < m3_y, is_identical_to(c(rep(TRUE, 3), rep(FALSE, 7))))
})

test_that("'>' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x > m1_y, is_identical_to(rep(FALSE, 10)))
  expect_that(m2_x > m2_y, is_identical_to(rep(FALSE, 10)))
  expect_that(m3_x > m3_y, is_identical_to(rep(FALSE, 10))) 
})

test_that("'!=' works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(m1_x != m1_y, is_identical_to(c(TRUE, rep(FALSE, 9))))
  expect_that(m2_x != m2_y, is_identical_to(c(rep(TRUE, 2), rep(FALSE, 8))))
  expect_that(m3_x != m3_y, is_identical_to(c(rep(TRUE, 3), rep(FALSE, 7))))
})

test_that("'m' check works", {
  m1_x <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(11:20, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(c(9:10, 13:30), ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:30, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(c(9:11, 14:40), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(11:40, ncol = 3))
  expect_that(compare(m1_x, m2_x), throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(compare(m1_x, m3_x), throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(m1_x <= m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(m1_x < m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(m1_x == m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(m1_x != m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
  expect_that(m1_x >= m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))  
  expect_that(m1_x > m2_x, throws_error(paste0("Cannot ", sQuote('compare'), " ",  sQuote('MTuples'), " objects with different ", sQuote('m'), ".")))
})

context("'findOverlaps'-based methods work")
test_that("'findOverlaps' with 'type = equal'", {
  m1_x <- MTuples("chr1", pos = matrix(1:6, ncol = 1))
  m1_y <- MTuples("chr1", pos = matrix(5:10, ncol = 1))
  m2_x <- MTuples("chr1", pos = matrix(1:12, byrow = TRUE, ncol = 2))
  m2_y <- MTuples("chr1", pos = matrix(11:20, byrow = TRUE, ncol = 2))
  m3_x <- MTuples("chr1", pos = matrix(1:18, byrow = TRUE, ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(13:24, byrow = TRUE, ncol = 3))
  
  expect_that(queryHits(findOverlaps(m1_x, m1_y, type = 'equal', select = 'all')), is_identical_to(5:6))
  expect_that(subjectHits(findOverlaps(m1_x, m1_y, type = 'equal', select = 'all')), is_identical_to(1:2))
  expect_that(queryHits(findOverlaps(m2_x, m2_y, type = 'equal', select = 'all')), is_identical_to(6L))
  expect_that(subjectHits(findOverlaps(m2_x, m2_y, type = 'equal', select = 'all')), is_identical_to(1L))
  expect_that(queryHits(findOverlaps(m3_x, m3_y, type = 'equal', select = 'all')), is_identical_to(5:6))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'equal', select = 'all')), is_identical_to(1:2))
})

test_that("'findOverlaps' with 'type = any, start, end or within'", {
  m3_x <- MTuples("chr1", pos = matrix(c(10, 20, 30), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(c(20, 30, 40), ncol = 3))
  expect_that(queryHits(findOverlaps(m3_x, m3_y, type = 'any')), is_identical_to(1L))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'any')), is_identical_to(1L))
  expect_that(queryHits(findOverlaps(m3_x, m3_y, type = 'start')), is_identical_to(integer(0)))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'start')), is_identical_to(integer(0)))  
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'end')), is_identical_to(integer(0)))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'within')), is_identical_to(integer(0)))  
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'equal')), is_identical_to(integer(0)))
  rm(m3_x, m3_y)
  m3_x <- MTuples("chr1", pos = matrix(c(10, 20, 30), ncol = 3))
  m3_y <- MTuples("chr1", pos = matrix(c(5, 30, 40), ncol = 3))
  expect_that(queryHits(findOverlaps(m3_x, m3_y, type = 'any')), is_identical_to(1L))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'any')), is_identical_to(1L))
  expect_that(queryHits(findOverlaps(m3_x, m3_y, type = 'start')), is_identical_to(integer(0)))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'end')), is_identical_to(integer(0)))
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'within')), is_identical_to(1L))  
  expect_that(subjectHits(findOverlaps(m3_x, m3_y, type = 'equal')), is_identical_to(integer(0)))
})

test_that("'findOverlaps' fails with different 'm'", {
  m1 <- MTuples("chr1", pos = matrix(1:6, ncol = 1))
  m2 <- MTuples("chr1", pos = matrix(1:12, byrow = TRUE, ncol = 2))
  expect_that(findOverlaps(m1, m2, type = 'any'), throws_error(paste0("Cannot ", sQuote("findOverlaps"), " between ", sQuote("MTuples"), " and ", sQuote("MTuples"), " if they have different ", sQuote("m"), ".")))
  expect_that(findOverlaps(m1, m2, type = 'equal'), throws_error(paste0("Cannot ", sQuote("findOverlaps"), " between ", sQuote("MTuples"), " and ", sQuote("MTuples"), " if they have different ", sQuote("m"), ".")))
})

context("mcols works")
test_that("'mcols' works", {
  m1 <- MTuples("chr1", pos = matrix(c(10, 12:20), ncol = 1))
  mcols(m1)$CGI <- c(rep(TRUE, 5), rep(FALSE, 5))
  expect_that(mcols(m1), is_identical_to(DataFrame(CGI = c(rep(TRUE, 5), rep(FALSE, 5)))))
})

context("duplicated works")
test_that("'fromLast = FALSE' works", {
  m1 <- MTuples("chr1", pos = matrix(rep(11:13, each = 3), ncol = 1), strand = rep(c('+', '-', '*'), times = 3))
  m2 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3)), ncol = 2), strand = rep(c('+', '-', '*'), times = 3))
  m3 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3), rep(17:19, each = 3)), ncol = 3), strand = rep(c('+', '-', '*'), times = 3))
  expect_that(duplicated(c(m1, m1[1, ])), is_identical_to(c(rep(FALSE, 9), TRUE)))
  expect_that(duplicated(c(m1, m1)), is_identical_to(c(rep(FALSE, 9), rep(TRUE, 9))))
  expect_that(duplicated(c(m2, m2[1, ])), is_identical_to(c(rep(FALSE, 9), TRUE)))
  expect_that(duplicated(c(m2, m2)), is_identical_to(c(rep(FALSE, 9), rep(TRUE, 9))))
  expect_that(duplicated(c(m3, m3[1, ])), is_identical_to(c(rep(FALSE, 9), TRUE)))
  expect_that(duplicated(c(m3, m3)), is_identical_to(c(rep(FALSE, 9), rep(TRUE, 9))))
  expect_that(duplicated(c(MTuples('chr1', pos = matrix(c(11, 14, 18, 11, 15, 17, 12, 14, 17), byrow = TRUE, ncol = 3)), m3)), is_identical_to(rep(FALSE, 12)))
})

test_that("'fromLast = TRUE' works", {
  m1 <- MTuples("chr1", pos = matrix(rep(11:13, each = 3), ncol = 1), strand = rep(c('+', '-', '*'), times = 3))
  m2 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3)), ncol = 2), strand = rep(c('+', '-', '*'), times = 3))
  m3 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3), rep(17:19, each = 3)), ncol = 3), strand = rep(c('+', '-', '*'), times = 3))
  expect_that(duplicated(c(m1, m1[1, ]), fromLast = TRUE), is_identical_to(c(TRUE, rep(FALSE, 9))))
  expect_that(duplicated(c(m1, m1), fromLast = TRUE), is_identical_to(c(rep(TRUE, 9), rep(FALSE, 9))))
  expect_that(duplicated(c(m2, m2[1, ]), fromLast = TRUE), is_identical_to(c(TRUE, rep(FALSE, 9))))
  expect_that(duplicated(c(m2, m2), fromLast = TRUE), is_identical_to(c(rep(TRUE, 9), rep(FALSE, 9))))
  expect_that(duplicated(c(m3, m3[1, ]), fromLast = TRUE), is_identical_to(c(TRUE, rep(FALSE, 9))))
  expect_that(duplicated(c(m3, m3), fromLast = TRUE), is_identical_to(c(rep(TRUE, 9), rep(FALSE, 9))))
  expect_that(duplicated(c(MTuples('chr1', pos = matrix(c(11, 14, 18, 11, 15, 17, 12, 14, 17), byrow = TRUE, ncol = 3)), m3), fromLast = TRUE), is_identical_to(rep(FALSE, 12)))
})

test_that("'anyDuplicated' works", {
  m1 <- MTuples("chr1", pos = matrix(rep(11:13, each = 3), ncol = 1), strand = rep(c('+', '-', '*'), times = 3))
  m2 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3)), ncol = 2), strand = rep(c('+', '-', '*'), times = 3))
  m3 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3), rep(17:19, each = 3)), ncol = 3), strand = rep(c('+', '-', '*'), times = 3))
  expect_that(anyDuplicated(m1), is_identical_to(0L))
  expect_that(anyDuplicated(m2), is_identical_to(0L))
  expect_that(anyDuplicated(m3), is_identical_to(0L))
  expect_that(anyDuplicated(c(m1, m1[1, ])), is_identical_to(10L))
  expect_that(anyDuplicated(c(m1, m1)), is_identical_to(10L))
  expect_that(anyDuplicated(c(m2, m2[1, ])), is_identical_to(10L))
  expect_that(anyDuplicated(c(m2, m2)), is_identical_to(10L))
  expect_that(anyDuplicated(c(m3, m3[1, ])), is_identical_to(10L))
  expect_that(anyDuplicated(c(m3, m3)), is_identical_to(10L))  
  expect_that(anyDuplicated(c(MTuples('chr1', pos = matrix(c(11, 14, 18, 11, 15, 17, 12, 14, 17), byrow = TRUE, ncol = 3)), m3)), is_identical_to(0L))
})

context("unique works")
test_that("'unique' works", {
  m1 <- MTuples("chr1", pos = matrix(rep(11:13, each = 3), ncol = 1), strand = rep(c('+', '-', '*'), times = 3))
  m2 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3)), ncol = 2), strand = rep(c('+', '-', '*'), times = 3))
  m3 <- MTuples("chr1", pos = matrix(c(rep(11:13, each = 3), rep(14:16, each = 3), rep(17:19, each = 3)), ncol = 3), strand = rep(c('+', '-', '*'), times = 3))
  expect_that(unique(c(m1, m1[1, ])), is_identical_to(m1))
  expect_that(unique(c(m1, m1)), is_identical_to(m1))
  expect_that(unique(c(m2, m2[1, ])), is_identical_to(m2))
  expect_that(unique(c(m2, m2)), is_identical_to(m2))
  expect_that(unique(c(m3, m3[1, ])), is_identical_to(m3))
  expect_that(unique(c(m3, m3)), is_identical_to(m3))
  expect_that(unique(c(MTuples('chr1', pos = matrix(c(11, 14, 18, 11, 15, 17, 12, 14, 17), byrow = TRUE, ncol = 3)), m3)), is_identical_to(c(MTuples('chr1', pos = matrix(c(11, 14, 18, 11, 15, 17, 12, 14, 17), byrow = TRUE, ncol = 3)), m3)))
})