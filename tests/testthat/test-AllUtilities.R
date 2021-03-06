context("'.make_m_tuple_names'")
test_that("'.make_m_tuple_names' works",{
  expect_that(.make_m_tuple_names(1L), is_identical_to(c('M', 'U')))
  expect_that(.make_m_tuple_names(2L), is_identical_to(c('MM', 'MU', 'UM', 'UU')))
  expect_that(.make_m_tuple_names(3L), is_identical_to(c('MMM', 'MMU', 'MUM', 'MUU', 'UMM', 'UMU', 'UUM', 'UUU')))
  expect_that(.make_m_tuple_names('A'), throws_error("'m' must be an int. 'm' must be greater than 0."))
  expect_that(.make_m_tuple_names(2), throws_error("'m' must be an int. 'm' must be greater than 0."))
  expect_that(.make_m_tuple_names(-1L), throws_error("'m' must be an int. 'm' must be greater than 0."))
})

context("'.zero_range'")
test_that("'.zero_range' works", {
  expect_that(.zero_range(rep(10L, 10)), is_true())
})

context("'.LOR'")
test_that(".LOR' works", {
  x <- DataFrame('MM' = 0, 'UU' = 0, 'MU' = 1, 'UM' = 0)
  expect_that(.LOR(x), is_identical_to((log2(0.5) + log2(0.5) - log2(1.5) - log2(0.5))))
  x <- DataFrame('MM' = 0, 'UU' = 0, 'MU' = 0, 'UM' = 0)
  expect_that(.LOR(x), is_identical_to(0))
})