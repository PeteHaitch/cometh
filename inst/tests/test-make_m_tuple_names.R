context("utils.r")
test_that("make_m_tuple_names works as expected",{
  expect_that(make_m_tuple_names(1L), is_identical_to(c('M', 'U')))
  expect_that(make_m_tuple_names(2L), is_identical_to(c('MM', 'MU', 'UM', 'UU')))
  expect_that(make_m_tuple_names(3L), is_identical_to(c('MMM', 'MMU', 'MUM', 'MUU', 'UMM', 'UMU', 'UUM', 'UUU')))
  expect_that(make_m_tuple_names('A'), throws_error("'m' must be an int. 'm' must be greater than 0."))
  expect_that(make_m_tuple_names(2), throws_error("'m' must be an int. 'm' must be greater than 0."))
  expect_that(make_m_tuple_names(-1L), throws_error("'m' must be an int. 'm' must be greater than 0."))
})