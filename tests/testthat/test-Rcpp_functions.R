context("'allRowsSortedCpp(NumericMatrix'")

test_that("'allRowsSortedCpp' works", {
  A <- matrix(1:100, byrow = TRUE, ncol = 10)
  expect_that(allRowsSortedCpp(A), is_true())
  expect_that(allRowsSortedCpp(A), is_identical_to(all(apply(X = A, FUN = function(x){!is.unsorted(x)}, MARGIN = 1))))
  B <- matrix(100:1, byrow = TRUE, ncol = 10)
  expect_that(allRowsSortedCpp(B), is_false())
  expect_that(allRowsSortedCpp(B), is_identical_to(all(apply(X = B, FUN = function(x){!is.unsorted(x)}, MARGIN = 1))))
  D <- matrix(c(1:3, 2:4), byrow = TRUE, ncol = 3)
  expect_that(allRowsSortedCpp(D), is_true())
  E <- matrix(c(1, 2, 3, 2, 3, 3), ncol = 3)
  expect_that(allRowsSortedCpp(E), is_false())
})

test_that("'rowDiffsCpp' works", {
  A <- matrix(1:100, byrow = TRUE, ncol = 10)
  expect_that(rowDiffsCpp(A), is_identical_to(matrix(1L, ncol = 9, nrow = 10)))
  B <- matrix(runif(100) + 100, ncol = 10)
  expect_that(rowDiffsCpp(B), throws_error("not compatible with requested type")) # Doesn't fail because floating point numbers are silently converted to integers
  D <- matrix(letters[1:20], ncol = 4)
  expect_that(rowDiffsCpp(D), throws_error("not compatible with requested type"))
})