### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test betaCor function.
###

context("Parameter checks")
test_that("Nothing yet...", {
  mls <- make_MethylationLociSet_data(n = 1000)
  mls <- MethylationLociSet(mls$seqnames, mls$ranges, mls$strand, mls$seqinfo, 
                            'CG')
  cometh <- make_CoMeth1_data(mls = mls, f = 0.7, s = 4)
})