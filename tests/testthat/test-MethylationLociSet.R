### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test MethylationLociSet constructor.
###

context("MethylationLociSet constructor")
test_that("Different methylation types", {
  mls <- make_MethylationLociSet_data(n = 10)
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = 'CG'), 
              is_a("MethylationLociSet"))
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = 'CHG'), 
              is_a("MethylationLociSet"))
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = 'CHH'), 
              is_a("MethylationLociSet"))
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = c('CG', 'CHG')), 
              is_a("MethylationLociSet"))
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = c('CNN')), 
              is_a("MethylationLociSet"))
})

test_that("Throws error if invalid methylation type", {
  mls <- make_MethylationLociSet_data(n = 10)
  expect_that(MethylationLociSet(seqnames = mls$seqnames, ranges = mls$ranges, 
                                 strand = mls$strand, seqinfo = mls$seqinfo, 
                                 methylation_type = c('CpG')), 
              throws_error(paste0("Invalid ", sQuote("methylation_type"))))
          })

## TODO: More tests of MethylationLociSet constructor and methods.