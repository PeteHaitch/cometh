### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Define test data.
###

pos <- DataFrame(seqnames = rep('chr1', 100), strand = rep('*', 100), pos1 = seq(from = 1000L, to = 1L, by = -10L), pos2 = seq(from = 1003L, to = 4L, by = -10L))

tuple <- Tuple(pos, seqinfo = Seqinfo(seqnames = 'chr1', seqlengths = 10000, genome = 'hg19'))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test Tuple constructor.
###
context("Tuple constructor")
expect_true(FALSE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test Tuple validity.
###
context("Tuple validity")
expect_true(FALSE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test Tuple methods.
###
context("Tuple methods")
expect_true(FALSE)