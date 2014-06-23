## TODO: It's possible to have two SummarizedExperiments (well, CoMeth) objects that have identical slots, all.equal returns TRUE but identical returns FALSE. 
## E.g. a <- make_MTuples_or_CoMeth_data(3L, 10L, 3L, sim_counts = TRUE), identical(combine(a[, 1], a[, 2], a[, 3]), a)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth constructor.
###
context('CoMeth constructor')

test_that("CoMeth constructor returns an CoMeth object when m = 1 for a single sample", {
  m1_s1 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 1L, 
                                       sim_counts = TRUE)
  m1_s1 <- CoMeth(assays = m1_s1$assays, rowData = m1_s1$rowData, 
                  colData = m1_s1$colData)
  expect_that(m1_s1, is_a("CoMeth"))
  expect_that(m1_s1, is_a("CoMeth1"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 1 for multiple samples", {
  m1_s3 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(assays = m1_s3$assays, rowData = m1_s3$rowData, 
                  colData = m1_s3$colData)
  expect_that(m1_s3, is_a("CoMeth"))
  expect_that(m1_s3, is_a("CoMeth1"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for a single sample", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  expect_that(m2_s1, is_a("CoMeth"))
  expect_that(m2_s1, is_a("CoMeth2"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for multiple samples", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  expect_that(m2_s3, is_a("CoMeth"))
  expect_that(m2_s3, is_a("CoMeth2"))
})

test_that("CoMeth constructor returns an CoMeth object when m > 2 for a single sample", {
  m3_s1 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 1L, sim_counts = TRUE)
  m3_s1 <- CoMeth(assays = m3_s1$assays, rowData = m3_s1$rowData, 
                  colData = m3_s1$colData)
  expect_that(m3_s1, is_a("CoMeth"))
  expect_that(m3_s1, is_a("CoMeth3Plus"))
  m4_s1 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 1L, sim_counts = TRUE)
  m4_s1 <- CoMeth(assays = m4_s1$assays, rowData = m4_s1$rowData, 
                  colData = m4_s1$colData)
  expect_that(m4_s1, is_a("CoMeth"))
  expect_that(m4_s1, is_a("CoMeth3Plus"))
})

test_that("CoMeth constructor returns an CoMeth object when m > 2 for multiple samples", {
  m3_s3 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(assays = m3_s3$assays, rowData = m3_s3$rowData, 
                  colData = m3_s3$colData)
  expect_that(m3_s3, is_a("CoMeth"))
  expect_that(m3_s3, is_a("CoMeth3Plus"))
  m4_s3 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(assays = m4_s3$assays, rowData = m4_s3$rowData, 
                  colData = m4_s3$colData)
  expect_that(m4_s3, is_a("CoMeth"))
  expect_that(m4_s3, is_a("CoMeth3Plus"))
})

## TODO: CoMeth parameter checking works

test_that("CoMeth parameter checking works: 'colData'", {
  ## TODO: 
})

test_that("CoMeth parameter checking works: 'exptData'", {
  ## TODO
})

test_that("CoMeth parameter checking works: '...'", {
  ## TODO
})

test_that("CoMeth parameter checking works: 'verbose'", {
  ## TODO
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth validity
###
context("CoMeth validity")

test_that("CoMeth validity checking works on good object", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  expect_that(validObject(m2_s3), is_true())
})

test_that("CoMeth validity checking works: '.valid.CoMeth.counts'", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  assay(m2_s3, 'MM') <- matrix(-1, ncol = ncol(m2_s3), nrow = nrow(m2_s3))
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('counts'), " cannot have negative entries.")))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  names(assays(m2_s3)) <- paste0((names(assays(m2_s3))), 'a')
  expect_that(validObject(m2_s3), throws_error(paste0("assay names must include: ", paste0(sQuote(.make_m_tuple_names(2L)), collapse = ', '), ".")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.methylation_type'", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  colData(m2_s3) <- colData(m2_s3)[, -1]
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('colData'), " of ", sQuote('CoMeth'), " must contain column ", sQuote('methylation_type'), " once and only once.")))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  colData(m2_s3)$methylation_type <- rep('CpG', 3)
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('methylation_type'), " for each sample must be ", sQuote('CG'), ", ", sQuote('CHG'), ", ", sQuote('CHH'), " or ", sQuote('CNN'), ", or some combination of these, e.g., ", sQuote("CG/CHG"), ".")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.rowData'", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  rowData(m2_s3) <- as(rowData(m2_s3), "GRanges")
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('rowData\\(CoMeth\\)'), " must be an ", sQuote('MTuples'), " object.")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.noDuplicates'", {
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  rowData(m2_s3)[2] <- rowData(m2_s3)[1]
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('CoMeth'), " object cannot contain duplicate m-tuples.")))
})

context("CoMeth methods")

test_that("'getPos' works", {
  ## TODO
  
#   m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
#   pos <- unlist(m2_s3$pos)
#   pos <- unique(pos)
#   row.names(pos) <- NULL
#   pos <- as.matrix(pos)
#   m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
#                   colData = m2_s3$colData)
#   ## NOTE: "expect_that(getPos(m2_s3), is_identical_to(pos))" will fail if 
#   ## positions are not unique across samples
#   ## The below test, which is a bit convoluted, should always work.
#   m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
#   m2_s3$pos[[2]] <- rbind(m2_s3$pos[[1]][4:10, ], m2_s3$pos[[2]][1:3, ])
#   pos <- cbind(as.vector(unlist(m2_s3$seqnames, use.names = FALSE)), as.matrix(unlist(m2_s3$pos)))
#   pos <- unique(pos)
#   row.names(pos) <- NULL
#   pos <- pos[, -1]
#   pos <- matrix(as.integer(pos), ncol = 2, dimnames = list(NULL, paste0('pos', 1:2)))
#   m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
#   expect_that(getPos(m2_s3)[order(getPos(m2_s3)[, 'pos1'], getPos(m2_s3)[, 'pos2']), ], is_identical_to(pos[order(pos[, 'pos1'], pos[, 'pos2']), ]))
})

test_that("'getM' works", {
  m1_s3 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(assays = m1_s3$assays, rowData = m1_s3$rowData, 
                  colData = m1_s3$colData)
  expect_that(getM(m1_s3), is_identical_to(1L))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  expect_that(getM(m2_s3), is_identical_to(2L))
  m3_s3 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(assays = m3_s3$assays, rowData = m3_s3$rowData, 
                  colData = m3_s3$colData)
  expect_that(getM(m3_s3), is_identical_to(3L))
  m4_s3 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(assays = m4_s3$assays, rowData = m4_s3$rowData, 
                  colData = m4_s3$colData)
  expect_that(getM(m4_s3), is_identical_to(4L))
})

test_that("'sampleNames' works", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', 1)))
  m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(assays = m2_s2$assays, rowData = m2_s2$rowData, 
                  colData = m2_s2$colData)
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', 1:2)))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', 1:3)))
  m2_s4 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(assays = m2_s4$assays, rowData = m2_s4$rowData, 
                  colData = m2_s4$colData)
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', 1:4)))
})

test_that("'sampleNames<-' works", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  sampleNames(m2_s1) <- paste0('sample', letters[1])
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', letters[1])))
  m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(assays = m2_s2$assays, rowData = m2_s2$rowData, 
                  colData = m2_s2$colData)
  sampleNames(m2_s2) <- paste0('sample', letters[1:2])
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', letters[1:2])))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  sampleNames(m2_s3) <- paste0('sample', letters[1:3])
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', letters[1:3])))
  m2_s4 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(assays = m2_s4$assays, rowData = m2_s4$rowData, 
                  colData = m2_s4$colData)
  sampleNames(m2_s4) <- paste0('sample', letters[1:4])
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', letters[1:4])))
})

test_that("'length' works", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  expect_that(length(m2_s1), is_identical_to(1L))
})

test_that("'nrow' works", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  expect_that(nrow(m2_s1), is_identical_to(10L))
})

test_that("'ncol' works", {
  m2_s1 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(assays = m2_s1$assays, rowData = m2_s1$rowData, 
                  colData = m2_s1$colData)
  expect_that(ncol(m2_s1), is_identical_to(1L))
  m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(assays = m2_s2$assays, rowData = m2_s2$rowData, 
                  colData = m2_s2$colData)
  expect_that(ncol(m2_s2), is_identical_to(2L))
  m2_s3 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(assays = m2_s3$assays, rowData = m2_s3$rowData, 
                  colData = m2_s3$colData)
  expect_that(ncol(m2_s3), is_identical_to(3L))
  m2_s4 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(assays = m2_s4$assays, rowData = m2_s4$rowData, 
                  colData = m2_s4$colData)
  expect_that(ncol(m2_s4), is_identical_to(4L))
})

test_that("'granges' is just an alias of 'rowData'", {
  ## granges(CoMeth) doesn't really make sense because the rowData is MTuples rather than GRanges.
  ## As of GenomicRanges_1.14.4, granges(SummarizedExperiment) is simply an alias of rowData(SummarizedExperiment), hence, by inheritance, granges(CoMeth) is just an alias of rowData(CoMeth).
  ## This test is in case this behaviour changes in future versions of GenomicRanges.
  m5_s4 <- make_MTuples_or_CoMeth_data(m = 5L, n = 10L, s = 4L, sim_counts = TRUE)
  m5_s4 <- CoMeth(assays = m5_s4$assays, rowData = m5_s4$rowData, 
                  colData = m5_s4$colData)
  expect_that(granges(m5_s4), is_identical_to(rowData(m5_s4)))
})

test_that("'cbind' works", {
  ## TODO
#   # TODO: Test cbind fails if objects have different seqinfo
#   # TODO: Test cbind for other failure modes
#   m1_s2 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
#   m1_s2 <- CoMeth(assays = m1_s2$assays, rowData = m1_s2$rowData, 
#                   colData = m1_s2$colData)
#   expect_that(cbind(m1_s2[, 1], m1_s2[, 2]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(cbind(m1_s2, m1_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
#   m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
#   m2_s2 <- CoMeth(assays = m2_s2$assays, rowData = m2_s2$rowData, 
#                   colData = m2_s2$colData)
#   expect_that(cbind(m2_s2[, 1], m2_s2[, 2]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(cbind(m2_s2, m2_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
#   m3_s2 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
#   m3_s2 <- CoMeth(assays = m3_s2$assays, rowData = m3_s2$rowData, 
#                   colData = m3_s2$colData)
#   expect_that(cbind(m3_s2[, 1], m3_s2[, 2]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(cbind(m3_s2, m3_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
#   m4_s2 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
#   m4_s2 <- CoMeth(assays = m4_s2$assays, rowData = m4_s2$rowData, 
#                   colData = m4_s2$colData)
#   expect_that(cbind(m4_s2[, 1], m4_s2[, 2]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(cbind(m4_s2, m4_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
#   expect_that(cbind(m2_s2[, 1], m4_s2[, 2]), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples, that is, different ", sQuote('m'), ".")))
})

test_that("'rbind' works", {
  # TODO
#   # TODO: Test rbind fails if objects have different seqinfo
#   # TODO: Test rbind for other failure modes
#   m1_s2 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
#   m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
#   expect_that(rbind(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(rbind(m1_s2, m1_s2), throws_error(paste0("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")))
#   m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
#   m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
#   expect_that(rbind(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m2_s2, m2_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
#   m3_s2 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
#   m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
#   expect_that(rbind(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m3_s2, m3_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
#   m4_s2 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
#   m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
#   expect_that(rbind(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(rbind(m4_s2, m4_s2), throws_error(paste0("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")))
})

test_that("'combine' works", {
  # TODO
#   m1_s2 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
#   m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
#   expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2))
#   expect_that(combine(m1_s2, m1_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
#   m2_s2 <- make_MTuples_or_CoMeth_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
#   m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
#   expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2))
#   expect_that(combine(m2_s2, m2_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
#   m3_s2 <- make_MTuples_or_CoMeth_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
#   m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
#   expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2))
#   expect_that(combine(m3_s2, m3_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
#   m4_s2 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
#   m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
#   expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
#   expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2))
#   expect_that(combine(m4_s2, m4_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
})

test_that("'getCoverage' works", {
  m1_s1 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 1L, sim_counts = TRUE)
  m1_s1 <- CoMeth(assays = m1_s1$assays, rowData = m1_s1$rowData, 
                  colData = m1_s1$colData)
  expect_that(getCoverage(m1_s1), is_identical_to(matrix(assay(m1_s1, 'M') + assay(m1_s1, 'U'), ncol = 1, dimnames = list(NULL, 'sample1'))))
  m1_s2 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(assays = m1_s2$assays, rowData = m1_s2$rowData, 
                  colData = m1_s2$colData)
  expect_that(getCoverage(m1_s2), is_identical_to(matrix(assay(m1_s2, 'M') + assay(m1_s2, 'U'), ncol = 2, dimnames = list(NULL, c('sample1', 'sample2')))))
  m1_s3 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(assays = m1_s3$assays, rowData = m1_s3$rowData, 
                  colData = m1_s3$colData)
  expect_that(getCoverage(m1_s3), is_identical_to(matrix(assay(m1_s3, 'M') + assay(m1_s3, 'U'), ncol = 3, dimnames = list(NULL, c('sample1', 'sample2', 'sample3')))))
  m1_s4 <- make_MTuples_or_CoMeth_data(m = 1L, n = 10L, s = 4L, sim_counts = TRUE)
  m1_s4 <- CoMeth(assays = m1_s4$assays, rowData = m1_s4$rowData, 
                  colData = m1_s4$colData)
  expect_that(getCoverage(m1_s4), is_identical_to(matrix(assay(m1_s4, 'M') + assay(m1_s4, 'U'), ncol = 4, dimnames = list(NULL, c('sample1', 'sample2', 'sample3', 'sample4')))))
})

## TODO: Test that duplicated works when m > = 3
test_that("'duplicated' works", {
  m4_s2 <- make_MTuples_or_CoMeth_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(assays = m4_s2$assays, rowData = m4_s2$rowData, 
                  colData = m4_s2$colData)
  start(m4_s2) <- 1L
  m4_s2@rowData@extraPos <- matrix(c(rep(3L, nrow(m4_s2)), rep(7L, nrow(m4_s2))), ncol = 2)
  end(m4_s2) <- 10L
  expect_that(any(duplicated(m4_s2)), is_true())
})