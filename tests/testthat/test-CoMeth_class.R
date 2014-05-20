## TODO: It's possible to have two SummarizedExperiments (well, CoMeth) objects that have identical slots, all.equal returns TRUE but identical returns FALSE. 
## E.g. a <- make_test_data(3L, 10L, 3L, sim_counts = TRUE), identical(combine(a[, 1], a[, 2], a[, 3]), a)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth constructor.
###
context('CoMeth constructor')

test_that("CoMeth constructor returns an CoMeth object when m = 1 for a single sample", {
  m1_s1 <- make_test_data(m = 1L, n = 10L, s = 1L, sim_counts = TRUE)
  m1_s1 <- CoMeth(sample_names = m1_s1$sample_names, methylation_type = m1_s1$methylation_type, counts = m1_s1$counts, seqnames = m1_s1$seqnames, pos = m1_s1$pos, seqinfo = m1_s1$seqinfo)
  expect_that(m1_s1, is_a("CoMeth"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 1 for multiple samples", {
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(m1_s3, is_a("CoMeth"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for a single sample", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(m2_s1, is_a("CoMeth"))
})

test_that("CoMeth constructor returns an CoMeth object when m = 2 for multiple samples", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(m2_s3, is_a("CoMeth"))
})

test_that("CoMeth constructor returns an CoMeth object when m >= 2 for a single sample", {
  m3_s1 <- make_test_data(m = 3L, n = 10L, s = 1L, sim_counts = TRUE)
  m3_s1 <- CoMeth(sample_names = m3_s1$sample_names, methylation_type = m3_s1$methylation_type, counts = m3_s1$counts, seqnames = m3_s1$seqnames, pos = m3_s1$pos, seqinfo = m3_s1$seqinfo)
  expect_that(m3_s1, is_a("CoMeth"))
  m4_s1 <- make_test_data(m = 4L, n = 10L, s = 1L, sim_counts = TRUE)
  m4_s1 <- CoMeth(sample_names = m4_s1$sample_names, methylation_type = m4_s1$methylation_type, counts = m4_s1$counts, seqnames = m4_s1$seqnames, pos = m4_s1$pos, seqinfo = m4_s1$seqinfo)
  expect_that(m4_s1, is_a("CoMeth"))
})

test_that("CoMeth constructor returns an CoMeth object when m >= 2 for multiple samples", {
  m3_s3 <- make_test_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(sample_names = m3_s3$sample_names, methylation_type = m3_s3$methylation_type, counts = m3_s3$counts, seqnames = m3_s3$seqnames, pos = m3_s3$pos, seqinfo = m3_s3$seqinfo)
  expect_that(m3_s3, is_a("CoMeth"))
  m4_s3 <- make_test_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(sample_names = m4_s3$sample_names, methylation_type = m4_s3$methylation_type, counts = m4_s3$counts, seqnames = m4_s3$seqnames, pos = m4_s3$pos, seqinfo = m4_s3$seqinfo)
  expect_that(m4_s3, is_a("CoMeth"))
})

test_that("CoMeth parameter checking works: 'sample_names'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote("sample_names"), " missing.")))
  expect_that(CoMeth(sample_names = unlist(m2_s3$sample_names), methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('sample_names'), " must be a ", sQuote('CharacterList'))))
  expect_that(CoMeth(sample_names = as(rep('sample_1', 3), "CharacterList"), methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Each element of ", sQuote("sample_names"), " must be unique.")))
})

test_that("CoMeth parameter checking works: 'methylation_type'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('methylation_type'), " missing.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = as.vector(m2_s3$methylation_type), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('methylation_type'), " must be a ", sQuote('CharacterList'))))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = unlist(m2_s3$methylation_type), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('methylation_type'), " must be a ", sQuote('CharacterList'))))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type[1:2], counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Names of ", sQuote('methylation_type'), " must match those in ", sQuote('sample_names'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = 'CpG', sample3 = "CHG"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('methylation_type'), " for each sample must be ", sQuote('CG'), ", ", sQuote('CHG'), ", ", sQuote('CHH'), " or ", sQuote('CNN'), ", or a vector of some combination of these")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = c('CG', 'CHH'), sample3 = "CHG"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), gives_warning(paste0("The supplied ", sQuote('methylation_type'), " parameter says that samples contain data for different methylation types. The union of these methylation types will be used as the ", sQuote('methylation_type'), " of the returned ", sQuote('CoMeth'), " object.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = CharacterList(sample1 = 'CG', sample2 = c('CG', 'CHH'), sample3 = "CHH"), counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), gives_warning(paste0("The supplied ", sQuote('methylation_type'), " parameter says that samples contain data for different methylation types. The union of these methylation types will be used as the ", sQuote('methylation_type'), " of the returned ", sQuote('CoMeth'), " object.")))
})

test_that("CoMeth parameter checking works: 'counts'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('counts'), " missing.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = as.list(m2_s3$counts), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('counts'), " must be a ", sQuote('DataFrameList'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts[1:2], seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Names of ", sQuote('counts'), " must match those in ", sQuote('sample_names'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){cbind(x[, -3])}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('ncol\\(counts\\)'), " must be identical for all elements of ", sQuote('counts'), " and should be a power of 2.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = c(m2_s3$counts[1:2], DataFrameList(sample3 = m2_s3$counts[[3]][, -3])), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('ncol\\(counts\\)'), " must be identical for all elements of ", sQuote('counts'), " and should be a power of 2.")))  
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){DataFrame(lapply(x, as.numeric))}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Class of columns in each ", sQuote('DataFrame'), " element of ", sQuote("counts"), " must be: ", paste0(rep(sQuote('integer'), 2 ^ 2), collapse = ", "), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = endoapply(m2_s3$counts, function(x){colnames(x) <- rev(colnames(x)); x}), seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('m'), " is set to ", 2L, " so colnames for all elements of ", sQuote('counts'), " must be: ", paste0(sQuote(.make_m_tuple_names(2L)), collapse = ", "), ".")))
})

test_that("CoMeth parameter checking works: 'seqnames'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('seqnames'), " missing.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = as.list(m2_s3$seqnames), pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('seqnames'), " must be an ", sQuote("RleList"), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames[1:2], pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Names of ", sQuote('seqnames'), " must match those in ", sQuote('sample_names'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = endoapply(m2_s3$seqnames, function(x){x[-1]}), pos = m2_s3$pos, seqinfo = m2_s3$seqinfo), throws_error(paste0("Length of each ", sQuote('Rle'), " element of ", sQuote("seqnames"), " must be identical to ", sQuote('nrow\\(pos\\)'), ".")))
})

test_that("CoMeth parameter checking works: 'pos'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('pos'), " missing.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = as.list(m2_s3$pos), seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('pos'), " must be a ", sQuote("DataFrameList"), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos[1:2], seqinfo = m2_s3$seqinfo), throws_error(paste0("Names of ", sQuote('pos'), " must match those in ", sQuote('sample_names'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = endoapply(m2_s3$pos, function(x){x[-1, ]}), seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('nrow\\(counts\\)'), " must be identical to ", sQuote('nrow\\(pos\\)'), ".")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = endoapply(m2_s3$pos, function(x){colnames(x) <- rev(colnames(x)); x}), seqinfo = m2_s3$seqinfo), throws_error(paste0(sQuote('m'), " is set to ", 2L, " so colnames for all elements of ", sQuote('pos'), " must be: ", paste0(sQuote(paste0('pos', seq_len(2L))), collapse = ', '), ".")))
})

test_that("CoMeth parameter checking works: 'seqinfo'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos), throws_error(paste0(sQuote('seqinfo'), " missing.")))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = as.data.frame(m2_s3$seqinfo)), throws_error(paste0(sQuote('seqinfo'), " must be a ", sQuote('Seqinfo'), " object.")))
})

test_that("CoMeth parameter checking works: 'strand'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = Rle()), throws_error(paste0(sQuote('strand'), " must be an ", sQuote('RleList'), ".")))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3 = Rle(factor('*'), length(m2_s3$seqnames[[3]])))
  strand <- endoapply(strand, function(x){levels(x) <- c('*', '+', '-', 'AAA'); x})
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = strand), throws_error("invalid strand levels in 'x': AAA"))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3 = Rle(factor('*'), length(m2_s3$seqnames[[3]])))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = endoapply(strand, function(x){x[-1]})), throws_error(paste0("Length of each ", sQuote('Rle'), " element in ", sQuote('strand'), " must equal the number of rows of its corresponding element in ", sQuote('pos'), ".")))
  strand <- RleList(sample1 = Rle(factor('*'), length(m2_s3$seqnames[[1]])), sample2 = Rle(factor('*'), length(m2_s3$seqnames[[2]])), sample3= Rle(factor('+'), length(m2_s3$seqnames[[1]])))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, strand = strand), gives_warning(paste0("The supplied ", sQuote('strand'), " argument contains ", sQuote('\\*'), " as well as at least one of ", sQuote('\\+'), " or ", sQuote('\\-'), ". m-tuples will not be combined across samples if they are on different strands.")))
})

test_that("CoMeth parameter checking works: 'colData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = data.frame()), throws_error(paste0(sQuote('colData'), ", if supplied, must be a ", sQuote('DataFrame'), ".")))
  colData <- DataFrame(a = c(1:3), row.names = paste0('sample', letters[1:3]))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = colData), throws_error(paste0(sQuote('row.names\\(colData\\)'), " must be identical to ", sQuote('sample_names'), ".")))
  colData <- DataFrame(cancer = c(TRUE, TRUE, FALSE), row.names = paste0('sample', 1:3))
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, colData = colData), is_a("CoMeth"))
})

test_that("CoMeth parameter checking works: 'exptData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, exptData = list()), throws_error(paste0(sQuote('exptData'), ", if supplied, must be a ", sQuote('SimpleList'), ".")))
  exptData <- SimpleList(assay_type = "methylC-seq")
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, exptData = exptData), is_a("CoMeth"))
})

test_that("CoMeth parameter checking works: '...'", {
  m2_s1 <-  make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  cgi <- rep(TRUE, 10)
  expect_that(CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo, cgi = cgi), is_a("CoMeth"))
  expect_that(mcols(CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo, cgi = cgi)), is_identical_to(DataFrame(cgi = cgi)))
})

test_that("CoMeth parameter checking works: 'verbose'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  expect_that(CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo, verbose = TRUE), is_a("CoMeth"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test CoMeth validity
###
context("CoMeth validity")

test_that("CoMeth validity checking works on good object", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(validObject(m2_s3), is_true())
})

test_that("CoMeth validity checking works: '.valid.CoMeth.counts'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  assay(m2_s3, 'MM') <- matrix(-1, ncol = ncol(m2_s3), nrow = nrow(m2_s3))
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('counts'), " cannot have negative entries.")))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  names(assays(m2_s3)) <- paste0((names(assays(m2_s3))), 'a')
  expect_that(validObject(m2_s3), throws_error(paste0("assay names must include: ", paste0(sQuote(.make_m_tuple_names(2L)), collapse = ', '), ".")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.methylation_type'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  colData(m2_s3) <- colData(m2_s3)[, -1]
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('colData'), " of ", sQuote('CoMeth'), " must contain column ", sQuote('methylation_type'), " once and only once.")))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  colData(m2_s3)$methylation_type <- rep('CpG', 3)
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('methylation_type'), " for each sample must be ", sQuote('CG'), ", ", sQuote('CHG'), ", ", sQuote('CHH'), " or ", sQuote('CNN'), ", or some combination of these, e.g., ", sQuote("CG/CHG"), ".")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.rowData'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  rowData(m2_s3) <- as(rowData(m2_s3), "GRanges")
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('rowData\\(CoMeth\\)'), " must be an ", sQuote('MTuples'), " object.")))
})

test_that("CoMeth validity checking works: '.valid.CoMeth.noDuplicates'", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  rowData(m2_s3)[2] <- rowData(m2_s3)[1]
  expect_that(validObject(m2_s3), throws_error(paste0(sQuote('CoMeth'), " object cannot contain duplicate m-tuples.")))
})

context("CoMeth methods")

test_that("'getPos' works", {
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  pos <- unlist(m2_s3$pos)
  pos <- unique(pos)
  row.names(pos) <- NULL
  pos <- as.matrix(pos)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  ## NOTE: "expect_that(getPos(m2_s3), is_identical_to(pos))" will fail if 
  ## positions are not unique across samples
  ## The below test, which is a bit convoluted, should always work.
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3$pos[[2]] <- rbind(m2_s3$pos[[1]][4:10, ], m2_s3$pos[[2]][1:3, ])
  pos <- cbind(as.vector(unlist(m2_s3$seqnames, use.names = FALSE)), as.matrix(unlist(m2_s3$pos)))
  pos <- unique(pos)
  row.names(pos) <- NULL
  pos <- pos[, -1]
  pos <- matrix(as.integer(pos), ncol = 2, dimnames = list(NULL, paste0('pos', 1:2)))
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(getPos(m2_s3)[order(getPos(m2_s3)[, 'pos1'], getPos(m2_s3)[, 'pos2']), ], is_identical_to(pos[order(pos[, 'pos1'], pos[, 'pos2']), ]))
})

test_that("'getM' works", {
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(getM(m1_s3), is_identical_to(1L))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(getM(m2_s3), is_identical_to(2L))
  m3_s3 <- make_test_data(m = 3L, n = 10L, s = 3L, sim_counts = TRUE)
  m3_s3 <- CoMeth(sample_names = m3_s3$sample_names, methylation_type = m3_s3$methylation_type, counts = m3_s3$counts, seqnames = m3_s3$seqnames, pos = m3_s3$pos, seqinfo = m3_s3$seqinfo)
  expect_that(getM(m3_s3), is_identical_to(3L))
  m4_s3 <- make_test_data(m = 4L, n = 10L, s = 3L, sim_counts = TRUE)
  m4_s3 <- CoMeth(sample_names = m4_s3$sample_names, methylation_type = m4_s3$methylation_type, counts = m4_s3$counts, seqnames = m4_s3$seqnames, pos = m4_s3$pos, seqinfo = m4_s3$seqinfo)
  expect_that(getM(m4_s3), is_identical_to(4L))
})

test_that("'sampleNames' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', 1)))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', 1:2)))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', 1:3)))
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', 1:4)))
})

test_that("'sampleNames<-' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  sampleNames(m2_s1) <- paste0('sample', letters[1])
  expect_that(sampleNames(m2_s1), is_identical_to(paste0('sample', letters[1])))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  sampleNames(m2_s2) <- paste0('sample', letters[1:2])
  expect_that(sampleNames(m2_s2), is_identical_to(paste0('sample', letters[1:2])))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  sampleNames(m2_s3) <- paste0('sample', letters[1:3])
  expect_that(sampleNames(m2_s3), is_identical_to(paste0('sample', letters[1:3])))
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  sampleNames(m2_s4) <- paste0('sample', letters[1:4])
  expect_that(sampleNames(m2_s4), is_identical_to(paste0('sample', letters[1:4])))
})

test_that("'length' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(length(m2_s1), is_identical_to(1L))
})

test_that("'nrow' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(nrow(m2_s1), is_identical_to(10L))
})

test_that("'ncol' works", {
  m2_s1 <- make_test_data(m = 2L, n = 10L, s = 1L, sim_counts = TRUE)
  m2_s1 <- CoMeth(sample_names = m2_s1$sample_names, methylation_type = m2_s1$methylation_type, counts = m2_s1$counts, seqnames = m2_s1$seqnames, pos = m2_s1$pos, seqinfo = m2_s1$seqinfo)
  expect_that(ncol(m2_s1), is_identical_to(1L))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(ncol(m2_s2), is_identical_to(2L))
  m2_s3 <- make_test_data(m = 2L, n = 10L, s = 3L, sim_counts = TRUE)
  m2_s3 <- CoMeth(sample_names = m2_s3$sample_names, methylation_type = m2_s3$methylation_type, counts = m2_s3$counts, seqnames = m2_s3$seqnames, pos = m2_s3$pos, seqinfo = m2_s3$seqinfo)
  expect_that(ncol(m2_s3), is_identical_to(3L))
  m2_s4 <- make_test_data(m = 2L, n = 10L, s = 4L, sim_counts = TRUE)
  m2_s4 <- CoMeth(sample_names = m2_s4$sample_names, methylation_type = m2_s4$methylation_type, counts = m2_s4$counts, seqnames = m2_s4$seqnames, pos = m2_s4$pos, seqinfo = m2_s4$seqinfo)
  expect_that(ncol(m2_s4), is_identical_to(4L))
})

test_that("'granges' is just an alias of 'rowData'", {
  ## granges(CoMeth) doesn't really make sense because the rowData is MTuples rather than GRanges.
  ## As of GenomicRanges_1.14.4, granges(SummarizedExperiment) is simply an alias of rowData(SummarizedExperiment), hence, by inheritance, granges(CoMeth) is just an alias of rowData(CoMeth).
  ## This test is in case this behaviour changes in future versions of GenomicRanges.
  m5_s4 <- make_test_data(m = 5L, n = 10L, s = 4L, sim_counts = TRUE)
  m5_s4 <- CoMeth(sample_names = m5_s4$sample_names, methylation_type = m5_s4$methylation_type, counts = m5_s4$counts, seqnames = m5_s4$seqnames, pos = m5_s4$pos, seqinfo = m5_s4$seqinfo)
  expect_that(granges(m5_s4), is_identical_to(rowData(m5_s4)))
})

test_that("'cbind' works", {
  # TODO: Test cbind fails if objects have different seqinfo
  # TODO: Test cbind for other failure modes
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(cbind(m1_s2[, 1], m1_s2[, 2]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m1_s2, m1_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(cbind(m2_s2[, 1], m2_s2[, 2]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m2_s2, m2_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(cbind(m3_s2[, 1], m3_s2[, 2]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m3_s2, m3_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(cbind(m4_s2[, 1], m4_s2[, 2]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(cbind(m4_s2, m4_s2), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with duplicate ", sQuote('sample_names'), ".")))
  expect_that(cbind(m2_s2[, 1], m4_s2[, 2]), throws_error(paste0("Cannot ", sQuote('cbind'), " ", sQuote('CoMeth'), " objects with the different sized m-tuples, that is, different ", sQuote('m'), ".")))
})

test_that("'rbind' works", {
  # TODO: Test rbind fails if objects have different seqinfo
  # TODO: Test rbind for other failure modes
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(rbind(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(rbind(m1_s2, m1_s2), throws_error(paste0("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(rbind(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m2_s2, m2_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(rbind(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).  expect_that(rbind(m3_s2, m3_s2), throws_error("Cannot ‘rbind’ ‘CoMeth’ object with any identical m-tuples."))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(rbind(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(rbind(m4_s2, m4_s2), throws_error(paste0("Cannot ", sQuote('rbind'), " ", sQuote('CoMeth'), " object with any identical m-tuples.")))
})

test_that("'combine' works", {
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m1_s2[seq.int(from = 1, to = floor(nrow(m1_s2) / 2), by = 1), ], m1_s2[seq.int(from = floor(nrow(m1_s2) / 2) + 1, to = nrow(m1_s2), by = 1), ]), is_equivalent_to(m1_s2))
  expect_that(combine(m1_s2, m1_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
  m2_s2 <- make_test_data(m = 2L, n = 10L, s = 2L, sim_counts = TRUE)
  m2_s2 <- CoMeth(sample_names = m2_s2$sample_names, methylation_type = m2_s2$methylation_type, counts = m2_s2$counts, seqnames = m2_s2$seqnames, pos = m2_s2$pos, seqinfo = m2_s2$seqinfo)
  expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m2_s2[seq.int(from = 1, to = floor(nrow(m2_s2) / 2), by = 1), ], m2_s2[seq.int(from = floor(nrow(m2_s2) / 2) + 1, to = nrow(m2_s2), by = 1), ]), is_equivalent_to(m2_s2))
  expect_that(combine(m2_s2, m2_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
  m3_s2 <- make_test_data(m = 3L, n = 10L, s = 2L, sim_counts = TRUE)
  m3_s2 <- CoMeth(sample_names = m3_s2$sample_names, methylation_type = m3_s2$methylation_type, counts = m3_s2$counts, seqnames = m3_s2$seqnames, pos = m3_s2$pos, seqinfo = m3_s2$seqinfo)
  expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m3_s2[seq.int(from = 1, to = floor(nrow(m3_s2) / 2), by = 1), ], m3_s2[seq.int(from = floor(nrow(m3_s2) / 2) + 1, to = nrow(m3_s2), by = 1), ]), is_equivalent_to(m3_s2))
  expect_that(combine(m3_s2, m3_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2)) # Can't use is_identical_to because the assays are stored as a ReferenceClass (at least, I think this is the reason).
  expect_that(combine(m4_s2[seq.int(from = 1, to = floor(nrow(m4_s2) / 2), by = 1), ], m4_s2[seq.int(from = floor(nrow(m4_s2) / 2) + 1, to = nrow(m4_s2), by = 1), ]), is_equivalent_to(m4_s2))
  expect_that(combine(m4_s2, m4_s2), throws_error(paste0(sQuote('combine'), " failed when trying to ", sQuote('rbind'), " the intermediate ", sQuote('CoMeth'), " objects.")))
})

test_that("'getCoverage' works", {
  m1_s1 <- make_test_data(m = 1L, n = 10L, s = 1L, sim_counts = TRUE)
  m1_s1 <- CoMeth(sample_names = m1_s1$sample_names, methylation_type = m1_s1$methylation_type, counts = m1_s1$counts, seqnames = m1_s1$seqnames, pos = m1_s1$pos, seqinfo = m1_s1$seqinfo)
  expect_that(getCoverage(m1_s1), is_identical_to(matrix(assay(m1_s1, 'M') + assay(m1_s1, 'U'), ncol = 1, dimnames = list(NULL, 'sample1'))))
  m1_s2 <- make_test_data(m = 1L, n = 10L, s = 2L, sim_counts = TRUE)
  m1_s2 <- CoMeth(sample_names = m1_s2$sample_names, methylation_type = m1_s2$methylation_type, counts = m1_s2$counts, seqnames = m1_s2$seqnames, pos = m1_s2$pos, seqinfo = m1_s2$seqinfo)
  expect_that(getCoverage(m1_s2), is_identical_to(matrix(assay(m1_s2, 'M') + assay(m1_s2, 'U'), ncol = 2, dimnames = list(NULL, c('sample1', 'sample2')))))
  m1_s3 <- make_test_data(m = 1L, n = 10L, s = 3L, sim_counts = TRUE)
  m1_s3 <- CoMeth(sample_names = m1_s3$sample_names, methylation_type = m1_s3$methylation_type, counts = m1_s3$counts, seqnames = m1_s3$seqnames, pos = m1_s3$pos, seqinfo = m1_s3$seqinfo)
  expect_that(getCoverage(m1_s3), is_identical_to(matrix(assay(m1_s3, 'M') + assay(m1_s3, 'U'), ncol = 3, dimnames = list(NULL, c('sample1', 'sample2', 'sample3')))))
  m1_s4 <- make_test_data(m = 1L, n = 10L, s = 4L, sim_counts = TRUE)
  m1_s4 <- CoMeth(sample_names = m1_s4$sample_names, methylation_type = m1_s4$methylation_type, counts = m1_s4$counts, seqnames = m1_s4$seqnames, pos = m1_s4$pos, seqinfo = m1_s4$seqinfo)
  expect_that(getCoverage(m1_s4), is_identical_to(matrix(assay(m1_s4, 'M') + assay(m1_s4, 'U'), ncol = 4, dimnames = list(NULL, c('sample1', 'sample2', 'sample3', 'sample4')))))
})

## TODO: Test that duplicated works when m > = 3
test_that("'duplicated' works", {
  m4_s2 <- make_test_data(m = 4L, n = 10L, s = 2L, sim_counts = TRUE)
  m4_s2 <- CoMeth(sample_names = m4_s2$sample_names, methylation_type = m4_s2$methylation_type, counts = m4_s2$counts, seqnames = m4_s2$seqnames, pos = m4_s2$pos, seqinfo = m4_s2$seqinfo)
  start(m4_s2) <- 1
  end(m4_s2) <- 10
  m4_s2@rowData@extraPos <- matrix(c(rep(3, nrow(m4_s2)), rep(7, nrow(m4_s2))), ncol = 2)
  expect_that(any(duplicated(m4_s2)), is_true())
})