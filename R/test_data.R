## TODO: Move to test/testthat

make_test_data <- function(m, n){
  test_data <- lapply(cbind(data.frame(chr = c(rep('chr1', 0.6 * n ), rep('chr2', 0.3 * n), rep('chrX', 0.1 * n)), stringsAsFactors = FALSE), as.data.frame(matrix(sort(sample(1:(n * m * 2), m * n, replace = FALSE)), ncol = m, byrow = T, dimnames = list(NULL, paste0('pos', 1:m)))), as.data.frame(matrix(rpois(2^m * n, 4), ncol = 2^m, dimnames = list(NULL, sort(do.call(paste0, expand.grid(lapply(seq_len(m), function(x){c('M', 'U')})))))))), function(x){x})
  seqnames <- Rle(as.character(test_data[['chr']]))
  pos <- DataFrame(lapply(test_data[grep('pos', names(test_data))], as.vector))
  counts <- DataFrame(lapply(test_data[grep('[MU]', names(test_data))], as.vector))
  return(list(seqnames = seqnames, pos = pos, counts = counts))
}

# Create some test data
set.seed(666)
m <- 3L
a <- make_test_data(m, 200)
b <- make_test_data(m, 100)
d <- make_test_data(m, 1000)
e <- list(pos = DataFrameList(a = a$pos, b = b$pos, d= d$pos), counts = DataFrameList(a = a$counts, b = b$counts, d = d$counts))
seqnames <- RleList(a = a$seqnames, b = b$seqnames, d = d$seqnames)
#seqnames <- RleList(a = a$seqnames)
pos <- DataFrameList(a = a$pos, b = b$pos, d = d$pos)
#pos <- DataFrameList(a = a$pos)
counts <-  DataFrameList(a = a$counts, b = b$counts, d = d$counts)
#counts <-  DataFrameList(a = a$counts)
strand <- RleList(a = Rle('*', 200), b = Rle('*', 100), d = Rle('*', 1000))
#strand <- RleList(a = Rle('*', 200))
sample_names <- as(c('a', 'b', 'd'), "CharacterList")
#sample_names <- as(c('a'), "CharacterList")
methylation_type <- CharacterList(a = 'CG', b = 'CG/CHH', d = 'CG')
#methylation_type <- CharacterList(a = 'CG')
seqinfo <- Seqinfo(seqnames = c('chr1', 'chr2', 'chrX'), seqlengths = c(249250621, 243199373, 155270560), genome = 'hg19')

z <- CoMeth(counts = counts, seqnames = seqnames, pos = pos, strand = strand, seqinfo = seqinfo, sample_names = sample_names, methylation_type = methylation_type)
w <- rowData(z)

