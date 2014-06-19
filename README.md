# cometh
[![Build Status](https://travis-ci.org/PeteHaitch/cometh.png?branch=CoMeth_VIRTUAL_class)](https://travis-ci.org/PeteHaitch/cometh)

An R package with tools for analysing, managing and visualising DNA methlyation m-tuples and co-methylation data. DNA methlyation m-tuples are the methylation patterns at sets of adjacent methylation loci on the same DNA fragment and, loosely speaking, co-methylation is the correlation structure of DNA methylation.

__This package is in early development. It requires the use of the development version of Bioconductor ([instructions here](http://bioconductor.org/developers/how-to/useDevel/)).__ It can be installed (provided it is passing the Travis CI build) with:

```R
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite(c('Biobase', 'S4Vectors', 'GenomeInfoDb', 'IRanges', 'GenomicRanges', 'BSgenome'))
devtools::install_github("PeteHaitch/cometh")
```