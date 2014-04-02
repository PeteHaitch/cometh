#### DESCRIPTION ####
# Peter Hickey
# 14/03/2014
# The "[" method is not working as expected for my simple class that extends theSummarizedExperiment class

library(GenomicRanges)

#### Define MyClass, basically a SummarizedExperiment with an "extra_slot", which is a numeric vector of the same length as the number of rows in the rowData ####
.MyClass <- setClass("MyClass", representation(extra_slot = "vector"), contains = "SummarizedExperiment")

MyClass <- function(assays, gr, extra_slot){
  se <- SummarizedExperiment(assays = assays, rowData = gr)
  mc <- .MyClass(se, extra_slot = extra_slot)
  mc
}


.MyClass2 <- setClass("MyClass2", contains = "SummarizedExperiment")

MyClass2 <- function(assays, gr){
  se <- SummarizedExperiment(assays = assays, rowData = gr)
  mc <- .MyClass2(se)
  mc
}

#### Define the "[" method for MyClass objects. Should basically do the same thing as "[" for a SummarizedExperiment but also subset the extra_slot ####
setMethod("[", c("MyClass", "ANY", "missing"),
          function(x, i, j, ..., drop = FALSE)
          {
            if (missing(i) && missing(j)){
              return(x)
            }
            if (missing(i)){
              i <- seq_len(nrow(x))
            }
            if (missing(j)){
              j <- seq_len(ncol(x))
            }                      
            initialize(x, as(x, "SummarizedExperiment")[i, j, ..., drop = drop], extra_slot = x@extra_slot[i]) # Is there any reason not to just use x[i, j, ..., drop = FALSE] instead of as(x, "SummarizedExperiment")[i, j, ..., drop = drop] ?
          })

#### Example data ####
nrows <- 20
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colnames(counts) <- letters[1:6]
gr <- GRanges('chr1', IRanges(floor(runif(nrows, 1e5, 1e6)), width=100))
extra_slot <- runif(n = length(gr))
gr2 <- gr
values(gr2)$pos2 <- extra_slot

ex <- MyClass(assays = SimpleList(counts), gr = gr, extra_slot = extra_slot)
ex2 <- MyClass2(assays = SimpleList(counts), gr = gr2)

#### Testing the "[" method on MyClass objects
ex[1:10, ]@extra_slot # As expected contains 10 elements
ex[1:10, 2]@extra_slot # Contains all 20 elements; why?

initialize(ex, as(ex, "SummarizedExperiment")[1:10, 2, drop = FALSE], extra_slot = ex@extra_slot[1:10])@extra_slot # As (un)expected contains 10 elements; why?
