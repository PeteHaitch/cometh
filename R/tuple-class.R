# Based on Martin Morgan's reply to my question on Bioc-Devel (https://stat.ethz.ch/pipermail/bioc-devel/2014-February/005242.html)
library(GenomicRanges)

.Tuple <- setClass("Tuple", representation(m="integer"), contains="SimpleList", )


Tuple <- function(seqnames, pos, count){
  ## note constructor: base class(es) as first and unnamed arg
  .Tuple(SimpleList(seqnames=seqnames, pos=pos, count=count), m=ncol(pos))
}

setMethod(show, signature(object = "Tuple"), function(object) {
  cat(length(object$seqnames), " x ", object@m, "-", sep="")
  callNextMethod()
})

#tuple <- Tuple(test_data$seqnames, matrix(unlist(test_data[2:4]), ncol=3),  matrix(unlist(test_data[5:12]), ncol=8))


.TupleList <- setClass("TupleList", contains="CompressedList", prototype=prototype(elementType="Tuple"))

setMethod("[", c("Tuple", "ANY", "missing"),
          function(x, i, j, ..., drop=FALSE)
          {
            if (!missing(drop))
              warning("'drop' ignored when subsetting ", sQuote(class(x)))
            initialize(x, SimpleList(seqnames=x$seqnames[i],
                                     pos=x$pos[i,,drop=FALSE], count=x$count[i,,drop=FALSE]))
          })