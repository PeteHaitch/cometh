% Generated by roxygen2 (4.0.1): do not edit by hand
\name{MethylationLociSet}
\alias{MethylationLociSet}
\title{The constructor function for MethylationLociSet objects.}
\usage{
MethylationLociSet(seqnames = Rle(), ranges = IRanges(), strand = Rle("*",
  length(seqnames)), seqinfo, methylation_type, ...)
}
\description{
A MethylationLociSet contains the positions of all methylation loci in the
sample. This is usually the set of all methylation loci in the reference
genome.
Most users will construct these objects using the \code{\link{makeMLS}}
function.
}

