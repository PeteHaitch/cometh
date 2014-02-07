\name{CoMeth}
\alias{CoMeth}
\title{The constructor function for cometh objects}
\usage{
CoMeth(seqnames, pos, counts, m, methylation_type, sample_name, strand = "*",
  seqlengths = NULL, seqinfo = NULL)
}
\arguments{
  \item{seqnames}{Rle object, character vector, or factor
  containing the sequence names.}

  \item{pos}{A numeric matrix storing the positions of each
  m-tuple as a row. Thus nrow(pos) = the number of m-tuples
  and ncol(pos) = m.}

  \item{counts}{A numeric matrix storing the number of
  times each co-methylation pattern is observed at each
  m-tuple. Thus nrow(pos) = the number of m-tuples and
  ncol(pos) = the number of possible co-methylation
  patterns = 2 ^ m.}

  \item{m}{An integer storing the size of the m-tuples,
  i.e. the "m" in m-tuple.}

  \item{methylation_type}{A character vector storing the
  type of methylation loci for these m-tuples. Possible
  values are "CG", "CHG", "CHH" or "CNN". Multiple values
  can be specified, e.g. c("CG", "CHG").}

  \item{sample_name}{The name of the sample.}

  \item{strand}{An (optional) Rle object, character vector,
  or factor containing the strand information.}

  \item{seqlengths}{The sequence lengths of the reference
  genome of the sample. Must be an integer vector named
  with the sequence names and containing the lengths (or
  NA) for each level(seqnames).}

  \item{seqinfo}{An (optional)
  \code{\link[GenomicRanges]{Seqinfo}} object containing
  information about the reference genome of the sample}
}
\value{
A \code{\link{CoMeth}} object
}
\description{
Construct a \code{CoMeth} object. Users will generally use
the \code{\link{read.comethylation}} to construct a
\code{CoMeth} object.
}
\examples{
cat("TODO")
}
\seealso{
\code{\link{read.comethylation}}
}

