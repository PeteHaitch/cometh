\name{read.comethylation}
\alias{read.comethylation}
\title{Parsing \code{tsv} output from \code{comethylation.py}}
\usage{
read.comethylation(file, m, methylation_type, sample_name, quiet = FALSE,
  seqlengths = NULL, seqinfo = NULL)
}
\arguments{
  \item{file}{The \code{.tsv} file created by
  comethylation.py. This file may be compressed with gzip
  or bzip2}

  \item{m}{The size of the m-tuples}

  \item{methylation_type}{A character vector with the type
  of methylation}

  \item{sample_name}{The sample name}

  \item{quiet}{logical: if \code{FALSE} (default),
  read.comethylation() will print a line, saying how many
  items (m-tuples) have been read.}

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
Read the \code{tsv} output file from
\code{comethylation.py} and construct a
\code{\link{CoMeth}} object.
}
\note{
The compression of \code{file} is determined by the file
extension: \code{.gz} for compressed with gzip and
\code{.bz2} for compressed with bzip2. Otherwise the file
is assumed to be uncompressed.
}
\examples{
cat("TODO")
}
\seealso{
\code{\link{CoMeth}}
}

