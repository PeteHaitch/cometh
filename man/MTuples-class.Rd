% Generated by roxygen2 (4.0.1): do not edit by hand
\docType{class}
\name{MTuples-class}
\alias{MTuples-class}
\title{An S4 class to represent m-tuples of genomic positions.}
\description{
An S4 class to represent m-tuples of genomic positions.
}
\details{
The \code{MTuples} class extends the \code{\link[GenomicRanges]{GRanges}}
class by adding the \code{extraPos} slot (see below for details). An m-tuple
is a tuple of individual basepairs that are on the same chromosome, where
'm' is the number of positions in the tuple. For example,
(chr1:30, chr1:33, chr1:40) is a 3-tuple of the positions on chromosome 1.
Note the strand of the m-tuple is optional.

Internally, this example 3-tuple is stored as a GRanges object with the first
and last positions of the m-tuple stored as the \code{start} and \code{end}
of the GRanges interval, respectively. That is,
\code{GRanges('chr1', IRanges(start = 30, end = 40))}. The "extra" position,
chr1:33, is stored in the \code{extraPos} matrix.
}
\section{Slots}{

\describe{
\item{\code{extraPos}}{A numeric matrix storing "extra" positions in m-tuples,
provided m >= 3. If m = 1 or 2, \code{extraPos} is a matrix of \code{NA}.
The \code{extraPos} matrix has as many rows as there are m-tuples in the
\code{MTuples} object.}
}}
\section{Constructor}{

\strong{TODO}: Insert help for constructor method.
}

\section{Coercion}{

\strong{TODO}: Insert help for any coerction methods.
}

\section{Accessors}{

\strong{TODO}: Insert help for any accessor methods.
}

\section{Splitting and combining}{

\strong{TODO}: Insert help for any splitting and combining methods.
}

\section{Subsetting}{

\strong{TODO}: Describe any subsetting methods.
}

\section{Filtering}{

\strong{TODO}: Describe any filtering methods.
}

\section{Methods based on findOverlaps}{

\strong{TODO} Insert help for any findOverlaps-based methods.
}

\section{Other methods}{

\strong{TODO}: Describe any other methods.
}

