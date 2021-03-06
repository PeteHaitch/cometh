% Generated by roxygen2 (4.0.1): do not edit by hand
\docType{class}
\name{CoMeth-class}
\alias{CoMeth-class}
\title{An S4 class to store co-methylation patterns at m-tuples of genomic
positions.}
\description{
An S4 class to store co-methylation patterns at m-tuples of genomic
positions.
}
\details{
The \code{CoMeth} class is based on the
\code{\link[GenomicRanges]{SummarizedExperiment}} class. The main difference
is that rather than using a \code{\link[GenomicRanges]{GRanges}} object as
the \code{rowData}, a \code{CoMeth} object uses an \code{\link{MTuples}}
object.

The assays of a \link{CoMeth} object are the counts of how many times each
co-methylation pattern is observed for that m-tuple in each sample.
For example, the possible co-methylation patterns of 2-tuples are 'MM', 'MU',
'UM' and 'UU' and thus there are four assays of the same names.
}
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

