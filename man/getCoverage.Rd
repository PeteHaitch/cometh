\name{getCoverage}
\alias{getCoverage}
\title{Get the \code{coverage} from a \code{\link{CoMeth}} object.}
\usage{
getCoverage(CoMeth)
}
\arguments{
  \item{CoMeth}{A \code{\link{CoMeth}} object}
}
\value{
A numeric vector
}
\description{
The \code{coverage} of a \code{\link{CoMeth}} object is the
number of reads for each m-tuple, i.e.
\code{rowSums(getCounts(CoMeth))}
}

