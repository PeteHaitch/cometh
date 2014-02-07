\name{make_m_tuple_names}
\alias{make_m_tuple_names}
\title{Make m-tuple names}
\usage{
make_m_tuple_names(m)
}
\arguments{
  \item{m}{The size of the m-tuple. Must be an int.}
}
\value{
A character vector
}
\description{
This helper function constructs m-tuple names in the
correct (alphabetical) order for a given value of m
}
\examples{
make_m_tuple_names(1L)
make_m_tuple_names(2L)
make_m_tuple_names(3L)
}
\keyword{helper}

