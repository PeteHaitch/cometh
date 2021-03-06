% Generated by roxygen2 (4.0.1): do not edit by hand
\name{CoMeth}
\alias{CoMeth}
\title{The constructor function for CoMeth objects}
\usage{
CoMeth(sample_names = CharacterList(), methylation_type = CharacterList(),
  counts = DataFrameList(), seqnames = RleList(), pos = DataFrameList(),
  seqinfo = Seqinfo(), strand = RleList(), colData = DataFrame(),
  exptData = SimpleList(), ..., verbose = FALSE)
}
\arguments{
\item{sample_names}{A \code{\link[IRanges]{CharacterList}} containing the
names of the samples. Sample names must be unique.}

\item{methylation_type}{A \code{\link[IRanges]{CharacterList}} storing the
type of methylation loci for these m-tuples. \code{methylation_type} must be
named and the names must match those given in \code{sample_names}. For each
sample the possible values are "CG", "CHG", "CHH" or "CNN" or multiple values
specified as a character vector, e.g. c("CG", "CHG") or c("CHG", "CG") for
"CG/CHG" methylation. The \code{methylation_type} of the resulting
\code{CoMeth} object is the union of this argument, e.g. if
\code{methylation_type = list('sample1' = 'CG', 'sample2' = 'CHH')} then the
\code{methylation_type} of the resulting \code{CoMeth} object is
\code{CG/CHH}; this also gives a warning.}

\item{counts}{A \code{\link[IRanges]{DataFrameList}}. \code{counts} must be
named and the names must match those given in \code{sample_names}. The
entries in each DataFrame corresponds to the number of times that particular
co-methylation pattern (columns) was observed for that particular m-tuple
(rows). Therefore, each DataFrame must have the same number of rows as its
corresponding DataFrame in \code{pos} and have \eqn{2 ^ \code{m}} columns.
All samples must have the same m, that is, all DataFrames
in \code{counts} must have the same number of columns.
The column names of each DataFrame must match those given by
\code{.make_m_tuple_names(m)}.}

\item{seqnames}{a \code{\link[IRanges]{RleList}}. \code{seqnames} must be
named and the names must match those given in \code{sample_names}. Each Rle
is the \code{\link[GenomicRanges]{seqnames}} of each m-tuple for that sample.}

\item{pos}{A \code{\link[IRanges]{DataFrameList}}. \code{pos} must be named
and the names must match those given in \code{sample_names}. The columns of
each DataFrame must be: \code{seqnames}, \code{pos1}, ..., \code{posm},
where, for example, \code{posm} is \code{pos3} if \code{m} = 3.
Therefore, the number of columns of each DataFrame is m, the size of the
m-tuples, and the number of rows is equal to the number of m-tuples for that
particular sample. All samples must have the same m, that is, all DataFrames
in \code{pos} must have the same number of columns.}

\item{seqinfo}{A \code{\link[GenomicRanges]{Seqinfo}} object containing
information about the reference genome of the samples. Only a single value is
accepted because all samples must be mapped against the same reference
genome. However, multiple genomes per sample are allowed in order to
accommodate a spiked-in unmethylated genome (normally lambda phage), which is
a common step in a bisulfite-sequencing protocol. \strong{NOTE}:
\code{seqinfo} is a required argument in order to construct a \code{CoMeth}
object, whereas for many Bioconductor objects it is an optional argument.}

\item{strand}{An optional \code{\link[IRanges]{RleList}} object containing
the strand information of each m-tuple. \strong{WARNING}: If \code{strand} is
not supplied, all m-tuples in the resulting \code{CoMeth} object will
have their strand set to \code{*} to signify that the strand is unknown or
irrelevant.
\strong{WARNING}: m-tuples will not be combined across samples if they are on
different strands.}

\item{colData}{An optional \code{\link[IRanges]{DataFrame}} describing the
samples. Row names must match the \code{sample_names}, otherwise an error is
returned.
\code{sample_names}.}

\item{exptData}{An optional \code{\link[IRanges]{SimpleList}} of arbitrary
content describing the overall experiment.}

\item{...}{Additional arguments passed to the internal call to the
\code{\link{MTuples}} constructor. See \code{\link{MTuples}} for details on
what these might be.}

\item{verbose}{A \code{logical(1)} indicating whether messages about data
coercion during construction should be printed.}
}
\value{
A \code{CoMeth1} (if \code{m} \eqn{= 1}), \code{CoMeth2} (if
\code{m} \eqn{= 2}) or \code{CoMeth3Plus} (if \code{m} \eqn{>= 3}) object.
All these are concrete subclasses of the VIRTUAL \code{\link{CoMeth}} class.


cat("TODO")
}
\description{
This is the general function to construct a \code{CoMeth} object. Users will generally use the \code{\link{read.comethylation}} function to construct a \code{CoMeth} object. The constructor allows for multiple samples in the same object by passing most arguments as a *List where each element of the *List corresponds to the arguments for a given sample. These multiple samples will then be appropriately combined into a single \code{CoMeth} object.
}
\seealso{
\code{\link{read.comethylation}} for a function to read in the
\code{.tsv} output file of \code{comethylation} and construct a \code{CoMeth}
object.

\code{\link[GenomicRanges]{SummarizedExperiment}} for the class that
\code{CoMeth} extends.
}

