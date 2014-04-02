### =========================================================================
### MTuple objects
### -------------------------------------------------------------------------

## An MTuple object is a SimpleList comprised of two elements: 'pos' and 'counts'. 
## The 'pos' is a Tuple object, which stores the genomic co-ordinates of the tuples.
## A tuple consists of a seqname, a strand, and the positions.
## All co-ordinates are 1-based.
## For example, chr1:+:(3, 7, 9) is a text-representation of tuple containing the three points (3, 7, 9) on the forward strand of chromosome 1. 
## Tuples are different to ranges (e.g. GRanges or IRanges) in that any point not explicity listed is not included in the tuple (e.g. chr1:+:4, chr1:+:5 or chr1:-:3 in the above example).
## The 'counts' is a SimpleList of matrices, where each matrix storing the counts of a co-methylation pattern for each sample at each tuple. The rows of the matrix are tuples and the columns are samples

## The class definition is adapted from Martin Morgan's reply to my question on Bioc-Devel (https://stat.ethz.ch/pipermail/bioc-devel/2014-February/005242.html) and Michael Lawerence's reply to my question on BioC-Help (https://stat.ethz.ch/pipermail/bioconductor/2014-March/058490.html)

.MTuple <- setClass("MTuple", 
                    representation(
                    exptData="SimpleList",                # overall description
                    pos = "Tuple",                        # tuples and their description
                    counts = "SimpleList"))
