## TODO: Add parallel-support either via parallel or BiocParallel. Should 
## parallelisation be by sample or by strata?
## TODO:  Feature should be a GRanges object; every loci in mls will be checked 
## whether it is in feature and this will be recorded as metadata of mls.
## TODO: Add option to specify the minimum number of observations per strata in
## order to compute a correlation for that strata.
## TODO: A "special value" of nil should mean that all pairs are used, 
## regardless of NIL, i.e. a traditional autocorrelation of the beta-values.

#' Compute correlations of pairs of beta-values.
#' 
#' Given a \code{\link{CoMeth1}} object, \code{betaCor} computes the 
#' correlations of pairs of beta-values. Correlations are computed per-sample 
#' and are stratified by the \code{ipd} and \code{features} of the pairs of 
#' methylation loci.
#' 
#' @param cometh A CoMeth1 object.
#' @param mls A \code{\link{MethylationLociSet}} object 
#' containing the positions of all methylation loci in the sample. This is 
#' usually the set of all methylation loci in the reference genome.
#' If a \code{\link{MethylationLociSet}} is unavailable, then set to \code{NA}, 
#' in which case \code{betaCor} effectively assumes that \code{cometh} contains 
#' all potential methylation loci (this is not recommended).
#' @param ipd An integer vector of IPDs. Correlations will only be computed 
#' using pairs with an IPD found in \code{ipd}. The default value, 
#' \code{NULL}, uses all observed IPDs. Setting \code{ipd = 0} will removes the 
#' IPD strata (see examples).
#' @param nil The number of intervening loci allowed in each pair. The default, 
#' \code{nil = 0}, means that only pairs made of adjacent loci are considered. 
#' To use all pairs, set \code{nil = Inf} but \emph{this will take a very long 
#' time (and probably crash) because there are so many combinations}.
#' @param feature A character. The name of the feature used to stratify the 
#' correlations. The default, \code{feature = NA} means that no feature is 
#' used. Otherwise, \code{feature} must be metadata of the m-tuples in the 
#' \code{cometh} object, that is, one of \code{colnames(mcols(cometh))}.
#' @param ignore_strand logical. Should the strand of the methylation loci be 
#' ignored (default: \code{FALSE})?
#' @param method A character string indicating which correlation coefficient 
#' is to be computed. One of "pearson" (default), "kendall", or "spearman", 
#' can be abbreviated. \emph{Currently, only "pearson" and "spearman" are 
#' supported because "kendall" is very slow (see Note in 
#' \code{\link[stats]{cor}}).}
#'
#' @note Correlations are computed using \code{\link[stats]{cor}}, 
#' with \code{use = "na.or.complete".}
#' 
#' @return A \code{data.frame}. While the exact format depends on the 
#' parameters used, the first column is the correlation, the second column is 
#' the sample name and the remaining columns are \code{IPD} and 
#' \code{features}, if relevant. Therefore, each row corresponds to the 
#' correlation for a particular "sample-IPD-features" combination.
#' 
betaCor <- function(cometh, mls, ipd, nil = 0L, feature = NA, 
                    ignore_strand = FALSE,
                    method = c("pearson", "kendall", "spearman")) {
  
  ## TODO: More argument checks
  
  # Check that all required arguments are not missing
  if (missing(cometh)) {
    stop(sQuote("cometh"), " missing.\nPlease see the help page for betaCor, ",
         "which can accessed by typing ", sQuote("?betaCor"), " at the R ",
         "prompt, for further details of this argument.")
  }
  if (missing(mls)) {
    stop(sQuote("mls"), " missing.\nPlease see the help page for betaCor, ", 
         "which can accessed by typing ", sQuote("?betaCor"), " at the R ", 
         "prompt, for further details of this argument.")
  }
  if (missing(ipd)) {
    message(sQuote('ipd'), " not specified, so will use all IPDs.")
  }
  
  # Check that all required arguments are of the correct class 
  if (!is(cometh, "CoMeth1")) {
    stop(sQuote("cometh"), " must be a ", sQuote("CoMeth1"), " object.")
  }
  if (!is(mls, "MethylationLociSet") & suppressWarnings(!is.na(mls))) {
    stop(sQuote("mls"), " must be a ", sQuote("MethylationLociSet"), " object ", 
         "or set as ", sQuote("NA"), " if not available.")
  }
  
  ## TODO: Remove? Probably not necessary.
  # Check that cometh contains the bare minimum number of methylation loci
  n <- nrow(cometh)
  if (n < 3){
    stop("Require at least 3 methylation loci to compute ", sQuote('betaCor'))
  }
  # Check that cometh and mls are defined for the same reference genome
  seqinfo <- try(merge(seqinfo(cometh), seqinfo(mls)), silent = TRUE)
  if (is(seqinfo, "try-error")){
    stop(sQuote("cometh"), " and ", sQuote("mls"), " have incompatible ", 
         sQuote("seqinfo"), ".")
  }
  # Check that cometh doesn't contain any loci not found in mls
  if (isTRUE(any(countOverlaps(cometh, mls) == 0L))) {
    stop("All loci in ", sQuote("cometh"), " must also be present in ", 
         sQuote("mls"), ".")
  }
  
  # Check all non-required arguments
  if (!missing(ipd)) {
    if (!is(ipd, "integer")){
      stop(sQuote('ipd'), " must be a vector of type ", sQuote('integer'), ".")
    }
  }
  if (nil != 0L){
    stop("Sorry, currently only supports ", sQuote("nil = 0"), ".")
  }
  if (all.equal(nil, floor(nil))){
    stop(sQuote("nil"), " must be an integer.")
  }
  if (!(feature %in% names(mcols(cometh)))) {
    stop(sQuote("feature"), " must match one of ", 
         sQuote("colnames(mcols(cometh))"), ".")
  }
  method <- match.arg(method, choices = c("pearson", "kendall", "spearman"), 
                      several.ok = FALSE)
  if (identical(method, "spearman")) {
    stop("Sorry, ", sQuote("method = spearman"), " is not yet supported.")
  }
  if (!isTRUEorFALSE(ignore_strand)){
    stop(sQuote('ignore_strand'), " must be ", sQuote("TRUE"), " or ", 
         sQuote("FALSE"), ".")
  }
  
  ## TODO: Dispatch on different methods depending on the value of NIL.
  ## See above for discussion of methods appropriate for different NIL-values.  
  
  ## TODO: Include methylation_type, NIL, IPD, feature, seqinfo, method etc. 
  ## in output. Probably want a BetaCor class for this job.
  
}

