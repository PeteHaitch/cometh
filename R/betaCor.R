### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bestCor: Compute within-sample correlations of beta-values.
###
### betaCor defers to specific methods depending on what parameters it is 
### passed.
### TODO: Add parallel-support either via parallel or BiocParallel. Should 
### parallelisation be by sample or by strata?
### TODO: Add option to specify the minimum number of observations per strata in
### order to compute a correlation for that strata?
### TODO: Extend to between-sample correlations of beta-values. This includes 
### between-sample correlations within each pair of loci and or between pairs of
### loci with the same feature-strand-IPD ombination (fsipdc).
### TODO: Make Spearman correlation the default because data are non-Normal?

#' Compute correlations of pairs of beta-values.
#' 
#' Given a \code{\link{CoMeth1}} object, \code{betaCor} computes the 
#' correlations of pairs of beta-values. Correlations are computed per-sample 
#' and are stratified by the \code{strand}, \code{features} and \code{ipd} of 
#' the pairs of methylation loci.
#' 
#' @param cometh A CoMeth1 object.
#' @param mls A \code{\link{MethylationLociSet}} object 
#' containing the positions of all methylation loci in the sample. This is 
#' usually the set of all methylation loci in the reference genome.
#' If a \code{\link{MethylationLociSet}} is unavailable, then set to \code{NA}, 
#' in which case \code{betaCor} effectively assumes that \code{cometh} contains 
#' all potential methylation loci (this is not recommended).
#' @param ipd An integer vector of IPDs. Correlations will only be computed 
#' using pairs with an IPD found in \code{ipd}. If \code{ipd} is specified, then
#' the number of intervening loci (NIL) is ignore. Only one of \code{ipd} or
#' \code{nil} should be specified.
#' @param nil The number of intervening loci (NIL) allowed in each pair. The 
#' default, \code{nil = 0}, means that only pairs made of adjacent loci are 
#' considered. See the \code{ipd} parameter for information on how to use all 
#' pairs, regardless of NIL. Only one of \code{ipd} or \code{nil} should be 
#' specified.
#' @param feature A \code{GRanges} object. If supplied, the correlations are 
#' stratified by whether the methylation loci overlap a range in the supplied 
#' \code{feature}. For example, correlations might be stratified by whether the
#' loci are from a CGI. The details of this stratification are controlled by the 
#' \code{same_strand} parameter.
#' @param same_feature logical. If supplying a \code{feature}, should pairs of
#' methylation loci be required to be in the same feature (TRUE) or just in the 
#' same feature class (FALSE). Two loci are in the same feature class if, for 
#' example, they are both in CGIs but each loci is in a \emph{different} CGI. 
#' This only affects pairs of loci that are in the feature, i.e., pairs in the 
#' "gaps" betweens features need not be in the same "gap".
#' @param ignore_strand logical. Should the strand of the methylation loci be 
#' ignored (default: \code{FALSE}).
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
betaCor <- function(cometh, mls, nil = 0L, ipd, feature, feature_name,
                    same_feature = FALSE,
                    ignore_strand = FALSE,
                    method = c("pearson", "kendall", "spearman")) {
  
  # Check that all required arguments are not missing.
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
  if ((missing(nil) & missing(ipd)) || (!missing(nil) & !missing(ipd))){
    stop("One, and only one, of ", sQuote("nil"), " and ", sQuote("ipd"), 
         " should be supplied.\nPlease see the help page for betaCor, ", 
         "which can accessed by typing ", sQuote("?betaCor"), " at the R ", 
         "prompt, for further details of these arguments.")
  }
  # Check that all required arguments are of the correct class 
  if (!is(cometh, "CoMeth1")) {
    stop(sQuote("cometh"), " must be a ", sQuote("CoMeth1"), " object.")
  }
  if (!is(mls, "MethylationLociSet") & suppressWarnings(!is.na(mls))) {
    stop(sQuote("mls"), " must be a ", sQuote("MethylationLociSet"), " object ", 
         "or set as ", sQuote("NA"), " if not available.")
  }
  if (missing(ipd)){
    if (length(nil) != 1L || !is.numeric(nil) || (nil != floor(nil))) {
      stop(sQuote("nil"), " must be an integer.")
    }
    nil <- NA_integer_
  } else if (missing(nil)){
    if (!is.numeric(ipd)){
      stop(sQuote('ipd'), " must be an integer vector.")
    }
    ipd <- NA_integer_
  }
  
  # Check that cometh and mls are defined on compatible reference genomes.
  seqinfo <- try(merge(seqinfo(cometh), seqinfo(mls)), silent = TRUE)
  if (is(seqinfo, "try-error")){
    stop(sQuote("cometh"), " and ", sQuote("mls"), " have incompatible ", 
         sQuote("seqinfo"), ".")
  }
  
  # Check that cometh doesn't contain any loci not found in mls.
  if (isTRUE(any(countOverlaps(cometh, mls) == 0L))) {
    stop("All loci in ", sQuote("cometh"), " must also be present in ", 
         sQuote("mls"), ".")
  }
  
  # Check that all non-required arguments, if supplied, are correct.
  # feature
  if (!inherits(feature, "GRanges")){
    stop(sQuote('feature'), " must be a ", sQuote("GRanges"), " object or ", 
         "inherit from the ", sQuote("GRanges"), " class.")
  }
  seqinfo <- try(merge(seqinfo(cometh), seqinfo(feature)), silent = TRUE)
  if (is(seqinfo, "try-error")){
    stop(sQuote("cometh"), " and ", sQuote("feature"), " have incompatible ", 
         sQuote("seqinfo"), ".")
  }
  if (!missing(feature)) {
    if (missing(feature_name)){
      stop("Please supply a ", sQuote('feature_name'), ".")
    }
    if (!isTRUEorFALSE(same_feature)) {
      stop(sQuote('same_feature'), " must be ", sQuote("TRUE"), " or ", 
           sQuote("FALSE"))
    }
  }
  if (!isTRUEorFALSE(ignore_strand)) {
    stop(sQuote('ignore_strand'), " must be ", sQuote("TRUE"), " or ", 
         sQuote("FALSE"))
  }
  method <- match.arg(method, choices = c("pearson", "kendall", "spearman"), 
                      several.ok = FALSE)
  if (identical(method, "spearman")) {
    stop("Sorry, ", sQuote("method = spearman"), " is not yet supported.")
  }
  
  # Remove the strand if ignore_strand = TRUE
  if (ignore_strand){
    cometh <- unstrand(cometh)
    mls <- unstrand(cometh)
  }
  
  # Extract rd_cometh and sort it. 
  # Also sort mls so that it is consistent with cometh.
  rd_cometh <- sort(rowData(cometh))
  mls <- sort(mls)
  
  # Annotate each methylation loci by feature (if feature supplied).
  if (!missing(feature)) {
    mcols(rd_cometh)$ifc <- overlapsAny(rd_cometh, feature)
    mcols(mls)$ifc <- overlapsAny(mls, feature)
  } else {
    mcols(rd_cometh)$ifc <- TRUE
    mcols(mls)$ifc <- TRUE
  }

  # Construct the index of "valid pairs".
  # Use .autoCorVP() if constructing pairs for a given IPD and ignoring NIL.
  # Use .generalNILVP() if contructing pairs for a given NIL and ignoring IPD.
  if (!is.na(ipd)) {
    vp_idx <- .autoCorVP(rd_cometh = rd_cometh, ipd = ipd)
  } else{
    vp_idx <- .generalNILVP(rd_cometh = rd_cometh, mls = mls, nil = nil)
  }
  
  # If same_feature = TRUE, then remove those pairs that are in the same 
  # feature class but not in the same feature.
  # Those pairs not in the feature class are unaffected since we already know 
  # that they cannot be in the same feature because neither of them is in the 
  # same feature class.
  if (same_feature) {
    # Construct candidate valid pairs, cvp.
    cvp <- GRanges(seqnames = seqnames(rd_cometh)[vp_idx[['xx']]],
                   ranges = IRanges(start = start(rd_cometh)[vp_idx[['xx']]],
                                    end = end(rd_cometh)[vp_idx[['yy']]]),
                   strand = strand(rd_cometh)[vp_idx[['xx']]],
                   ifc = mcols(rd_cometh)[['ifc']][vp_idx[['xx']]])
    
    # Remove those pairs not within a single feature (nwasf).
    # The ifc rd_cometh metadata then means that the pair is within a single 
    # feature (TRUE) or not in the feature class (FALSE).
    nwasf <- !overlapsAny(cvp, feature, type = 'within') & 
      (mcols(cvp)[['ifc']] == TRUE)
    fscipdc_df_to_remove <- vp_idx[['fsipdc_df']][['idx']] %in%
      sort(unique(vp_idx[['vp2fsipdc']][!nwasf]))
    vp_idx <- list(xx = vp_idx[['xx']][!nwasf], yy = vp_idx[['yy']][!nwasf], 
                   vp2fsipdc = vp_idx[['vp2fsipdc']][!nwasf], 
                   fsipdc_df = vp_idx[['fsipdc_df']][fscipdc_df_to_remove, ])
  }
  
  # Compute correlations of "valid pairs", stratified by vp2fsipdc.
  ## TODO: Is lapply+by the best way to compute these?
  ## TODO: Could be parallelised; by sample or by index or by both?
  ## TODO: Should warnings produced by cor be suppressed?
  corr <- lapply(sampleNames(cometh), function(sn, beta, xx, yy, vp2fsipdc,
                                              method) {
    beta <- cbind(beta[xx, sn, drop = TRUE], beta[yy, sn, drop = TRUE])
    val <- by(data = beta, INDICES = vp2fsipdc,
              FUN = function(beta, use, method) {
                cor(x = beta[, 1], y = beta[, 2], use = use, method = method)
              }, use = "na.or.complete", method = method, simplify = TRUE)
    return(val)
  }, beta = assay(cometh, "beta"), xx = vp_idx$xx, yy = vp_idx$yy, 
  vp2fsipdc = vp_idx$vp2fsipdc, method = method)
  
  ## TODO: Refine output. Should include methylation_type, NIL, IPD, feature, 
  ## feature_name, seqinfo, method etc. in output. 
  ## Columns of output should include sampleName (only for within-sample 
  ## correlatons), co-ordinates (only for between-sample correlations of 
  ## specific pairs), IFC (if feature supplied), strand, IPD and cor.
  df <- DataFrame(sample_name = Rle(sampleNames(cometh), sapply(corr, length)),
                  in_feature = rep(vp_idx[['fsipdc_df']][['feature']], 
                                   ncol(cometh)),
                  strand = rep(vp_idx[['fsipdc_df']][['strand']], 
                               ncol(cometh)),
                  IPD = rep(vp_idx[['fsipdc_df']][['ipd']], 
                            ncol(cometh)),
                  corr = unlist(corr, use.names = FALSE)
  )
  beta_cor <- new("BetaCor", df, methylation_type = getMethylationType(cometh), 
             NIL = nil, IPD = ipd, feature_name = feature_name, 
             same_feature = same_feature, ignore_strand = ignore_strand, 
             seqinfo = seqinfo(cometh), method = method)

}
