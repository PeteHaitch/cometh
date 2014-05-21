#' Compute correlations of pairs of beta-values.
#' 
#' Given a \code{\link{CoMeth1}} object, \code{betaCor} computes the 
#' correlations of pairs of beta-values. Correlations are computed per-sample 
#' and are stratified by the \code{ipd} and \code{features} of the pairs of 
#' methylation loci.
#' 
#' @param cometh A CoMeth1 object.
#' @param features The names of features used to stratify the correlations. The
#' default, \code{features = NULL} means that no features are used. Otherwise, 
#' \code{features} must be metadata of the m-tuples in the \code{BetaVals} 
#' object, that is, one-or-more-of \code{colnames(mcols(rowData(BetaVals)))}.
#' @param methylation_loci_set A \code{\link{MethylationLociSet}} object 
#' containing the positions of all methylation loci in the sample. This is 
#' usually the set of all methylation loci in the reference genome.
#' @param ipd An integer vector of IPDs. Correlations will only be computed 
#' using pairs with an IPD found in \code{ipd}. The default value, 
#' \code{sort(unique(width(beta_vals)))}, uses all observed IPDs. 
#' Setting \code{ipd = 0} will removes the IPD strata (see examples).
#' @param nil The number of intervening loci allowed in each pair. The default, 
#' \code{nil = 0}, means that only pairs made of adjacent loci are considered. 
#' To use all pairs, set \code{nil = Inf} but \emph{this will take a very long 
#' time (and probably crash) because there are so many combinations}.
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
betaCor <- function(cometh, methylation_loci_set, features = NULL, 
                    ipd = sort(unique(width(beta_vals))), nil = 0, 
                    method = c("pearson", "kendall", "spearman")){
  if (!is(cometh, "CoMeth1")){
    stop(sQuote("cometh"), " must be a ", sQuote("CoMeth1"), " object.")
  }
  
}