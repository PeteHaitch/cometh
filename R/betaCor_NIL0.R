## DEPRECATED: Deprecated in favour of betaCor_NILGeneral.R.
## TODO: Add some sort of parallelisation
#' Compute betaCor when NIL = 0.
#' 
#' This function should only be called by betaCor and not directly by the user.
#' @keywords internal
#' 
.betaCorNIL0 <- function(cometh, mls, ipd, ignore_strand){
  # Construct the "valid pairs" of methylation loci.
  # A pair is valid if both methylation loci are on the same seqlevel and 
  # strand, as well as within the same feature class and satisfying the NIL
  # constraint.
  # NOTE: Being in the same feature class is not the same as being in the 
  # same feature. E.g. Two CpGs in two distinct CGIs are both in the same 
  # feature class (CGI) but not the same feature.
  ## TODO: What if there is no feature?
  ## TODO: What if cometh or mls contain all three strands (+, - and *)
  if (ignore_strand){
    cometh <- unstrand(cometh)
  }
  cometh <- sort(cometh)
  # x, y index the first and second methlyation loci in each pair, respectively.
  # NOTE: It is _much_ faster to extract-then-subset than subset-then-extract,
  # e.g. start(cometh)[y] is much faster than start(cometh[y]).
  # This is because subsetting the large CoMeth1 object is slow.
  n <- nrow(cometh)
  x <- seq.int(from = 1, to = n - 1, by = 1)
  y <- seq.int(from = 2, to = n, by = 1)
  ipd0 <- start(cometh)[y] - start(cometh)[x]
  same_chrom <- seqnames(cometh)[x] == seqnames(cometh)[y]
  same_strand <- strand(cometh[x]) == strand(cometh[y])
  same_feature_class <- (mcols(cometh)[[feature]][x] == 
                           mcols(cometh)[[feature]][y])
  # vp is the index of valid pairs (no NIL constraint, yet).
  vp <- same_chrom & same_strand & same_feature_class
  # vp_gr is a GRanges object of the valid pairs (no NIL constraint, yet).
  vp_gr <- GRanges(seqnames = seqnames(cometh[x])[vp], 
                   ranges = IRanges(start = start(cometh)[x[as.vector(vp)]],
                                    end = start(cometh)[y[as.vector(vp)]]),
                   strand = strand(cometh)[x[as.vector(vp)]],
                   ifc = mcols(cometh)[[feature]][x[as.vector(vp)]])

  # nil0 is the actual number of intervening loci for each pair.
  # -2L to remove the "outer" loci from each pair and thus only count "internal"
  # loci.
  nil0 <- Rle(countOverlaps(query = vp_gr, subject = mls, type = 'any')) - 2L  

  # Refine the valid pairs to include only those pairs that satisfy the NIL
  # constraint.
  vp[vp] <- (nil0 == Rle(nil, length(nil0)))
  vp_gr <- vp_gr[nil0 == Rle(nil, length(nil0))]
  
  # Construct all feature-strand combinations, fsc.
  if (ignore_strand) {
    fsc <- expand.grid(feature = c(TRUE, FALSE), strand = c('+', '-'))
  } else{
    fsc <- expand.grid(feature = c(TRUE, FALSE), strand = c('*'))
  }
  
  # If ipd wasn't defined, compute the possible IPDs for each feature-strand 
  # combination.
  # Using artificial GRanges for fast matching (using findOverlaps) of each pair 
  # to a feature-strand cobination.
  # The are artificial because the ranges are identically (1, 1) and the 
  # seqnames are TRUE/FALSE.
  ## TODO: Probably a bit convoluted; use a simpler way to find matches.
  if (missing(ipd)) {
    fsc_gr <- GRanges(seqnames = fsc$feature, 
                      ranges = IRanges(start = 1, width = 1),
                      strand = fsc$strand)
    # vp_fsc = valid pairs' feature-strand combination.
    vp_fsc_gr <- GRanges(seqnames = mcols(vp_gr)$ifc, 
                         ranges = IRanges(start = 1, width = 1),
                         strand = strand(vp_gr))
    # vp2fsc = matches each vp to its correponding fsc.
    vp2fsc <- subjectHits(findOverlaps(vp_fsc_gr, fsc_gr))
    # Find the observed IPDs, i.e. the (unique) widths of the valid pairs, 
    # for each fsc.
    ipd <- tapply(X = width(vp_gr), INDEX = vp2fsc, FUN = function(x){
      unique(sort(x))
    }, simplify = FALSE)
  } else{
    ipd <- lapply(seq_len(nrow(fsc)), function(x){ipd})
  }
  
  # Now, create an index for every valid pair (vp2fsipdc) mapping each valid 
  # pair to its feature-strand-IPD combination, fsipdc.
  # Using artificial GRanges for fast matching (using findOverlaps) of each pair 
  # to a feature-strand combination
  # The are artificial because the ranges are identically (IPD, IPD) and the
  # seqnames are TRUE/FALSE.
  fsipdc_gr <- GRanges(seqnames = Rle(fsc$feature, lengths = sapply(ipd, length)),
                       ranges = IRanges(start = unlist(ipd), width = 1),
                       strand = Rle(fsc$strand, lengths = sapply(ipd, length)))
  # vp_fsipdc = valid pairs' feature-strand-IPD combination.
  vp_fsipdc_gr <- GRanges(seqnames = mcols(vp_gr)$ifc,
                          ranges = IRanges(start = width(vp_gr), width = 1),
                          strand = strand(vp_gr))
  vp2fsipdc <- subjectHits(findOverlaps(vp_fsipdc_gr, fsipdc_gr))
  
  # Now, for each sample, compute correlations of beta values for all pairs
  # with the same vp2fsipdc.
  # Basically, cor(beta_x, beta_y) ~ sample + vp2fsipdc
  # This part could be parallel-ised (by sample or vp2fsipdc?)
  ## TODO: Is by the right function to use? It has some overhead.
  lapply(sampleNames(cometh), function(sn, beta, x_vp, y_vp, vp2fsipdc, 
                                       method) {
    beta <- cbind(beta[x_vp, sn , drop = TRUE], beta[y_vp, sn , drop = TRUE])
    val <- by(data = beta, INDICES = vp2fsipdc, 
              FUN = function(beta, use, method) {
                cor(x = beta[, 1], y = beta[, 2], use = use, method = method)
              }, use = "na.or.complete", method = method, simplify = TRUE)
    return(val)
  }, beta = assay(cometh, "beta"), x_vp = x[as.vector(vp)], 
  y_vp = y[as.vector(vp)], vp2fsipdc = vp2fsipdc, method = method)
}