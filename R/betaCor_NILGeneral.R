## TODO: Add some sort of parallelisation
## TODO: Should (am) I demean(-ing) the beta values before computing cor?

#' Compute betaCor for general NIL.
#' 
#' This function should only be called by betaCor and not directly by the user.
#' @param cometh Passed from betaCor.R
#' @param mls Passed form betaCor.R
#' @param ipd Passed from betaCor.R.
#' @param nil Passed from betaCor.R.
#' @param ignore_strand Passed from betaCor.R.
#' @keywords internal
#' 
.betaCorNILGeneral <- function(cometh, mls, ipd, nil, ignore_strand){
  
  # First, sort both cometh and mls. 
  # This ensures that the first loci in each pair comes before the second.
  cometh <- sort(cometh)
  mls <- sort(mls)
  
  # Construct the "valid pairs" of methylation loci.
  # A pair is valid if both methylation loci are on the same seqlevel and 
  # strand, as well as within the same feature class and satisfying the NIL
  # constraint.
  # We construct valid pairs from mls rather than cometh; this is in contrast
  # to .betaCorNIL0.
  # NOTE: Being in the same feature class is not the same as being in the 
  # same feature. E.g. Two CpGs in two distinct CGIs are both in the same 
  # feature class (CGI) but not the same feature.
  ## TODO: What if there is no feature?
  ## TODO: What if cometh or mls contain all three strands (+, - and *)
  if (ignore_strand){
    mls <- unstrand(mls)
  }

  # x and y index the first and second methlyation loci in each pair, 
  # respectively, constructed from mls.
  # NOTE: It is _much_ faster to extract-then-subset than subset-then-extract,
  # e.g. start(cometh)[y] is much faster than start(cometh[y]).
  # This is because subsetting the large CoMeth1 object is slow.
  n <- length(mls)
  x <- seq.int(from = 1, to = n - nil - 1, by = 1)
  y <- seq.int(from = nil + 2, to = n, by = 1)
  ipd0 <- start(mls)[y] - start(mls)[x]
  same_chrom <- seqnames(mls)[x] == seqnames(mls)[y]
  same_strand <- strand(mls[x]) == strand(mls[y])
  # same_feature_class could be generalised from 2 states to 4 states, i.e.,
  # from TRUE, FALSE 
  # to BOTH_IN, BOTH_OUT, FIRST_IN_SECOND_OUT, FIRST_OUT_SECOND_IN
  # Currently throwing away all pairs with either FIRST_IN_SECOND_OUT or 
  # FIRST_OUT_SECOND_IN.
  same_feature_class <- mcols(mls)[['ifc']][x] == mcols(mls)[['ifc']][y]
  # mls_vp is the index of valid pairs in mls (no NIL constraint, yet).
  mls_vp <- same_chrom & same_strand & same_feature_class
  # mls_vp_gr is a GRanges object of the valid pairs in mls 
  # (no NIL constraint, yet).
  mls_vp_gr <- GRanges(seqnames = seqnames(mls[x])[mls_vp], 
                   ranges = IRanges(start = start(mls)[x[as.vector(mls_vp)]],
                                    end = start(mls)[y[as.vector(mls_vp)]]),
                   strand = strand(mls)[x[as.vector(mls_vp)]],
                   ifc = mcols(mls)[['ifc']][x[as.vector(mls_vp)]])
  
  # This is just a sanity check and isn't required.
  # nil0 is the actual number of intervening loci for each pair.
  nil0 <- Rle(countOverlaps(query = mls_vp_gr, subject = mls, type = 'any')) - 2L
  stopifnot(.zero_range(runValue(nil0)) || as.vector(nil0[1] != nil))
  
  # Find which valid pairs can be constructed from cometh.
  # Do this by matching cometh to start and end of mls_vp_gr.
  # Then, remove those loci where only the start or end matches.
  # xx and yy index the first and second methlyation loci in each pair, 
  # respectively, constructed from cometh.
  s <- findOverlaps(query = cometh, subject = mls_vp_gr, type = 'start')
  e <- findOverlaps(query = cometh, subject = mls_vp_gr, type = 'end') 
  xx <- queryHits(s)[!is.na(match(subjectHits(s), subjectHits(e)))]
  yy <- queryHits(e)[!is.na(match(subjectHits(e), subjectHits(s)))]
    
  stopifnot(identical(seqnames(cometh)[xx], seqnames(cometh)[yy]))
  stopifnot(identical(strand(cometh)[xx], strand(cometh)[yy]))
  
  # Construct valid pairs from cometh.
  cometh_vp_gr <- GRanges(seqnames = seqnames(cometh)[xx],
                          ranges = IRanges(start = start(cometh)[xx], 
                                           end = end(cometh[yy])),
                          strand = strand(cometh)[xx])
  ol <- findOverlaps(cometh_vp_gr, mls_vp_gr, type = 'equal')
  
  # This is a sanity check and isn't strictly required.
  # Check that all pairs construced from cometh actually exist in vp_gr.
  sanity_check <- as.logical(countQueryHits(ol))
  stopifnot(all(sanity_check))  
  
  # Add ifc to cometh_vp_gr
  mcols(cometh_vp_gr)[['ifc']] <- mcols(mls_vp_gr)[['ifc']][subjectHits(ol)]
  
  # Construct all feature-strand combinations, fsc.
  fsc <- expand.grid(ifc = unique(mcols(mls)$ifc), strand = unique(strand(mls)))
  
  # If ipd wasn't defined, compute the possible IPDs for each feature-strand 
  # combination.
  # Using artificial GRanges for fast matching (using findOverlaps) of each pair 
  # to a feature-strand cobination.
  # The are artificial because the ranges are identically (1, 1) and the 
  # seqnames are TRUE/FALSE.
  ## TODO: Probably a bit convoluted; use a simpler way to find matches.
  if (missing(ipd)) {
    fsc_gr <- GRanges(seqnames = fsc[['ifc']], 
                      ranges = IRanges(start = 1, width = 1),
                      strand = fsc$strand)
    # cometh_vp_fsc = cometh valid pairs' feature-strand combination.
    cometh_vp_fsc_gr <- GRanges(seqnames = mcols(cometh_vp_gr)[['ifc']], 
                         ranges = IRanges(start = 1, width = 1),
                         strand = strand(cometh_vp_gr))
    # vp2fsc = matches each vp to its correponding fsc.
    cometh_vp2fsc <- subjectHits(findOverlaps(cometh_vp_fsc_gr, fsc_gr))
    # Find the observed IPDs, i.e. the (unique) widths of the valid pairs, 
    # for each fsc.
    ipd <- tapply(X = width(cometh_vp_gr), INDEX = cometh_vp2fsc, 
                  FUN = function(x){
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
  fsipdc_gr <- GRanges(seqnames = Rle(fsc[['ifc']], 
                                      lengths = sapply(ipd, length)),
                       ranges = IRanges(start = unlist(ipd), width = 1),
                       strand = Rle(fsc$strand, lengths = sapply(ipd, length)))
  # cometh_vp_fsipdc = cometh valid pairs' feature-strand-IPD combination.
  cometh_vp_fsipdc_gr <- GRanges(seqnames = mcols(cometh_vp_gr)[['ifc']],
                                 ranges = IRanges(start = width(cometh_vp_gr), 
                                                  width = 1),
                                 strand = strand(cometh_vp_gr))
  vp2fsipdc <- subjectHits(findOverlaps(cometh_vp_fsipdc_gr, fsipdc_gr))
  
  # Now, for each sample, compute correlations of beta values for all pairs
  # with the same vp2fsipdc.
  # Basically, cor(beta_x, beta_y) ~ sample + vp2fsipdc
  # This part could be parallel-ised (by sample or vp2fsipdc?)
  ## TODO: Is by the right function to use? It has some overhead.
  val <- lapply(sampleNames(cometh), function(sn, beta, xx, yy, vp2fsipdc, 
                                       method) {
    beta <- cbind(beta[xx, sn , drop = TRUE], beta[yy, sn , drop = TRUE])
    val <- by(data = beta, INDICES = vp2fsipdc, 
              FUN = function(beta, use, method) {
                cor(x = beta[, 1], y = beta[, 2], use = use, method = method)
              }, use = "na.or.complete", method = method, simplify = TRUE)
    return(val)
  }, beta = assay(cometh, "beta"), xx = xx, yy = yy, vp2fsipdc = vp2fsipdc, 
  method = method)
  
  df <- DataFrame(sample = Rle(sampleNames(cometh), 
                               rep(length(fsipdc_gr), ncol(cometh))),
                  IFC = rep(seqnames(fsipdc_gr), ncol(cometh)),
                  strand = rep(strand(fsipdc_gr), ncol(cometh)),
                  IPD = rep(start(fsipdc_gr), ncol(cometh)),
                  cor = unlist(val)
                  )
  
  return(df)
}