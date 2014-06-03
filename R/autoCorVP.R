# rd_cometh = rowData(cometh). unstrand(rowData(cometh)) if ignore_strand = TRUE
# ipd = IPD-vector
.autoCorVP <- function(rd_cometh, ipd) {
  
  # Sort rd_cometh without regard to strand.
  # rd_cometh_idx allows us to map back to the origin rd_cometh object.
  # z is the sorted version of rd_cometh
  rd_cometh_idx <- order(unstrand(rd_cometh))
  z <- rd_cometh[rd_cometh_idx]
  
  # p is a list of candidate pairs.
  # p$x (resp. p$y) indexes the first (resp. last) element of each pair. 
  # where the index is with respect to z.
  p <- .makeAutoCorVP(s = start(z), ipd = ipd)
  
  # Determine which pairs are "valid pairs", i.e. satisfy 
  # the seqlevel-feature-strand (sfsc) constraint.
  same_seqlevel <- seqnames(z)[p$x] == seqnames(z)[p$y]
  same_strand <- strand(z[p$x]) == strand(z[p$y])
  same_feature_class <- mcols(z)$ifc[p$x] == mcols(z)$ifc[p$y]
  # z_vp_idx is the index of valid pairs in (p$x, p$y).
  z_vp_idx <- as.vector(same_seqlevel & same_strand & same_feature_class)
  # Remove those pairs that don't satisfy the sfsc contraint.
  x <- p$x[z_vp_idx]
  y <- p$y[z_vp_idx]
  
  # xx and yy translate x and y from z to rd_cometh
  # THIS IS THE INDEX WE CARE ABOUT!
  # (xx, yy) are on the same space as the original cometh object whereas 
  # (x, y) are on the space of the z object.
  xx <- rd_cometh_idx[x]
  yy <- rd_cometh_idx[y]
  
  # Construct the pairs from rd_cometh.
  rd_cometh_vp <- GRanges(seqnames(rd_cometh)[xx], 
                          IRanges(start(rd_cometh)[xx], end(rd_cometh)[yy]),
                          strand(rd_cometh)[xx],
                          ifc = mcols(rd_cometh)[['ifc']][xx])
  # Sanity check: check that there are no duplicates.
  # NOTE: Using any(duplicated()) because there is no anyDuplicated-method 
  # defined for GRanges objects.
  stopifnot(!any(duplicated(rd_cometh_vp)))
  
  # Construct all feature-strand combinations, fsc.
  fsc <- expand.grid(ifc = unique(mcols(z)[['ifc']]),
                     strand = unique(strand(z)))
  
  # Now, create an index (vp2fsipdc) for every valid pair in rd_cometh.
  # vp2fsipdc maps each valid  pair to its feature-strand-IPD combination 
  # (fsipdc). Using "artificial" GRanges for fast matching (via findOverlaps) of 
  # each pair to a feature-strand-IPD combination.
  # They are "artificial" because the ranges are identically (IPD, IPD) and the
  # seqnames are the levels of the feature.
  fsipdc_gr <- GRanges(seqnames = Rle(fsc[['ifc']], 
                                      lengths = rep(length(ipd), nrow(fsc))),
                       ranges = IRanges(start = rep(ipd, nrow(fsc)), width = 1),
                       strand = Rle(fsc[['strand']], 
                                    lengths = rep(length(ipd), nrow(fsc))))
  # rd_cometh_vp_fsipdc = rd_cometh valid pairs' feature-strand-IPD combination.
  rd_cometh_vp_fsipdc_gr <- GRanges(seqnames = mcols(rd_cometh_vp)[['ifc']],
                                    ranges = IRanges(start = 
                                                       width(rd_cometh_vp) - 1,
                                                     width = 1),
                                    strand = strand(rd_cometh_vp))
  vp2fsipdc <- subjectHits(findOverlaps(rd_cometh_vp_fsipdc_gr, fsipdc_gr))
  
  return(list(xx = xx, yy = yy, vp2fsipdc = vp2fsipdc))
}