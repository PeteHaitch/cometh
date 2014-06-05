# rd_cometh = rowData(cometh). unstrand(rowData(cometh)) if ignore_strand = TRUE
# ipd = IPD-vector
.autoCorVP <- function(rd_cometh, ipd) {
  
  # p is a list of candidate pairs.
  # p$xx (resp. p$yy) indexes the first (resp. last) element of each pair. 
  # where the index is with respect to rd_cometh.
  p <- .makeAutoCorVP(s = start(rd_cometh), ipd = ipd)
  
  # Determine which pairs are "valid pairs", i.e. satisfy 
  # the seqlevel-feature-strand (sfsc) constraint.
  same_seqlevel <- seqnames(rd_cometh)[p$x] == seqnames(rd_cometh)[p$y]
  same_strand <- strand(rd_cometh[p$x]) == strand(rd_cometh[p$y])
  same_feature_class <- mcols(rd_cometh)$ifc[p$x] == mcols(rd_cometh)$ifc[p$y]
  # z_vp_idx is the index of valid pairs in (p$x, p$y).
  rd_cometh_vp_idx <- as.vector(same_seqlevel & same_strand & 
                                  same_feature_class)
  # Remove those pairs that don't satisfy the sfsc contraint.
  xx <- p$xx[rd_cometh_vp_idx]
  yy <- p$yy[rd_cometh_vp_idx]
  
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
  fsc <- expand.grid(ifc = unique(mcols(rd_cometh)[['ifc']]),
                     strand = unique(strand(rd_cometh)))
  
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
  
  # fsipdc_df maps vp2fsipdc to the actual fsipdc.
  fsipdc_df <- DataFrame(idx = seq_along(fsipdc_gr), 
                         feature = seqnames(fsipdc_gr), 
                         strand = strand(fsipdc_gr), 
                         ipd = start(fsipdc_gr))
  # Only include those fsipdc that are actually observed.
  fsipdc_df <- fsipdc_df[sort(unique(vp2fsipdc)), ]
  
  return(list(xx = xx, yy = yy, vp2fsipdc = vp2fsipdc, fsipdc_df = fsipdc_df))
}