# rd_cometh = sort(rowData(cometh)). sort(unstrand(rowData(cometh))) if ignore_strand = TRUE.
# mls = sort(mls). sort(unstrand(rowD)) if ignore_strand = TRUE.
# ipd = IPD-vector
.generalNILVP <- function(rd_cometh, mls, nil) {
  
  # Both rd_cometh and mls should already be sorted (and unstranded, if reqrd.)
  
  # x and y index the first and second methlyation loci in each pair, 
  # respectively, constructed from mls.
  # NOTE: It is _much_ faster to extract-then-subset than subset-then-extract,
  # e.g. start(cometh)[y] is much faster than start(cometh[y]).
  # This is because subsetting the large CoMeth1 object is slow.
  n <- length(mls)
  x <- seq.int(from = 1, to = n - nil - 1, by = 1)
  y <- seq.int(from = nil + 2, to = n, by = 1)
  ipd0 <- start(mls)[y] - start(mls)[x]
  same_seqlevel <- seqnames(mls)[x] == seqnames(mls)[y]
  same_strand <- strand(mls[x]) == strand(mls[y])
  same_feature_class <- mcols(mls)[['ifc']][x] == mcols(mls)[['ifc']][y]
  # mls_vp is the index of valid pairs in mls (no NIL constraint, yet).
  mls_vp_idx <- as.vector(same_seqlevel & same_strand & same_feature_class)
  # Remove those pairs that don't satisfy the sfsc contraint.
  x <- x[mls_vp_idx]
  y <- y[mls_vp_idx]
  
  # mls_vp_gr is a GRanges object of the valid pairs in mls 
  # (no NIL constraint, yet).
  mls_vp_gr <- GRanges(seqnames = seqnames(mls)[x], 
                       ranges = IRanges(start = start(mls)[x],
                                        end = start(mls)[y]),
                       strand = strand(mls)[x],
                       ifc = mcols(mls)[['ifc']][x])
  
  # This is just a sanity check and isn't required.
  # nil0 is the actual number of intervening loci for each pair.
  nil0 <- Rle(countOverlaps(query = mls_vp_gr, subject = mls, type = 'any')) - 2L
  stopifnot(.zero_range(runValue(nil0)) || as.vector(nil0[1] != nil))
  
  # Find which valid pairs can be constructed from cometh.
  # Do this by matching rd_cometh to start and end of mls_vp_gr.
  # Then, remove those loci where only the start or end matches.
  # xx and yy index the first and second methlyation loci in each pair, 
  # respectively, constructed from rd_cometh
  s <- findOverlaps(query = rd_cometh, subject = mls_vp_gr, type = 'start')
  e <- findOverlaps(query = rd_cometh, subject = mls_vp_gr, type = 'end') 
  xx <- queryHits(s)[!is.na(match(subjectHits(s), subjectHits(e)))]
  yy <- queryHits(e)[!is.na(match(subjectHits(e), subjectHits(s)))]
  
  # Sanity check
  stopifnot(identical(seqnames(rd_cometh)[xx], seqnames(rd_cometh)[yy]))
  stopifnot(identical(strand(rd_cometh)[xx], strand(rd_cometh)[yy]))
  
  # Construct valid pairs from rd_cometh.
  rd_cometh_vp <- GRanges(seqnames = seqnames(rd_cometh)[xx], 
                          ranges = IRanges(start(rd_cometh)[xx], 
                                           end(rd_cometh)[yy]),
                          strand(rd_cometh)[xx])
  # Add ifc to rd_cometh_vp
  ol <- findOverlaps(rd_cometh_vp, mls_vp_gr, type = 'equal')
  mcols(rd_cometh_vp)[['ifc']] <- mcols(mls_vp_gr)[['ifc']][subjectHits(ol)]
  # Sanity check: Check that all pairs construced from rd_cometh actually exist 
  # in mls_vp_gr.
  sanity_check <- as.logical(countQueryHits(ol))
  stopifnot(all(sanity_check)) 
  
  # Construct all feature-strand combinations, fsc.
  fsc <- expand.grid(ifc = unique(mcols(mls)[['ifc']]), 
                     strand = unique(strand(mls)))
  
  # Compute the observed IPDs for each feature-strand combination.
  # Using "artificial" GRanges for fast matching (via findOverlaps) of each pair 
  # to a feature-strand cobination.
  # The are "artificial" because the ranges are identically (1, 1) and the 
  # seqnames are levels of the feature.
  fsc_gr <- GRanges(seqnames = fsc[['ifc']], 
                    ranges = IRanges(start = 1, width = 1),
                    strand = fsc[['strand']])
  # rd_cometh_vp_fsc = rd_cometh valid pairs' feature-strand combination.
  rd_cometh_vp_fsc_gr <- GRanges(seqnames = mcols(rd_cometh_vp)[['ifc']], 
                                 ranges = IRanges(start = 1, width = 1),
                                 strand = strand(rd_cometh_vp))
  # vp2fsc = matches each vp to its correponding fsc.
  rd_cometh_vp2fsc <- subjectHits(findOverlaps(rd_cometh_vp_fsc_gr, fsc_gr))
  # Compute the observed IPDs, i.e. the (unique) widths of the rd_cometh valid 
  # pairs, for each fsc.
  ipd <- tapply(X = width(rd_cometh_vp), INDEX = rd_cometh_vp2fsc, 
                FUN = function(x){
                  unique(sort(x))
                }, simplify = FALSE)
  
  # Now, create an index (vp2fsipdc) for every valid pair in rd_cometh.
  # vp2fsipdc maps each valid  pair to its feature-strand-IPD combination 
  # (fsipdc). Using "artificial" GRanges for fast matching (via findOverlaps) of 
  # each pair to a feature-strand-IPD combination.
  # They are "artificial" because the ranges are identically (IPD, IPD) and the
  # seqnames are the levels of the feature.
  fsipdc_gr <- GRanges(seqnames = Rle(fsc[['ifc']], 
                                      lengths = sapply(ipd, length)),
                       ranges = IRanges(start = unlist(ipd), width = 1),
                       strand = Rle(fsc[['strand']], 
                                    lengths = sapply(ipd, length)))
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