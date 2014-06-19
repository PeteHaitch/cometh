## The case of NIL = Inf is best solved by "traditional" autocorrelation 
## methods, i.e. treat the beta-values as a vector (perhaps suitably padded 
## by NAs) for all positions, not just methylation loci, without an observed 
## beta-value. I would still need a way to combine data from different 
## chromosomes. An existing algorithm will surely do a better job of computing
## fast and accurate autocorrelations in this case than anything I write 
## myself.
## Unfortunately, the genome that I work with (human, mouse) are longer than the
## longest allowable vector in R (2^32 - 1), which is too small to record the 
## every position in these genomes. There is some support for long vectors 
## (2^64 - 1) but it is limited. In any case, this is an incredibly inefficient 
## way of storing beta-values since 99% of sites have beta-value set to NA.
## Computing per-sample-chromosome acfs seems doable (at least, for lag.max no 
## bigger than a few thousand). But there's no way to post-hoc combine these.
## TODO: Write the internal function to compute autocorrelation of beta-values.
## Will use acf (or some other fast autocorrelation function).
## TODO: Read Venables and Ripley for details of acf
## TODO: Add some sort of parallelisation
## TODO: Warn the user that this will use __LOTS__ of memory for large (i.e. 
## mammalian) chromosomes.

#' Compute autocorrelations of beta-values.
#' 
#' This function should only be called by betaCor and not directly by the user.
#' @param cometh Passed from betaCor.R
#' @param max_lag Passed form betaCor.R
#' @param ignore_strand Passed from betaCor.R.
#' @keywords internal
#' 
.betaAutoCor <- function(cometh, max_lag, ignore_strand){
  
  cometh <- sort(cometh)
  
  if (ignore_strand){
    unstrand(cometh)
  }

  # Construct all seqlevels-feature-strand combinations, sfsc
  sfsc <- expand.grid(seqlevels = seqlevels(cometh),
                      ifc = unique(mcols(cometh)$ifc),
                      strand = unique(strand(cometh)))  
  
  # Option 1: acf-based. WON'T WORK/SCALE BASED ON CURRENT (ad hoc) TESTING.
  # Create a LC-vector for each chromosome, where LC = length of chromosome.
  # Record the beta-value or NA, appropriate.
  # Convert the LC-vectors to a multivariate time-series, z.
  # Compute acf(z).
  # Unfortunately, the genome that I work with (human, mouse) are longer than 
  # the longest allowable vector in R (2^32 - 1), which is too small to record 
  # the every position in these genomes. There is some support for long vectors 
  # (2^64 - 1) but it is limited. In any case, this is an incredibly inefficient 
  # way of storing beta-values since 99% of sites have beta-value set to NA.
  # Computing per-sample-chromosome acfs seems doable (at least, for lag.max no 
  # bigger than a few thousand). But there's no way to post-hoc combine these.
  beta <- lapply(sampleNames(cometh), function(sn, cometh, sfsc){
    
    lapply(seq_len(nrow(sfsc), function(i, sfsc, cometh, s) {
      
      # WARNING: This is a _really_ inefficient way to make an Rle of the 
      # NA-padded beta-values.
      # Ideally, we'd just initiate the Rle without inflating/deflating the 
      # NA entries (which should be doable).
      beta <- rep(NA_real_, seqlengths(cometh)[sfsc$seqlevels[i]])
      idx <- as.vector(seqnames(cometh) == sfsc$seqlevels[i] & 
                         strand(cometh) == sfsc$strand[i])
      beta[start(cometh)[idx]] <- assay(cometh, 'beta')[idx, sn]

      ## Want to compute acf now. However, need all chromosomes worth of data.
      
    }, sfsc = sfsc, cometh = cometh, sn = sn))
  }, cometh = cometh, sfsc = sfsc)
  
  # Option 2: for-loop-based.
  # Loop over (sorted) positions of 1-tuples and record pairs that are within 
  # max_lag nt (or in the requested IPD-vector). Then compute correlations 
  # using cor.
  # How long should x and y be?
  # The current allocation is __incredibly__ conservative.
  # Dynamic allocation would be better but is expensive in R.
  # Perhaps move to Rcpp and use STL vectors, which I think have cheap vector 
  # extension.
  
  # Need to temporarily ignore strand of rowData _or_ do separately for each strand
  
  # Option A: Split by strand
  rdss <- split(rowData(cometh), strand(cometh))
  # x and y index the first and last element of each pair, respectively.
  endoapply(rdss, FUN = function(z, ipd) {
    
    s <- start(z)
    # x and y index the first and last element of each pair, respectively.
    n <- length(z)
    x <- y <- rep(NA_integer_, n * max(ipd))
    k <- 1
    for (i in seq_len(nrow(cometh) - 1)) {
      j <- i + 1
      while(j <= n){
        if((s[j] - s[i]) %in% ipd) {
          #print('blonko')
          x[k] <- i
          y[k] <- j
          j <- j + 1
          k <- k + 1
        } else if (s[j] - s[i] < max(ipd)){
          #print('still hope')
          j <- j + 1
        } else {
          #print('too big')
          break
        }
      }
    }
    
    # Drop NA-padding.
    ii <- !is.na(x) & !is.na(y)
    x <- x[ii]
    y <- y[ii]
    
    # Determine which pairs are "valid pairs", i.e. satisfy sfsc constraint.
    same_chrom <- seqnames(z)[x] == seqnames(z)[y]
    same_strand <- strand(z[x]) == strand(z[y])
    # same_feature_class could be generalised from 2 states to 4 states, i.e.,
    # from TRUE, FALSE 
    # to BOTH_IN, BOTH_OUT, FIRST_IN_SECOND_OUT, FIRST_OUT_SECOND_IN
    # Currently throwing away all pairs with either FIRST_IN_SECOND_OUT or 
    # FIRST_OUT_SECOND_IN.
    same_feature_class <- mcols(z)$ifc[x] == mcols(z)$ifc[y]
    # mls_vp is the index of valid pairs in mls (no NIL constraint, yet).
    z_vp <- as.vector(same_chrom & same_strand & same_feature_class)
    
    # Construct valid pairs from z
    z_vp_gr <- GRanges(seqnames = seqnames(z)[x][z_vp],
                       ranges = IRanges(start = start(z)[x][z_vp],
                                        end = end(z)[y][z_vp]),
                       strand = strand(z)[x][z_vp],
                       ifc = mcols(z)$ifc[x][z_vp]
    )
    
    # Find which valid pairs can be constructed from cometh.
    # Do this by matching cometh to start and end of z_vp_gr.
    # Then, remove those loci where only the start or end matches.
    # xx and yy index the first and second methlyation loci in each pair, 
    # respectively, constructed from cometh.
    s <- findOverlaps(query = cometh, subject = z_vp_gr, type = 'start')
    #   ss <- findOverlaps(query = cometh, 
    #                      subject = resize(z_vp_gr, width = 1, fix = 'start', 
    #                                       ignore.strand = TRUE),
    #                      type = 'any')
    #   stopifnot(identical(s, ss))
    e <- findOverlaps(query = cometh, subject = z_vp_gr, type = 'end')
    ## Might be unnecessary because queryHits(s) == x and queryHits(e) == y.
    xx <- queryHits(s)[!is.na(match(subjectHits(s), subjectHits(e)))]
    yy <- queryHits(e)[!is.na(match(subjectHits(e), subjectHits(s)))]
    
    stopifnot(identical(seqnames(cometh)[xx], seqnames(cometh)[yy]))
    stopifnot(identical(strand(cometh)[xx], strand(cometh)[yy]))
    
    # Construct valid pairs from cometh.
    cometh_vp_gr <- GRanges(seqnames = seqnames(cometh)[xx],
                            ranges = IRanges(start = start(cometh)[xx], 
                                             end = end(cometh[yy])),
                            strand = strand(cometh)[xx],
                            ifc = mcols(cometh)$ifc[xx])
    stopifnot(!any(duplicated(cometh_vp_gr)))
    ol <- findOverlaps(cometh_vp_gr, rdis_vp_gr, type = 'equal')
    
  }, ipd = seq_len(max_lag))

  
  
  
  # Option B: Temporarily ignoring strand
  rdis <- sort(rowData(cometh), ignore.strand = TRUE)
  s <- start(rdis)
  ipd <- seq_len(max_lag)
  # x and y index the first and last element of each pair, respectively.
  n <- length(rdis)
  x <- y <- rep(NA_integer_, n * max(ipd))
  k <- 1
  max_ipd <- max(ipd)
  for (i in seq_len(nrow(cometh) - 1)) {
    j <- i + 1
    #print(paste0('i = ', i))
    #print(paste0('j = ', j))
    while(j <= n){
      if((s[j] - s[i]) %in% ipd) {
        #print('blonko')
        x[k] <- i
        y[k] <- j
        j <- j + 1
        k <- k + 1
      } else if (s[j] - s[i] < max_ipd){
        #print('still hope')
        j <- j + 1
      } else {
        #print('too big')
        break
      }
    }
  }

  # Drop NA-padding.
  ii <- !is.na(x) & !is.na(y)
  x <- x[ii]
  y <- y[ii]
  
  # Determine which pairs are "valid pairs", i.e. satisfy sfsc constraint.
  same_chrom <- seqnames(rdis)[x] == seqnames(rdis)[y]
  same_strand <- strand(rdis[x]) == strand(rdis[y])
  # same_feature_class could be generalised from 2 states to 4 states, i.e.,
  # from TRUE, FALSE 
  # to BOTH_IN, BOTH_OUT, FIRST_IN_SECOND_OUT, FIRST_OUT_SECOND_IN
  # Currently throwing away all pairs with either FIRST_IN_SECOND_OUT or 
  # FIRST_OUT_SECOND_IN.
  same_feature_class <- mcols(rdis)$ifc[x] == mcols(rdis)$ifc[y]
  # mls_vp is the index of valid pairs in mls (no NIL constraint, yet).
  rdis_vp <- as.vector(same_chrom & same_strand & same_feature_class)
  # Remove those pairs that don't satisfy the sfsc contraint.
  x <- x[rdis_vp]
  y <- y[rdis_vp]
  
  # Construct valid pairs from rdis
  rdis_vp_gr <- GRanges(seqnames = seqnames(rdis)[x],
                        ranges = IRanges(start = start(rdis)[x],
                                         end = end(rdis)[y]),
                        strand = strand(rdis)[x],
                        ifc = mcols(rdis)$ifc[x]
                        )
  stopifnot(!any(duplicated(rdis_vp_gr)))

#### EVERYTHING FINE UP TO HERE; BUT CONSTRUCTING VALID PAIRS FROM COMETH IS WRONG ####
  
  # Find which valid pairs can be constructed from cometh.
  # Do this by matching cometh to start and end of rdis_vp_gr.
  # Then, remove those loci where only the start or end matches.
  # xx and yy index the first and second methlyation loci in each pair, 
  # respectively, constructed from cometh.
  stopifnot(!any(duplicated(cometh)))
  stopifnot(!any(duplicated(rdis_vp_gr)))
  s <- findOverlaps(query = cometh, subject = rdis_vp_gr, type = 'start', 
                    select = 'all')

  e <- findOverlaps(query = cometh, subject = rdis_vp_gr, type = 'end', 
                    select = 'all')

  ## Might be unnecessary because queryHits(s) == x and queryHits(e) == y.
  xx <- queryHits(s)[!is.na(match(subjectHits(s), subjectHits(e)))]
  stopifnot(identical(xx, queryHits(s)))
  yy <- queryHits(e)[!is.na(match(subjectHits(e), subjectHits(s)))]
  stopifnot(identical(yy, queryHits(e)))  

  stopifnot(identical(seqnames(cometh)[xx], seqnames(cometh)[yy]))
  stopifnot(identical(strand(cometh)[xx], strand(cometh)[yy]))
  
  # Construct valid pairs from cometh.
  ## TODO: This is producing duplicate GRanges; why? It shouldn't.
  cometh_vp_gr <- GRanges(seqnames = seqnames(cometh)[xx],
                          ranges = IRanges(start = start(cometh)[xx], 
                                           end = end(cometh)[yy]),
                          strand = strand(cometh)[xx],
                          ifc = mcols(cometh)$ifc[xx])

  stopifnot(!any(duplicated(cometh_vp_gr)))
  # d_cometh_vp_gr indexes the duplicate elements in cometh_vp_gr
  d_cometh_vp_gr <- which(duplicated(cometh_vp_gr) | duplicated(cometh_vp_gr, fromLast = TRUE))
  cometh_vp_gr[d_cometh_vp_gr]
  # qh is the correponding element in rdis_vp_gr
  qh <- queryHits(findOverlaps(rdis_vp_gr, cometh_vp_gr[duplicated(cometh_vp_gr)], type = 'equal'))
  xx[d_cometh_vp_gr]
  yy[d_cometh_vp_gr]  
  x[sort(c(qh, qh + 1))]
  y[sort(c(qh, qh + 1))]
  cbind(xx[d_cometh_vp_gr], yy[d_cometh_vp_gr], x[sort(c(qh, qh + 1))], y[sort(c(qh, qh + 1))])





  ol <- findOverlaps(cometh_vp_gr, rdis_vp_gr, type = 'equal')
  
  # This is a sanity check and shouldn't strictly be required.
  # Check that cometh_vp_gr and rdis_vp_gr have identical pairs.
  sanity_check <- countQueryHits(ol)
  stopifnot(all(sanity_check == 1L)) 
  sanity_check2 <- countSubjectHits(ol)
  stopifnot(all(sanity_check2 == 1L))
}