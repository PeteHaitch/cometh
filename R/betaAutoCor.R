## The case of NIL = Inf is best solved by "traditional" autocorrelation 
## methods, i.e. treat the beta-values as a vector (perhaps suitably padded) 
## by NAs for all positions, not just methylation loci, without an observed 
## beta-value. I would still need a way to combine data from different 
## chromosomes. An existing algorithm will surely do a better job of computing
## fast and accurante autocorrelations in this case than anything I write 
## myself.
## TODO: Write the internal function to compute autocorrelation of beta-values.
## Will use acf (or some other fast autocorrelation function).
## TODO: Read Venables and Ripley for details of acf
