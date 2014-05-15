#include <Rcpp.h>
using namespace Rcpp;

//' An internal helper function used in finding duplicate m-tuples.
//'
//' The m-tuple should have already been converted to an integer representation,
//' namely an integer vector for the chromosome, an integer vector 
//' for the strand and an integer matrix for the  positions.
//' 
//' @param a An integer vector of length n. An integer representation of the 
//' chromosome of each m-tuple.
//' @param b An integer vector of length n. An integer representation of the 
//' strand of each m-tuple.
//' @param C An integer matrix with n rows. Each row represents the positions 
//' of an m-tuple.
//' 
//' @return A logical vector. TRUE if the corresponding row of \code{x} is a 
//' potential duplicate and FALSE otherwise. 
//' \emph{NOTE:} TRUE does not mean that it is a duplicate, merely that it is a 
//' candidate. Further checking is required of these candidates, e.g. using the 
//' \code{\link[base]{duplicated.array}} method.
//' 
// [[Rcpp::export(".candidateDuplicateMTuplesCpp")]]
LogicalVector candidateDuplicateMTuplesCpp(IntegerVector a, IntegerVector b, IntegerMatrix C){
  
  int nr = C.nrow();
  int nc = C.ncol();
  LogicalVector val(nr);
  IntegerVector hash(nr);
  
  // Compute a hash, which is just the "rowSums"
  for (int i = 0; i < nr; i++){
    int row_total = a[i] + b[i];
    for (int j = 0; j < nc; j++){
      row_total += C(i, j);
    }
    hash[i] = row_total;
  }
  
  // Check for duplicates in hash
  // rev(duplicated(rev(hash))) is equivalent to duplicated(hash, fromLast = TRUE)
  val = duplicated(hash) | rev(duplicated(rev(hash)));
  
  // NOTE: Still need to check candidate duplicates to see if they are truly 
  // duplicate m-tuples. This is currently done outside of this function.
  
  return val;
}
