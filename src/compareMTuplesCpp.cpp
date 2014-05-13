#include <Rcpp.h>
using namespace Rcpp;

//' An internal function used to compare m-tuples.
//' 
//' @details 
//' The m-tuple should have already been converted to an integer representation,
//' namely an integer vector for the difference in chromosome, an integer vector 
//' for the difference in strand and an integer matrix for the difference in 
//' positions.
//' 
//' @param a An integer vector of length n. An integer representation of the difference in chromosomes of each m-tuple.
//' @param b An integer vector of length n. An integer representation of the difference in strand of each m-tuple.
//' @param C An integer matrix with n rows. Each row represents the difference in positions of each m-tuple. 
//' 
//' 
//' @return An integer vector where each element is the comparison of a pair
//' of m-tuples.
//' If the first m-tuple in the pair is "<" than the second m-tuple then the
//' return value for that element is < 0, if the first m-tuple in the pair is 
//' "==" the second m-tuple then the return value is 0, and if the first
//' m-tuple is ">" that the second m-tuple then the return value is > 0.
//' 
//' @keywords internal
//'
//' 
// [[Rcpp::export]]
IntegerVector compareMTuplesCpp(IntegerVector a, IntegerVector b, IntegerMatrix C){
  int n = a.size();
  IntegerVector val(n);
  int nc = C.ncol();
  
  for (int i = 0; i < n; i++){
    if (a[i] != 0){
      val[i] = a[i];
      continue;
    }
    if (b[i] != 0){
      val[i] = b[i];
      continue;
    }
    // Have to do it in this weird order otherwise can't "continue" the outer loop
    val[i] = C(i, nc - 1); 
    for (int j = 0; j < (nc - 1); j++){
      if (C(i, j) != 0){
        val[i] = C(i, j);
        break;
      }
    }
  }
  return val;
}
