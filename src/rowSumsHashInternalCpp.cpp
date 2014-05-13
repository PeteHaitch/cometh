#include <Rcpp.h>
using namespace Rcpp;

//' An internal helper function used in finding duplicate rows of an integer matrix.
//' 
//' @details
//' Not actually used in cometh package but written during development of mTuplesHashCpp.
//' 
//' @param x An integer matrix.
//' 
//' @return A logical vector. TRUE if the corresponding row of \code{x} is a 
//' potential duplicate and FALSE otherwise. 
//' \emph{NOTE:} TRUE does not mean that it is a duplicate, merely that it could be.
//' Further checking is required of elements with TRUE, e.g. using the 
//' \code{\link[base]{duplicated.array}} method.
//' 
//' @keywords internal
//' 
//' 
// [[Rcpp::export]]
LogicalVector rowSumsHashInternalCpp(IntegerMatrix x){
  
  int nr = x.nrow();
  int nc = x.ncol();
  LogicalVector val(nr);
  IntegerVector hash(nr);
  
  // Compute a hash, which is just the rowSums
  for (int i = 0; i < nr; i++){
    int row_total = 0;
    for (int j = 0; j < nc; j++){
      row_total += x(i, j);
    }
    hash[i] = row_total;
  }
  
  // Check for duplicates in hash
  // rev(duplicated(rev(hash))) is equivalent to 
  // duplicated(hash, fromLast = TRUE)
  val = duplicated(hash) | rev(duplicated(rev(hash)));
  
  // NOTE: Still need to check duplicate hashes to see if they are truly 
  // duplicate elements
  
  return val;
}