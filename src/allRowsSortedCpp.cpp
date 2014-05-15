#include <Rcpp.h>
using namespace Rcpp;

//' An internal helper to check that each row of a numeric matrix is sorted in 
//' strictly increasing order
//' 
//' @param A A numeric matrix
//' 
//' @return TRUE if each row is sorted in strictly increasing order, FALSE otherwise
//' 
//' @details 
//' .allRowsSortedCpp is adapted from http://stackoverflow.com/a/7601857.
//' Strict inequalities are required in order to avoid m-tuples with repeated 
//' positions (https://github.com/PeteHaitch/cometh/issues/8).
//' 
//' @keywords internal
//' 
// [[Rcpp::export(".allRowsSortedCpp")]]
bool allRowsSortedCpp(NumericMatrix A){
  int nrow = A.nrow(), ncol = A.ncol();
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol - 1; j++){
      if (A(i, j) >= A(i, j + 1)){
        return false; 
      }
    }
  }
  return true;
}