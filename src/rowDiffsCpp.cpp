#include <Rcpp.h>
using namespace Rcpp;

//' An internal function to compute differences for each row in a matrix.
//' 
//' @param A An integer \eqn{N \times K} matrix.
//' 
//' @keywords internal
//'
//' @return An integer \eqn{N \times (K - 1)} matrix.
//' 
//' @note This function silently coerces numeric matrices to integer matrices
//' and does integer subtraction. \emph{This will give unexpected results but
//' it's not a problem for me since I only use it on integer matrices.}
// [[Rcpp::export(".rowDiffsCpp")]]
IntegerMatrix rowDiffsCpp(IntegerMatrix A){
  int nrow = A.nrow(), ncol = A.ncol();
  IntegerMatrix D(nrow, ncol - 1);
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol - 1; j++){
      D(i, j) = A(i, j + 1) - A(i, j);
    }
  }
  return D;
}