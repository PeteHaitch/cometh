#include <Rcpp.h>
using namespace Rcpp;

// allRowsSorted is adapted from http://stackoverflow.com/a/7601857
// [[Rcpp::export]]
bool allRowsSorted(NumericMatrix A){
  int nrow = A.nrow(), ncol = A.ncol();
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol - 1; j++) {
      if (A(i, j) > A(i, j + 1)) { 
        return false; 
      }
    }
  }
  return true;
}