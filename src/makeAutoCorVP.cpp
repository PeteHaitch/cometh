#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".makeAutoCorVP")]]

List makeAutoCorVP(IntegerVector s, IntegerVector ipd) {
  
  // Initiate vectors to store results.
  std::vector<int> x;
  std::vector<int> y;

  // Create variables used in the for-loop.
  int n = s.size();
  int j = 0;
  int k = 1;
  int max_ipd = max(ipd);
  
  // Loop over start positions, s, and find pairs the the specified ipd.
  for (int i = 0; i < (n - 1); i++) {
    j = i + 1;
    while (j <= (n - 1)) {
      /* 
      Ideally, would use (s[j] - s[i]) %in% ipd but not %in% sugar operator.
      Instead, use this monstrosity. Note need to wrap in is_true(); 
      see http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-May/005807.html
      */
      if (is_true(any(ipd == (s[j] - s[i])))) {
        
        /* 
        (i, j) are a pair, so add to x and y.
        Record (i + 1, j + 1) because R vectors are 1-based and C++ vectors 
        are 0-based. 
        */
        x.push_back(i + 1); 
        y.push_back(j + 1);
        j += 1;
        k += 1;
      } else if ((s[j] - s[i]) < max_ipd) {
        // (i, j) aren't a pair but (i, j + 1) might be, so keep looking.
        j += 1;
      } else {
        /*
        (i, j) aren't a pair and (i, j + 1) can't be, so move onto pairs 
        with x = i + 1.
        */
        break;
      }
    }
  }
  
  return List::create(_["x"] = x, _["y"] = y);
  
}