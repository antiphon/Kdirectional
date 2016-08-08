#include <Rcpp.h>
#include "helpers.h"
#include "epa_integral.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector intensity_at_other_points_c(NumericMatrix x, 
                                NumericMatrix other,
                                NumericMatrix bbox, 
                                double bw_r,
                                int n=0){
  
  NumericVector out(other.nrow());
  
  NumericVector epaw = epa_integral(other, bbox, bw_r, n);
  
  int i,j,k;
  
  double dim = x.ncol();
  double d, w;
  
  for(i=0; i < x.nrow(); i++) {
    for(j=0; j < other.nrow(); j++) {
      d = 0;
      for(k=0; k < dim; k++) d += pow( x(i,k)-other(j,k), 2);
      d = sqrt(d);
      w = kernel_epa(d, bw_r);
      if(w>0) out(j) += w/epaw(j);
    }
  }
  return out;
}

