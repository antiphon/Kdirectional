#include <Rcpp.h>
#include "helpers.h"
#include "epa_integral.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector intensity_at_points_c(NumericMatrix x, 
                                  NumericMatrix bbox, 
                                  double bw_r,
                                  int b=1){
  
  NumericVector out(x.nrow());
  
  int i,j,k;
  
  double dim = x.ncol();
  double d, w;
  
  NumericVector epaw = epa_integral(x, bbox, bw_r, b);
  
  //return epaw;
  
  for(i=0; i < x.nrow()-1; i++) {
    for(j=i+1; j < x.nrow(); j++) {
      d = 0;
      for(k=0; k < dim; k++) d += pow( x.at(i,k)-x.at(j,k), 2);
      d = sqrt(d);
      w = kernel_epa(d, bw_r);
      if(w >0) {
        out.at(i) += w/epaw.at(i);
        out.at(j) += w/epaw.at(j);
      }
    }
  }
  
  return out;
  
}

