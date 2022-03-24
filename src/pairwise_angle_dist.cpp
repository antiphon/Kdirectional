
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix c_pairwise_dist_angle(NumericMatrix x) {
  int i, j, ind, l;
  double dx, dy, dz;
  int n = x.nrow();
  int dim = x.ncol();
  double d, ang;
  // Upper triangle of angles and distances:
  NumericMatrix val( n*(n-1) * 0.5, dim);
  ind = 0;
  for(i=0; i < n - 1; i++) {
    for(j=i+1; j < n; j++) {
      d=0;
      for(l=0; l < dim; l++)  d += pow(x(j,l)-x(i,l), 2);
      d = sqrt(d);
      dx = x(j,0) - x(i,0);
      dy = x(j,1) - x(i,1);
      ang = atan2(dy, dx);
      if(ang<0) ang = 2*M_PI+ang;
      // upper triangle location      
      val(ind,0) = d;
      val(ind,1) = ang;
      if(dim==3){
        dz = x(j,2) - x(i,2);
        ang = acos(dz/d);
        val(ind,2) = ang;
      }
      ind++;
    }
  }
  return val;
}
