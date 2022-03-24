
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix c_pairwise_dist_angle_subset(NumericMatrix x, IntegerVector from, IntegerVector to) {
  int i, j, ind, l, t, f;
  double dx, dy, dz;
  int dim = x.ncol();
  double d, ang;
  int nfrom = from.length();
  int nto = to.length();
  // Upper triangle no longer meaningful structure, so just a storage.
  NumericMatrix val( nfrom * nto, dim);
  ind = 0;
  for(i=0; i < nfrom; i++) {
    f = from[i]-1;
    for(j=0; j < nto; j++) {
      t = to[j]-1;
      if(f!=t){
        d=0;
        for(l=0; l < dim; l++)  d += pow(x(f,l)-x(t,l), 2);
        d = sqrt(d);
        dx = x(f,0) - x(t,0);
        dy = x(f,1) - x(t,1);
        ang = atan2(dy, dx);
        if(ang<0) ang = 2*M_PI+ang;
        // upper triangle location      
        val(ind,0) = d;
        val(ind,1) = ang;
        if(dim==3){
          dz = x(f,2) - x(t,2);
          ang = acos(dz/d);
          val(ind,2) = ang;
        }
        ind++;
      }
    }
  }
  return val;
}
