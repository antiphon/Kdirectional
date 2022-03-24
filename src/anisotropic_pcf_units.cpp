// Unit vector version
#include <Rcpp.h>
#include <vector>
#include "helpers.h"

using namespace Rcpp;

//#define PI (3.141592653589793)


// [[Rcpp::export]]
List c_anisotropic_unit_pcf(NumericMatrix x,
                                     NumericMatrix r, 
                                     NumericVector h, 
                                     NumericMatrix bbox, 
                                     int correction) {
    int i, j, k, l;
    double d, ang, v, w, dot, diff, W;
    int nr, nvals;
    
    int n = x.nrow();
    int dim = x.ncol();
    
    
    // number of grid vectors
    nr = r.nrow();
    nvals  = nr;
    // pre-compute the lengths
    NumericVector rlength(nr);
    for(k=0; k < nr; k++) {
      d=0;
      for(i=0; i < dim; i++) d += r(k,i)*r(k,i);
      rlength(k) = sqrt(d);
    }
    IntegerVector counts(nr);
    NumericVector val(nvals);
    
    NumericVector boxlen(dim);
    for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
    
    for(i=0; i < n-1; i++) {
      for(j=i+1; j < n; j++) {
        // distance
        d=0;
        for(l=0; l < dim; l++)  d += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        d = sqrt(d);
        // border correction
        w = 1;
        if(correction>1){
          for(k=0; k < dim; k++) w *= boxlen(k) - fabs(x(i,k)-x(j,k));
        }
        // contribution to each direction
        for(k=0; k < nr; k++) {
          dot = 0;
          for(l=0; l < dim; l++)  dot += (x(j,l)-x(i,l))*r(k,l);
          dot = dot/(d * rlength(k));
          ang = acos(dot);
          ang = fmin(ang, M_PI-ang); // antipodal
          diff = fabs( d - rlength(k) );
          if(ang < h(1) & diff < h(0)){
            //W = pow(rlength(k), dim-1);
            W = pow(d, dim-1);
            v = kernel_epa(diff, h(0)) * kernel_epa(ang, h(1)) / W;
            val[k] += v/w;
            counts(k)++;
          }
        }
      }
    }
    List res(2);
    res(0) = wrap(val);
    res(1) = wrap(counts);
    return res;
}

