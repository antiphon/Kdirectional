
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;


// nearest neighbour vectors

// [[Rcpp::export]]
List c_angles(NumericMatrix x, IntegerVector from, IntegerVector to) {
    int i, j, l, f, t;
    int n = x.nrow();
    int dim = x.ncol();
    double dx, dy, dz, d, ang, nnd;
    NumericVector angles(n);
    NumericVector angles2(n);
    NumericVector nndists(n);
    int nn;
    
    for(i=0; i < from.length(); i++) {
      f = from[i]-1;
      nnd = DBL_MAX;
      for(j=0; j < to.length(); j++) {
        t = to[j]-1;
        if(t!=f) {
          d=0;
          for(l=0; l < dim; l++)  d += pow(x(f,l)-x(t,l), 2);
          d = sqrt(d);
          if (d < nnd){
            nn = t;
            nnd = d;
          }
        }
      }
      t = nn;
      dx = x(t,0) - x(f,0);
      dy = x(t,1) - x(f,1);
      ang = atan2(dy, dx);
      if(ang<0) ang = 2*PI+ang;
      angles[f] = ang;
      nndists[f] = nnd;
      if(dim==3){
        dz = x(t,2) - x(f,2);
        ang = acos(dz/nnd);
        angles2[f] = ang;
      }
    }
    if(dim==3) return List::create(angles, angles2, nndists);
    return List::create(angles, nndists);
}
