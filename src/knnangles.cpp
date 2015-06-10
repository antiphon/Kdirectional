
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List c_knnangles(NumericMatrix x, int k, IntegerVector from, IntegerVector to) {
    int i, j, l, m, pos, f, t;
    int n = x.nrow();
    int dim = x.ncol();
    double dx, dy, dz, d, ang, nnd;
    NumericVector angles(n);
    NumericVector angles2(n);
    NumericVector nndists(n);
    NumericVector knnd(k);
    IntegerVector knni(k);
    
    for(i=0; i < from.length(); i++) {
      f = from[i]-1;
      for(l=0; l< k; l++){
        knnd[l] = DBL_MAX;
        knni[l] = -1;
      }
      for(j=0; j < to.length(); j++) {
        t = to[j]-1;
        if(t!=f) {
          d=0;
          for(l=0; l < dim; l++)  d += pow(x(f,l)-x(t,l), 2);
          pos = k;
          for(l=0; l < k; l++){
            if (d < knnd[l]){
              pos = l;
              break;
            }
          }
          if(pos < k){
            for(m=k-1; m > pos; m--){
                knnd[m] = knnd[m-1];
                knni[m] = knni[m-1];
              }
              knnd[pos] = d;
              knni[pos] = t;
          }
        }
      }
      t = knni[k-1];
      nnd = sqrt(knnd[k-1]);
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
