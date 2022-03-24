
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List c_angles_in_a_cone(NumericMatrix x, NumericVector unit, double theta, 
                        IntegerVector from, IntegerVector to, bool antipodal = false) {
    int i, j, l, f, t;
    int n = x.nrow();
    int dim = x.ncol();
    double dx, dy, dz, d, ang, dang, nnd, dot;
    NumericVector angles(n);
    NumericVector angles2(n);
    NumericVector nndists(n);
    int nn;
    int found = 0;
    for(i=0; i < from.length(); i++) {
      f = from(i)-1;
      nnd = DBL_MAX;
      found = 0; 
      for(j=0; j < to.length(); j++) {
        t = to(j)-1;
        if(t!=f) {
          d=0;
          for(l=0; l < dim; l++)  d += pow(x(f,l)-x(t,l), 2);
          d = sqrt(d);
          // in the cone?
          dot = 0;
          for(l=0; l < dim; l++)  dot += (x(t,l)-x(f,l)) * unit(l);
          dang = acos(dot/d);
          if(antipodal) dang = min(dang, M_PI - dang);
          if (dang < theta & d < nnd){
            nn = t;
            nnd = d;
            found = 1;
          }
        }
      }
      if(found){
        t = nn;
        dx = x(t,0) - x(f,0);
        dy = x(t,1) - x(f,1);
        ang = atan2(dy, dx);
        if(ang<0) ang = 2*M_PI+ang;
        angles[f] = ang;
        nndists[f] = nnd;
        if(dim==3){
          dz = x(t,2) - x(f,2);
          ang = acos(dz/nnd);
          angles2[f] = ang;
        }
      }
      else{
        nndists[f] = nnd;
      }
    }
    if(dim==3) return List::create(angles, angles2, nndists);
    return List::create(angles, nndists);
}
