
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

//#define PI (3.141592653589793)

//double min(double, double);

// [[Rcpp::export]]
List c_directed_geom(NumericMatrix x, NumericVector u, double theta,  double r, IntegerVector from, IntegerVector to) {
    int i, j, l, f, t;
    int n = x.nrow();
    int dim = x.ncol();
    double dot, d, d2, ang, r2 = r*r;
    List nlist(n);
    std::vector<int> neighs;
    for(i=0; i < from.length(); i++) {
      neighs.clear();
      f = from[i]-1;
      for(j=0; j < to.length(); j++) {
        t = to[j]-1;
        if(t!=f) {
            // first check distance
          d2=0;
          for(l=0; l < dim; l++)  d2 += pow(x(f,l)-x(t,l), 2);
          if(d2 < r2){
            d = sqrt(d2);
            //then check angle
            dot = 0;
            for(l=0; l < dim; l++)  dot += (x(t,l)-x(f,l))*u[l];
            ang = acos(dot/d);
            ang = min(ang, M_PI-ang);
            if(ang < theta)
              neighs.push_back(t+1); 
          }
        }
      }
      nlist(f) = wrap( neighs );
    }
    return nlist;
}

//double min(double a, double b){
//  if(a < b ) return a;
//  return b;
//}
