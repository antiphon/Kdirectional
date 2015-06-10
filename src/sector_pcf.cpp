
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

//#define PI (3.141592653589793)

double kernel(double, double);
double abss(double);

// [[Rcpp::export]]
NumericVector c_sector_pcf(NumericMatrix x, NumericVector u, double theta,  
                           NumericVector r, double h, NumericMatrix bbox, int correction) {
    int i, j, k, l;
    double dot, d, ang, v, w;
    
    int n = x.nrow();
    int dim = x.ncol();
    
    int nr = r.length();
    
    NumericVector val(nr);
    
    NumericVector boxlen(dim);
    for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
    
    int count = 0;
    
    for(i=0; i < n-1; i++) {
      for(j=i+1; j < n; j++) {
        // distance
        d=0;
        for(l=0; l < dim; l++)  d += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        d = sqrt(d);
        //check angle
        dot = 0;
        for(l=0; l < dim; l++)  dot += (x(i,l)-x(j,l))*u[l];
        ang = acos(dot/d);
        ang = min(ang, PI-ang);
        if(ang < theta){ // this pair is aligned
          // translation correction?
          w = 1;
          if(correction>1){
            for(k=0; k < dim; k++) w *= boxlen(k) - abss(x(i,k)-x(j,k));
          }
          for(k=0; k < nr; k++){
            v = kernel(abss(d-r[k]), h);
            val[k] += v/w;
          }
          count++;
        }
      }
    }
    return val;
}

double kernel(double r, double h){
  if(r > h) return 0;
  //return 0.5/h;
  // Epanechnikov
  return 0.75 * (1-r*r/(h*h))/h;
}

double abss(double x){
  if(x < 0) return -x;
  return x;
}

//double min(double a, double b){
//  if(a < b ) return a;
//  return b;
//}
