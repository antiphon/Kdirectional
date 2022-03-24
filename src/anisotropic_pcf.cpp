
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

//#define PI (3.141592653589793)

double kernell(double, double);
double abss_a(double);

// [[Rcpp::export]]
NumericVector c_anisotropic_pcf(NumericMatrix x, List theta,  
                                NumericVector r, NumericVector h, NumericMatrix bbox, int correction) {
    int i, j, k, l, m, ind;
    double d, ang, ang2, v, w;
    int nang, nang2;
    
    int n = x.nrow();
    int dim = x.ncol();
    
    // number of ranges
    int nr = r.length();
    int nvals = nr;
    
    // number of angles
    NumericVector angs = theta(0);
    NumericVector angs2;
    
    nang = angs.length();
    nvals *= nang;
    
    if(dim==3) {
      angs2 = theta(1);
      nang2 = angs2.length();
      nvals *= nang2;
    }
    NumericVector val(nvals);
    
    NumericVector boxlen(dim);
    for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
    
    for(i=0; i < n-1; i++) {
      for(j=i+1; j < n; j++) {
        // distance
        d=0;
        for(l=0; l < dim; l++)  d += (x(i,l)-x(j,l))*(x(i,l)-x(j,l));
        d = sqrt(d);
        // translation correction
        w = 1;
        if(correction>1){
          for(k=0; k < dim; k++) w *= boxlen(k) - abss_a(x(i,k)-x(j,k));
        }
        // angle
        ang = atan2(x(i,1)-x(j,1), x(i,0)-x(j,0));
        if(dim==2){
          if(ang < 0) ang = M_PI + ang;
          // contribute
          for(k=0; k < nr; k++){
            for(l=0; l < nang; l++){
             v = kernell(abss_a(d-r[k]), h(0)) * kernell(abss_a(ang - angs(l)), h(1));
             v += kernell(abss_a(d-r[k]), h(0)) * kernell(abss_a(M_PI - ang - angs(l)), h(1));
             ind = k + l*nr;
             val[ind] += 0.5*v/w; 
            }
          }
        }
        else{
          ang2 = acos( (x(i,2) - x(j,2)) /d);
          if(ang < 0){ // flip to positive side other side
            ang = M_PI + ang;
            ang2 = M_PI - ang2;
          }
          // contribute
          for(k=0; k < nr; k++){
            for(l=0; l < nang; l++){
              for(m=0; m < nang2; m++){
                v = kernell(abss_a(d-r(k)), h(0)) * kernell(abss_a(ang-angs(l)), h(1)) *  kernell(abss_a(ang2-angs2(m)), h(1));
                v += kernell(abss_a(d-r(k)), h(0)) * kernell(abss_a(M_PI + ang-angs(l)), h(1)) *  kernell(abss_a(M_PI-ang2-angs2(m)), h(1));
                ind = k + l*nr + m * nr * nang;
                val[ind] += 0.5*v/w;
              }
            }
          }
        }
      }
    }
    return val;
}

double kernell(double r, double h){
  if(r > h) return 0;
  //return 0.5/h;
  // Epanechnikov
  return 0.75 * (1-r*r/(h*h))/h;
}

double abss_a(double x){
  if(x < 0) return -x;
  return x;
}

//double min(double a, double b){
//  if(a < b ) return a;
//  return b;
//}
