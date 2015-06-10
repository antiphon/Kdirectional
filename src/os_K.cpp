// Ohser-Stoyan directed K-function

#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

//#define PI (3.141592653589793)

//double kernel(double, double);
double absss(double);

// [[Rcpp::export]]
NumericVector c_oh_K(NumericMatrix x, List theta, NumericVector r, NumericMatrix bbox) {
    int i, j, k, ri, ai, a2i, nang2, ii, jj, kk;
    double d, ang, ang2, v, w;
    
    int n = x.nrow();
    int dim = x.ncol();
    
    // number of ranges
    int nr = r.length();
    int nvals = nr;
    
    // number of angles
    NumericVector angs = theta(0);
    NumericVector angs2;
    
    int nang = angs.length();
    nvals *= nang;
    
    if(dim==3) {
      angs2 = theta(1);
      nang2 = angs2.length();
      nvals *= nang2;
    }
    printf("nr %i, nang %i, nval %i\n", nr, nang, nvals);
    NumericVector vals(nvals);
    for(i=0; i<nvals; i++) vals[i]=0.0;
    NumericVector boxlen(dim);
    for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
    
    int cr = 0, ca = 0, ca2 = 0;
    
    
    // main loop
    for(i=0; i < n-1; i++) {
      for(j=i+1; j < n; j++) {
        // distance
        d=0;
        for(k=0; k < dim; k++)  d += (x(i,k)-x(j,k))*(x(i,k)-x(j,k));
        d = sqrt(d);
        // angle
        ang = atan2(x(i,1)-x(j,1), x(i,0)-x(j,0));
        if(dim==3){
          ang2 = acos( (x(i,2) - x(j,2)) /d);
        }
        if(ang < 0) ang = PI+ang;
        //
        // Then attribute to correct grid slot
        // range:
        ri = nr;
        for(k=0; k < nr; k++) if(r(k) < d) ri = k;
        // pair is within range
        if(ri < nr - 1){
          cr++;
          // check angle 1:
          ai = nang;
          for(k=0; k < nang; k++) if(angs(k) < ang) ai = k;
          // within first angles
          if(ai < nang - 1 ) {
            ca++;
            if(dim==2){ // 2D
              //printf("%f/%f, %f/%f \n", d, r(ri), ang, angs(ai));
              // Translation correction
              w = 1;
              for(k=0; k < dim; k++) w *= boxlen(k) - absss(x(i,k)-x(j,k));          
              v = 1/w;
              for(ii=ri; ii < nr; ii++)
                for(jj=ai; jj < nang; jj++){
                  k = ii + jj * nr;
                  vals[k] += v;
                }
            }
            else{ // 3D.
              // check angle 2:
              a2i = nang2;
              for(k=0; k < nang2; k++) if(angs2(k) < ang2) a2i = k;
              if(a2i < nang2){
                ca2++;
                // Translation correction
                w = 1;
                for(k=0; k < dim; k++) w *= boxlen(k) - absss(x(i,k)-x(j,k));          
                v = 1/w;
                for(ii=ri; ii < nr; ii++)
                  for(jj=ai; jj < nang; jj++)
                    for(kk=a2i; kk < nang2; kk++){
                      k = ii + jj * nr + kk * nr * nang;
                      vals[k] += v;
                    }
                
              }
            }
          }
        }
      }
    }
    //printf("r ok %i, a ok %i, a2 ok %i\n", cr, ca, ca2);
    return vals;
}


double absss(double x){
  if(x < 0) return -x;
  return x;
}

//double min(double a, double b){
//  if(a < b ) return a;
//  return b;
//}
