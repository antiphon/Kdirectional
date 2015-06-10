
#include <Rcpp.h>
#include <vector>
#include "helpers.h"
using namespace Rcpp;

//#define PI (3.141592653589793)

//double min(double, double);


// Upper triangle vector

// [[Rcpp::export]]
NumericVector c_translation_weights(NumericMatrix x, NumericMatrix bbox) {
    
    int i, j, k, ind;
    double wij = 0;
    int n = x.nrow();
    int dim = x.ncol();
    NumericVector boxlen(dim);
    NumericVector w( n*(n-1) * 0.5 );
    
    ind = 0;
    
    for(i=0; i < dim; i++) boxlen[i] = bbox(1,i) - bbox(0,i);
    
    for(i=0; i < n-1; i++) {
      for(j=i+1; j < n; j++) {
        wij = 1;
        for(k=0; k < dim; k++) wij *= boxlen(k) - abs(x(i,k)-x(j,k));
        w(ind) = wij;
        ind ++ ; 
        }
      }
    return w;
}

//double min(double a, double b){
//  if(a < b ) return a;
//  return b;
//}
