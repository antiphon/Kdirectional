#include <Rcpp.h>
#include "helpers.h"
// #include "epa_integral.h"
using namespace Rcpp;


/* Utility functions for intensity estimation
*/

#define SQRT_TWOPI 2.506628

// all kernels 1D. Use product all D.

double epakern(double r){
  if(r < 0) return epakern(-r);
  if(r > 1) return 0;
  return 0.75 * (1 - r*r);
}

double epaCumu(double u){
  if(u < -1) return 0;
  if(u > 1) return 1;
  return  0.75 * (u - u*u*u/3.0 + 2.0/3.0);
}

double epaRange(double h) {
  return h;
}

double gausskern(double r){
  return exp(-0.5*(r*r))/SQRT_TWOPI;
}

double gaussCumu(double u){
  return 0.5 + 0.5 * tanh(u * 0.85);// approx by logistic
}

double gaussRange(double h) {
  return 3.0 * h;
}

double bfun_local(double a, double b){
  return a;
}

double bfun_global(double a, double b){
  return b;
}

double bfun_none(double a, double b){
  return 1.0;
}

/* end utility functions */

// predefine
NumericVector int_at_anywhere_toroidal_c(NumericMatrix x,
                                         NumericMatrix loc,
                                         NumericMatrix bbox,
                                         double bw,
                                         int kernel);

NumericVector int_at_points_toroidal_c(NumericMatrix x,
                                       NumericMatrix bbox,
                                       double bw,
                                       int kernel,
                                       int loo);
  

/*
intensity at data points for several bandwidths. no edge correction, no leave one out.
*/

// [[Rcpp::export]]

NumericMatrix int_bw_c(NumericMatrix x,
                       NumericVector h,
                       int kernel = 0) { // 0 = Gaussian 1 = epanechnikov
  int i, j, k;
  int nh = h.size();
  int dim = x.ncol();
  int n = x.nrow();
  double d, w;
  
  NumericMatrix out(n, nh);
  
  double (*kern)(double);
  double hmax = max(h) * 9;
  if(kernel == 0){
    kern = &gausskern;
    hmax = gaussRange(max(h));
  }
  else {
    kern = &epakern;
    hmax = epaRange(max(h));
  }
  
  
  // compute
  for(i = 0; i < n; i++){
    for(j = i; j < n; j++) {
      d = 0;
      for(k = 0; k < dim; k++) d += pow(x(i,k)-x(j,k),2);
      d = sqrt(d);
      if(d < hmax){
        for(k = 0; k < nh; k++) {
          w = kern(d/h(k)) / pow(h(k), dim);
          out(i,k) += w;
          if(i < j) out(j,k) += w;
        }
      }
    }
  }
  return out;
}


/* Kernel intensity estimation at data points */

// [[Rcpp::export]]
NumericVector int_at_points_c(NumericMatrix x,
                              NumericMatrix bbox,
                              double bw,
                              int kernel = 0, //0 gauss, 1 prod epa
                              int border = 1, // 0 none, 1 local, 2 global 3 toroidal
                              int loo = 0){ // 1 leave one out, 0 not
  
  NumericVector out(x.nrow());
  
  int i,j,k;
  
  double dim = x.ncol();
  double w, wbi, wbj;
  double a, b;
  
  double (*bcf)(double ai, double bi) = bfun_none;
  // border = 0 no edge correction
  if(border == 1){
    bcf = bfun_local;
  }
  else if(border == 2){
    bcf = bfun_global;
  }
  else if(border == 3){ // toroidal
    return int_at_points_toroidal_c(x, bbox, bw, kernel, loo);
  }
  
  double (*kernCumu)(double d) = &gaussCumu;
  double (*kern)(double d) = &gausskern;
  double (*kernRange)(double d) = &gaussRange;
  if(kernel == 1) {
    kernCumu = &epaCumu;
    kern = &epakern;
    kernRange = &epaRange;
  }
  
  double bwd = pow(bw, dim);
  
  for(i=0; i < x.nrow(); i++) {
    wbi = 1; //
    for(k = 0; k < dim; k++) {
      a = bbox(0,k)-x(i,k);
      b = bbox(1,k)-x(i,k);
      wbi *= kernCumu(b/bw) - kernCumu(a/bw);
    }
    for(j=i+loo; j < x.nrow(); j++) {
      {
        wbj = 1;
        w = 1/bwd;
        for(k=0; k < dim; k++){
          a = x(i,k)-x(j,k);
          w *= kern(a/bw);
          a = bbox(0,k)-x(j,k);
          b = bbox(1,k)-x(j,k);
          wbj *= kernCumu(b/bw) - kernCumu(a/bw);
        }
        out(i) += w / bcf(wbi, wbj);
        if(i<j) out(j) += w / bcf(wbj, wbi);
      }
    }
  }
  
  return out;
  
}



NumericVector int_at_points_toroidal_c(NumericMatrix x,
                                       NumericMatrix bbox,
                                       double bw,
                                       int kernel = 0, //0 gauss 1 prod epa
                                       int loo = 0){ // leave one out
  
  NumericVector out(x.nrow());
  
  int i,j,k;
  
  double dim = x.ncol();
  double w;
  double a, b;
  double (*kern)(double d) = &gausskern;
  if(kernel == 1) {
    kern = &epakern;
  }
  NumericVector blen(dim);
  for(k = 0; k < dim; k++) blen(k) = bbox(1,k) - bbox(0,k);
  
  double bwd = pow(bw, dim);
  
  for(i=0; i < x.nrow(); i++) {
    for(j=i+loo; j < x.nrow(); j++) {
      w = 1.0/bwd;
      for(k=0; k < dim; k++){
        a = abs(x(i,k)-x(j,k));
        b = min(blen(k)-a, a);
        w *= kern(b/bw);
      }
      out(i) += w;
      if(i<j) out(j) += w;
      
    }
  }
  
  return out;
  
}




/* Kernel intensity estimation at non-data points */

// [[Rcpp::export]]
NumericVector int_at_anywhere_c(NumericMatrix x,
                                NumericMatrix loc,
                                NumericMatrix bbox,
                                double bw,
                                int kernel = 0, //0 gauss, 1 prod epa
                                int border = 1){ // 0 none, 1 local, 2 global 3 toroidal
  
  NumericVector out(loc.nrow());
  
  int i,j,k;
  
  double dim = x.ncol();
  double w, wbi, wbj;
  double a, b;
  
  double (*bcf)(double ai, double bi) = bfun_none;
  if(border == 1){
    bcf = bfun_local;
  }
  else if(border == 2){
    bcf = bfun_global;
  }
  else if(border == 3){
    return int_at_anywhere_toroidal_c(x, loc, bbox, bw, kernel);
  }
  
  double (*kernCumu)(double d) = &gaussCumu;
  double (*kern)(double d) = &gausskern;
  double (*kernRange)(double d) = &gaussRange;
  if(kernel == 1) {
    kernCumu = &epaCumu;
    kern = &epakern;
    kernRange = &epaRange;
  }
  

  double bwd = pow(bw, dim);
  
  for(i=0; i < loc.nrow(); i++) {
    wbi = 1;
    for(k = 0; k < dim; k++) {
      a = bbox(0,k)-loc(i,k);
      b = bbox(1,k)-loc(i,k);
      wbi *= kernCumu(b/bw) - kernCumu(a/bw);
    }
    for(j = 0; j < x.nrow(); j++) {
      wbj = 1;
      w = 1/bwd;
      for(k=0; k < dim; k++){
        a = loc(i,k)-x(j,k);
        w *= kern(a/bw);
        a = bbox(0,k)-x(j,k);
        b = bbox(1,k)-x(j,k);
        wbj *= kernCumu(b/bw) - kernCumu(a/bw);
      }
      out(i) += w / bcf(wbi, wbj);
    }
  }
  return out;
}



// toroidal case separately

NumericVector int_at_anywhere_toroidal_c(NumericMatrix x,
                                         NumericMatrix loc,
                                         NumericMatrix bbox,
                                         double bw,
                                         int kernel){
  
  NumericVector out(loc.nrow());
  
  int i,j,k;
  
  double dim = x.ncol();
  double w;
  double a, b;
  // kernel 0 = gaussian
  double (*kern)(double d) = &gausskern;
  if(kernel == 1) {
    kern = &epakern;
  }
  
  NumericVector blen(dim);
  for(k = 0; k < dim; k++) blen(k) = bbox(1,k) - bbox(0,k);
  
  double bwd = pow(bw, dim); // kernel h^d
  
  for(i = 0; i < loc.nrow(); i++) {
    for(j = 0; j < x.nrow(); j++) {
      w = 1.0/bwd;
      for(k=0; k < dim; k++){
        a = abs(loc(i,k)-x(j,k));
        b = min(blen(k)-a, a);
        w *= kern(b/bw);
      }
      out(i) += w;
    }
  }
  
  return out;
  
}





/* Old versions

 
 // [[Rcpp::export]]
 NumericVector intensity_at_points_c(NumericMatrix x, 
                                     NumericMatrix bbox, 
double bw_r,
int b=1){

NumericVector out(x.nrow());

int i,j,k;

double dim = x.ncol();
double d, w;

NumericVector epaw = epa_integral(x, bbox, bw_r, b);

//return epaw;

for(i=0; i < x.nrow()-1; i++) {
for(j=i+1; j < x.nrow(); j++) {
d = 0;
for(k=0; k < dim; k++) d += pow( x(i,k)-x(j,k), 2);
d = sqrt(d);
w = kernel_epa(d, bw_r);
if(w >0) {
out(i) += w/epaw(i);
out(j) += w/epaw(j);
}
}
}

return out;

}


 
// [[Rcpp::export]]
NumericVector intensity_at_other_points_c(NumericMatrix x, 
                                NumericMatrix other,
                                NumericMatrix bbox, 
                                double bw_r,
                                int n=0){
  
  NumericVector out(other.nrow());
  
  NumericVector epaw = epa_integral(other, bbox, bw_r, n);
  
  int i,j,k;
  
  double dim = x.ncol();
  double d, w;
  
  for(i=0; i < x.nrow(); i++) {
    for(j=0; j < other.nrow(); j++) {
      d = 0;
      for(k=0; k < dim; k++) d += pow( x(i,k)-other(j,k), 2);
      d = sqrt(d);
      w = kernel_epa(d, bw_r);
      if(w>0) out(j) += w/epaw(j);
    }
  }
  return out;
}
*/ // eof old versions

 
