#include "epa_integral.h"

// [[Rcpp::export]]
NumericVector epa_integral(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n=1){
  // choose approximation
  if(n==1) return box_integral(x, bbox, bw, n);
  if(n==2) return epa_integral_biased(x, bbox, bw, 0);
  return epa_integral_grid(x, bbox, bw, n);
}
  

// [[Rcpp::export]]
NumericVector epa_integral_biased(NumericMatrix x,
                                NumericMatrix bbox,
                                double bw,
                                int n=21) {
  
  int dim = x.ncol();
  // this is bad in 3D
  // biased box integral
  NumericVector out(x.nrow());
  int i,k;
  NumericMatrix lims(2, dim);
  double s, a;
  double h = bw, w=1;
  
  if(dim==2) w = h;
  
  for(i=0; i < x.nrow(); i++) {
    // integration box
    for(k=0; k < dim; k++) {
      lims.at(0,k) = fmin(h, x.at(i,k)-bbox.at(0,k));
      lims.at(1,k) = fmin(h, bbox.at(1,k)-x.at(i,k));
    }
    s = 1;
    for(k=0; k < dim; k++) s*= (lims.at(0,k)+lims.at(1,k));
    s *= 0.75/h;
    
    if(dim==2){
      a  = (pow(lims.at(0,0),3) + pow(lims.at(1,0), 3) ) * (lims.at(0,1)+lims.at(1,1));
      a += (lims.at(0,0)+lims.at(1,0)) * (pow(lims.at(0,1),3) + pow(lims.at(1,1), 3) );
    }
    else{
      a  = (pow(lims.at(0,0), 3) + pow(lims.at(1,0), 3) ) * (lims.at(0,1)+lims.at(1,1)) *
           (lims.at(0,2) + lims.at(1,2));
      a += (pow(lims.at(0,1), 3) + pow(lims.at(1,1), 3) ) * (lims.at(0,0)+lims.at(1,0)) *
           (lims.at(0,2) + lims.at(1,2));
      a += (pow(lims.at(0,2), 3) + pow(lims.at(1,2), 3) ) * (lims.at(0,1)+lims.at(1,1)) *
           (lims.at(0,0) + lims.at(1,0));
    }
    s  -= a / (4*pow(h, 3));
    out.at(i) = s/w;
  }
  return out;
}



/*--------------------------------------------------*/
/* Fine grid approximation  */

// [[Rcpp::export]]
NumericVector epa_integral_grid(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n) {
  
  // fine grid approximation
  NumericVector out(x.nrow());
  int i,j,k;
  int dim = x.ncol();
  double c, s, d;
  bool inside;
  // make sure n not even
  if(n%2 == 0) n+=1;
  // must be n > 1
  if(n==1) n = 3;
  
  // make a mask
  NumericVector vec(n);
  int m = (n-1)/2;
  double dx = bw/m;
  for(k=0; k < n; k++) {
    vec.at(k) = (double)(k-m) * dx;
  }
  
  // we need all n^dim combinations
  NumericMatrix mask(pow(n, dim), dim+1);
  int idx;
  s=0;
  if(dim==2){
    for(i=0; i < n; i++)
      for(j=0; j < n; j++){
        idx = i*n+j;
        mask.at(idx,0) = vec.at(i);
        mask.at(idx,1) = vec.at(j);
        d = sqrt(vec.at(i)*vec.at(i) + vec.at(j)*vec.at(j));
        mask.at(idx,2) = kernel_epa(d, bw); // mask value
        s+=mask.at(idx,2);
      }
  }
  else if(dim==3) {
    for(i=0; i < n; i++)
      for(j=0; j < n; j++)
        for(k=0; k < n; k++){
          idx = i*n*n + j*n + k;
          mask.at(idx,0) = vec.at(i);
          mask.at(idx,1) = vec.at(j);
          mask.at(idx,2) = vec.at(k);
          d = sqrt(vec.at(i)*vec.at(i)+vec.at(j)*vec.at(j)+vec.at(k)*vec.at(k));
          mask.at(idx,3) = kernel_epa(d, bw);
          s+=mask.at(idx,3);
      }
  }
  else {
    Rprintf("Dimension not 2 or 3, will fail.\n");
  }
  double tot = s;
  
  // normalize so that the mask sums to 1
  for(i=0; i < mask.nrow(); i++) mask.at(i,dim) /= tot;
  
  // run integration: Sum all values of mask inside the window
  for(i=0; i < x.nrow(); i++) {
    s = 0.0;
    for(j=0; j < mask.nrow(); j++) {
      inside = true;
      for(k=0; k < dim; k++) {
        c = x.at(i,k) + mask.at(j,k); // shifted mask point's k'th coordinate
        if( c < bbox.at(0,k) | c > bbox.at(1,k) ) {
          inside = false;
          break; // no need to check the rest of dimensions
        }
      }
      if(inside){
        s += mask.at(j, dim);
      }
    }
    out.at(i) = s;
  }
  return out;
}