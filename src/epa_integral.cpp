#include "epa_integral.h"

// [[Rcpp::export]]
NumericVector epa_integral(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n=0){
  // choose approximation
  if(n < 0){
    NumericVector xx(x.nrow(), 1.0);
    return xx;
  }
  if(n==0) return epa_integral_2d(x, bbox, bw, n);
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
      lims(0,k) = fmin(h, x(i,k)-bbox(0,k));
      lims(1,k) = fmin(h, bbox(1,k)-x(i,k));
    }
    s = 1;
    for(k=0; k < dim; k++) s*= (lims(0,k)+lims(1,k));
    s *= 0.75/h;
    
    if(dim==2){
      a  = (pow(lims(0,0),3) + pow(lims(1,0), 3) ) * (lims(0,1)+lims(1,1));
      a += (lims(0,0)+lims(1,0)) * (pow(lims(0,1),3) + pow(lims(1,1), 3) );
    }
    else{
      a  = (pow(lims(0,0), 3) + pow(lims(1,0), 3) ) * (lims(0,1)+lims(1,1)) *
           (lims(0,2) + lims(1,2));
      a += (pow(lims(0,1), 3) + pow(lims(1,1), 3) ) * (lims(0,0)+lims(1,0)) *
           (lims(0,2) + lims(1,2));
      a += (pow(lims(0,2), 3) + pow(lims(1,2), 3) ) * (lims(0,1)+lims(1,1)) *
           (lims(0,0) + lims(1,0));
    }
    s  -= a / (4*pow(h, 3));
    out(i) = s/w;
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
    vec(k) = (double)(k-m) * dx;
  }
  
  // we need all n^dim combinations
  NumericMatrix mask(pow(n, dim), dim+1);
  int idx;
  s=0;
  if(dim==2){
    for(i=0; i < n; i++)
      for(j=0; j < n; j++){
        idx = i*n+j;
        mask(idx,0) = vec(i);
        mask(idx,1) = vec(j);
        d = sqrt(vec(i)*vec(i) + vec(j)*vec(j));
        mask(idx,2) = kernel_epa(d, bw); // mask value
        s+=mask(idx,2);
      }
  }
  else if(dim==3) {
    for(i=0; i < n; i++)
      for(j=0; j < n; j++)
        for(k=0; k < n; k++){
          idx = i*n*n + j*n + k;
          mask(idx,0) = vec(i);
          mask(idx,1) = vec(j);
          mask(idx,2) = vec(k);
          d = sqrt(vec(i)*vec(i)+vec(j)*vec(j)+vec(k)*vec(k));
          mask(idx,3) = kernel_epa(d, bw);
          s+=mask(idx,3);
      }
  }
  else {
    Rprintf("Dimension not 2 or 3, will fail.\n");
  }
  double tot = s;
  
  // normalize so that the mask sums to 1
  for(i=0; i < mask.nrow(); i++) mask(i,dim) /= tot;
  
  // run integration: Sum all values of mask inside the window
  for(i=0; i < x.nrow(); i++) {
    s = 0.0;
    for(j=0; j < mask.nrow(); j++) {
      inside = true;
      for(k=0; k < dim; k++) {
        c = x(i,k) + mask(j,k); // shifted mask point's k'th coordinate
        if( c < bbox(0,k) | c > bbox(1,k) ) {
          inside = false;
          break; // no need to check the rest of dimensions
        }
      }
      if(inside){
        s += mask(j, dim);
      }
    }
    out(i) = s;
  }
  return out;
}



// for the intergration
double Iy(double y){
  return (y*sqrt(1-y*y)*(5-2*y*y)+3*asin(y))/8.0;
}

double IIy(double a, double b){
  return Iy(b)-Iy(a);
}



// [[Rcpp::export]]
NumericVector epa_integral_2d(NumericMatrix x,
                              NumericMatrix bbox,
                              double bw,
                              int n=0) {
  // fine grid?
  if(n > 3) return epa_integral_grid(x, bbox, bw, n);
  
  NumericVector out(x.nrow());
  int i;
  int dim = x.ncol();
  double s, h = bw;
  //double h2 = h*h;
  
  if(dim >2) {
    Rprintf("3D analytical border correction for Epanechnikov not implemented; use other border methods.\n");
    for(i=0; i < out.size();i++) out(1)=1;
    return out; 
  }
  
  // 2D only.
  
  double w =  PI/2.0;
  double a, b, c, d;
  double vh, ch, dh;
  //
  for(i=0; i < x.nrow(); i++) {
    // possibly hitting the edge:
    a = fmin(h, x(i,0)-bbox(0,0))/h;
    b = fmin(h, bbox(1,0)-x(i,0))/h;
    c = fmin(h, x(i,1)-bbox(0,1))/h;
    d = fmin(h, bbox(1,1)-x(i,1))/h; 
    // not hitting the edge?
    if(a==1 & b==1 & c==1 & d==1) s=w;
    else{
      s=0;
      // hit edge
      // part 1
      vh = sqrt(1-a*a);
      ch = fmin(vh, c);
      dh = fmin(vh, d);
      s+= a*(ch+dh)-a*a*a*(ch+dh)/3.0-a*(ch*ch*ch+dh*dh*dh)/3.0 + (2.0/3.0)*(IIy(-c,-ch)+IIy(dh,d));
      // part 2
      vh = sqrt(1-b*b);
      ch = fmin(vh, c);
      dh = fmin(vh, d);
      s+= b*(ch+dh)-b*b*b*(ch+dh)/3.0-b*(ch*ch*ch+dh*dh*dh)/3.0 + (2.0/3.0)*(IIy(-c,-ch)+IIy(dh,d));
    }
    out(i) = s/w;
  }
  return out;
}

