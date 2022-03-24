#include <Rcpp.h>
#include "box_integral.h" 


// integrand for 2D 
double Iy(double y, double bw){
  return 0.5 * ( y * sqrt(bw*bw-y*y)+ bw *bw * asin(y/bw));
}



// [[Rcpp::export]]
NumericVector box_integral(NumericMatrix x,
                                NumericMatrix bbox,
                                double bw,
                                int n=0) {
  // fine grid?
  if(n > 3) return box_integral_grid(x, bbox, bw, n);
  
  NumericVector out(x.nrow());
  int i;
  int dim = x.ncol();
  double s, h = bw;
  double h2 = h*h;
  
  if(dim >2) {
    Rprintf("3D not implemented.\n");
    return out; 
  }
  
  // 2D only.
  
  double w =  M_PI * h2;
  double a, b, c, d;
  double vh, ch, dh;
  //
  for(i=0; i < x.nrow(); i++) {
    // possibly hitting the edge:
    a = fmin(h, x(i,0)-bbox(0,0));
    b = fmin(h, bbox(1,0)-x(i,0));
    c = fmin(h, x(i,1)-bbox(0,1));
    d = fmin(h, bbox(1,1)-x(i,1)); 
    // not hitting the edge?
    if(a==h & b==h & c==h & d==h) s=w;
    else{
      s=0;
      // hit edge
      // part 1
      vh = sqrt(h2-a*a);
      ch = fmin(vh, c);
      dh = fmin(vh, d);
      s+= a * (ch+dh) + Iy(-ch, h)-Iy(-c,h) + Iy(d,h)-Iy(dh,h);
      // part 2
      vh = sqrt(h2-b*b);
      ch = fmin(vh, c);
      dh = fmin(vh, d);
      s+= b * (ch+dh) + Iy(-ch, h)-Iy(-c,h) + Iy(d,h)-Iy(dh,h);
    }
    out(i) = s/w;
  }
  return out;
}



/* Fine grid approximation  */

// [[Rcpp::export]]
NumericVector box_integral_grid(NumericMatrix x,
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
        mask(idx,2) = kernel_box(d, bw); // mask value
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
