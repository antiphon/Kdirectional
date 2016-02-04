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
  
  double w =  PI * h2;
  double a, b, c, d;
  double vh, ch, dh;
  //
  for(i=0; i < x.nrow(); i++) {
    // possibly hitting the edge:
    a = fmin(h, x.at(i,0)-bbox.at(0,0));
    b = fmin(h, bbox.at(1,0)-x.at(i,0));
    c = fmin(h, x.at(i,1)-bbox.at(0,1));
    d = fmin(h, bbox.at(1,1)-x.at(i,1)); 
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
    out.at(i) = s/w;
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
        mask.at(idx,2) = kernel_box(d, bw); // mask value
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
