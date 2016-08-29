#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
List markcor_anin_c(NumericMatrix coord, 
                    NumericVector marks, 
                    NumericVector lambda, 
                    NumericMatrix bbox, 
                    NumericVector r, 
                    NumericMatrix directions,
                    double bw_r, double bw_a, 
                    int border, Function f, int divisor_i) {
  Pp pp(coord, marks, bbox); 
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix upper(nr, ndir);
  NumericMatrix lower(nr, ndir);
  
  NumericVector fvec;
  int i, j, l, rc, ri, ui;
  
  double dr = r(1)-r(0);
  
  int iw = 2 + (int)(bw_r/dr);
  double rmax = r(nr-1) + bw_r;
  double rmaxc = rmax + iw*dr;
  
  double d, w, ka, kr, dot, ang, div;
  
  double (*divisorf)(double d, double r) = &divisor_d;
  if(divisor_i == 0) divisorf = &divisor_r;
  
  for(i=0; i < pp.size()-1 ; i++) {
    fvec =  f(marks(i), marks[seq(i+1, pp.size()-1)]  ); // semi vectorized, memory requirement ~n
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); // distance(i,j)
      if(d < rmaxc) {
        w = 1;
        if(border >0) w *= pp.getWeight(&i, &j);
        w *= lambda(i) * lambda(j); // inhomogeneous
        for(ui=0; ui < ndir; ui++){ 
          dot = 0;
          for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          rc =  floor(d/rmax * nr);
          //ka = kernel_epa(ang, bw_a);
          ka = 1;
          //if(ka > 0){
          if(ang < bw_a){
            //for(ri = max(0, rc-iw); ri < min(rc+iw, nr); ri++) {
            for(ri=0; ri < nr; ri++){
              div  = pow(divisorf(r(ri), d), dim-1);
              kr = kernel_epa(fabs(d-r(ri)), bw_r);
              upper(ri,ui) += fvec(j-i-1) * kr * ka / (w * div);
              lower(ri,ui) += kr * ka / (w * div);
            }
          }
        }
      }
    }
  }
  return List::create(upper=upper, lower=lower);
}
