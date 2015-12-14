#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;


// double kernel_epa(double d, double bw){
//   if(d > bw) return 0;
//   return (1-(d/bw)*(d/bw))*3.0/(4*bw);
// }

int max(int i, int j){
  if(i > j) return i;
  return j;
}

int min(int i, int j){
  if(i < j) return i;
  return j;
}


// [[Rcpp::export]]
List anisotropic_markcor_c(NumericMatrix coord, NumericVector marks, NumericMatrix bbox, 
                        NumericVector r, NumericMatrix directions,
                        double bw_r, double bw_a, Function f) {
  
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
  double rmax = r(nr-1);
  double rmaxc = rmax + iw*dr;
  
  double d, w, ka, kr, dot, ang;
  for(i=0; i < pp.size()-1 ; i++) {
    fvec =  f(marks(i), marks[seq(i+1, pp.size()-1)]  ); // semi vectorized, memory requirement ~n
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); // distance(i,j)
      if(d < rmaxc) {
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          w = pp.getWeight(&i, &j);
          rc =  floor(d/rmax * nr);
          ka = kernel_epa(ang, bw_a);
          if(ka > 0){
            for(ri = max(0, rc-iw); ri < min(rc+iw, nr); ri++) {
              kr = kernel_epa(fabs(d-r(ri)), bw_r);
              upper(ri,ui) += fvec(j-i-1) * kr * ka/ w;
              lower(ri,ui) += kr * ka / w;
            }
          }
        }
      }
    }
  }
  return List::create(upper=upper, lower=lower);
}


// [[Rcpp::export]]
List anisotropic_markcor_c_d(NumericMatrix coord, NumericVector marks, NumericMatrix bbox, 
                           NumericVector r, NumericMatrix directions,
                           double bw_r, double bw_a, Function f) {
  
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
  double rmax = r(nr-1);
  double rmaxc = rmax + ((double)iw)*dr;
  
  double d, w, ka, kr, dot, ang, pd;
  for(i=0; i < pp.size()-1 ; i++) {
    fvec =  f(marks(i), marks[seq(i+1, pp.size()-1)]  ); // semi vectorized, memory requirement ~n
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j);
      if(d < rmaxc) {
        pd = pow(d, dim-1); 
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  
            dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          if(ang < bw_a){
            w = pp.getWeight(&i, &j) * pd;
            rc =  floor(d/rmax * nr);
            ka = kernel_epa(ang, bw_a);
            for(ri = max(0, rc-iw); ri < min(rc+iw, nr); ri++) {
              kr = kernel_epa(fabs(d-r(ri)), bw_r);
              upper(ri,ui) += fvec(j-i-1) * kr * ka/ w;
              lower(ri,ui) += kr * ka / w;
            }
          }
        }
      }
    }
  }
  return List::create(upper=upper, lower=lower);
}
