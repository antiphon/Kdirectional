#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Kest_anin_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, 
                           NumericVector r, NumericMatrix directions,
                           double epsilon, int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);

  int i,j,l, ui, ri;
  double d, w, dot, ang;
  double rmax = r.at(nr-1);
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w = 1/ (pp.getMark(&i) * pp.getMark(&j) );
        if(border==1) w/=pp.getWeight(&i, &j) ;
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          if(ang < epsilon){
            for(ri=0; ri < nr; ri++){
              if(d < r.at(ri)) out.at(ri, ui) += w;
            }
          }
        }
      }
    }
  }
  return out;
}
