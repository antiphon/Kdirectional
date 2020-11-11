#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

double gkern(double d) {
  return exp(-0.5 * d*d);
}

// [[Rcpp::export]]
NumericMatrix Kest_gaussian_c(NumericMatrix coord, 
                              NumericVector lambda, 
                              NumericMatrix bbox, 
                              NumericVector r, 
                              NumericMatrix directions,
                              double kappa, 
                              int border=1) {
  
  Pp pp(coord, lambda, bbox); // mark with lambda
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, w, dot, dot2, v, a;
  double rmax, normaliser;
  double quan = 2.0; // corresponds to about 95% reach in x-direction
  // the compression diagonal matrix. We have Sigma = TT^T, with T=RC. u defines R; kappa defines C.
  NumericVector C(dim);
  C(0) = 1; // sqrt of diagonal of the covariance of the axis-oriented kernel, major axis longer
  for(i=1; i < dim; i++) C(i) = pow(kappa, 1.0/(dim-1.0)); // in case 1D or.
  //
  normaliser = 1.0 / ( sqrt(pow(2*PI, dim))  ); // basic gaussian constant
  // max relevant range is in direction u 
  rmax =  C(0) * r(nr-1) * 2  ; 
  //
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w = 1/ (pp.getMark(&i) * pp.getMark(&j) ); // should be lambdas
        if(border==1) w/= pp.getWeight(&i, &j) ; // translation correction...
        w *= normaliser;
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  {
            dot +=  (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l); // length of fry vector in direction u
          }
          dot2 = sqrt(d*d - dot*dot); // length of fry in any direction perp to u
          if(dot < rmax && dot2 < rmax){
            for(ri=0; ri < nr; ri++){
              a = r(ri) / quan;
              v =       gkern( dot / (a*C(0)) ) / (a*C(0)) ; // major axis 
              v *= pow( gkern(dot2 / (a*C(1)) ) / (a*C(1)) , dim-1); // minor axes
              v *= w;
              out(ri, ui) += v;
              //Rprintf("d1:%f, d2:%f: v:%f\n", dot, dot2, v);
            }
          }
        }
      }
    }
  }
  return out;
}
