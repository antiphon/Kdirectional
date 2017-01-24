#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pcf_anin_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, 
                           NumericVector r, double r_h, NumericMatrix directions,
                           double epsilon, int border=1, int divisor_i = 0) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);

  int i,j,l, ui, ri;
  double d, w, dot, ang;
  double k;
  double rmax = r(nr-1)+r_h;
  
  double (*divisorf)(double d, double r) = &divisor_d;
  if(divisor_i == 0) divisorf = &divisor_r;
  double div;
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w =  1.0/(pp.getMark(&i) * pp.getMark(&j) );
        if(border==1) w /= pp.getWeight(&i, &j) ;
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          if(ang < epsilon){
            for(ri=0; ri < nr; ri++){
              div  = pow(divisorf(r(ri), d), dim-1);
              k = kernel_epa(abs(d-r(ri)), r_h);
              out(ri, ui) += k*w / div;
            }
          }
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix pcf_anin_cylindrical_c(NumericMatrix coord, 
                                     NumericVector lambda, 
                                     NumericMatrix bbox, 
                                     NumericVector r, 
                                     double r_h, 
                                     NumericMatrix directions,
                                    double epsilon, // half width of cylinder (radius)
                                    int border=1, int divisor_i = 0) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  int i,j,l, ui, ri;
  double d, w, dot, dist;
  double k;
  double rmax = r(nr-1)+r_h;
  double dif[dim];
  //double (*divisorf)(double d, double r) = &divisor_d; // not applicable
  //if(divisor_i == 0) divisorf = &divisor_r;
  double div = 1; 
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w =  1.0/(pp.getMark(&i) * pp.getMark(&j) );
        if(border==1) w /= pp.getWeight(&i, &j) ;
        for(ui=0; ui < ndir; ui++){
          // distance from the cross section hitting i.e. projection plane going through 0
          dot = 0;
          for(l=0; l < dim; l++)  {
            dif[l] = pp.getCoord(&j,&l)-pp.getCoord(&i,&l);
            dot += dif[l] * directions(ui, l);
          }
          d = abs(dot); 
          // distance from the line in direction u
          dist  = 0;
          for(l = 0; l < dim; l++) dist += pow( dif[l] - dot*directions(ui,l), 2);
          dist = sqrt(dist);
          if(dist < epsilon){ // inside the infinite cylinder
            for(ri=0; ri < nr; ri++){
              k = kernel_epa(abs(d-r(ri)), r_h);
              out(ri, ui) += k*w / div;
            }
          }
        }
      }
    }
  }
  return out;
}




/// obsolete below

/*

// [[Rcpp::export]]
NumericMatrix pcf_anin_c_d(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, 
                         NumericVector r, double r_h, NumericMatrix directions,
                         double epsilon, int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, w, dot, ang;
  double k, pd;
  double rmax = r(nr-1)+r_h;
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        pd = pow(d, dim-1); 
        w = 1/ (pd * pp.getMark(&i) * pp.getMark(&j) );
        if(border==1) w/=pp.getWeight(&i, &j) ;
        for(ui=0; ui < ndir; ui++){
          dot = 0;
          for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
          ang = acos(dot/d);
          ang = fmin(ang, PI-ang);
          if(ang < epsilon){
            for(ri=0; ri < nr; ri++){
              k = kernel_epa(abs(d-r(ri)), r_h);
              out(ri, ui) += k*w;
            }
          }
        }
      }
    }
  }
  return out;
}


*/