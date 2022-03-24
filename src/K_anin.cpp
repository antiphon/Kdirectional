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
  double rmax = r(nr-1);
  
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
          ang = fmin(ang, M_PI-ang);
          if(ang < epsilon){
            for(ri=0; ri < nr; ri++){
              if(d < r(ri)) out(ri, ui) += w;
            }
          }
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix Kest_anin_border_c(NumericMatrix coord, NumericVector lambda, 
                                 NumericMatrix bbox, NumericVector bdist, 
                                 NumericVector r, NumericMatrix directions,
                                 double epsilon, int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, w, dot, ang;
  double rmax = r(nr-1);
  
  for(i=0; i < pp.size(); i++) {
    for(j= 0; j < pp.size(); j++) {
      if(i!=j){
        d = pp.getDist(&i, &j); 
        if(d < rmax) {
          w = 1 / (pp.getMark(&i) * pp.getMark(&j) );
          for(ui=0; ui < ndir; ui++){
            dot = 0;
            for(l=0; l < dim; l++)  dot += (pp.getCoord(&j,&l)-pp.getCoord(&i,&l)) * directions(ui, l);
            ang = acos(dot/d);
            ang = fmin(ang, M_PI-ang);
            if(ang < epsilon){
              for(ri=0; ri < nr; ri++){
                if(d < r(ri) && (bdist(i) > r(ri) | border==0) ) out(ri, ui) += w;
              }
            }
          }
        }
      }
    }
  }
  return out;
}





// [[Rcpp::export]]
NumericMatrix Kest_anin_cylinder_c(NumericMatrix coord, 
                                   NumericVector lambda,
                                   NumericMatrix bbox, 
                                   NumericVector r, NumericMatrix directions,
                                   NumericVector epsilon, int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, w, dot, dist;
  double rmax = r(nr-1);
  double dif[dim];
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w = 1/ (pp.getMark(&i) * pp.getMark(&j) );
        if(border==1) w/=pp.getWeight(&i, &j) ;
        for(ui=0; ui < ndir; ui++){
          // distance from the cross section hitting i.e. projection plane going through 0
          dot = 0;
          for(l=0; l < dim; l++)  {
            dif[l] = pp.getCoord(&j,&l)-pp.getCoord(&i,&l);
            dot += dif[l] * directions(ui, l);
          }
          d = fabs(dot); 
          // distance from the line in direction u
          dist  = 0;
          for(l = 0; l < dim; l++) dist += pow( dif[l] - dot*directions(ui,l), 2);
          dist = sqrt(dist);
          //if(dist < epsilon){ // inside the infinite cylinder
          {
            for(ri=0; ri < nr; ri++){
              if(d < r(ri) & dist < epsilon(ri)) out(ri, ui) += w;
            }
          }
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix Kest_anin_cylinder_border_c(NumericMatrix coord, NumericVector lambda, 
                                 NumericMatrix bbox, NumericVector bdist, 
                                 NumericVector r, NumericMatrix directions,
                                 NumericVector epsilon, int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, w, dot, dist;
  double rmax = r(nr-1);
  double dif[dim];
  
  for(i=0; i < pp.size(); i++) {
    for(j= 0; j < pp.size(); j++) {
      if(i!=j){
        d = pp.getDist(&i, &j); 
        if(d < rmax) {
          w = 1 / (pp.getMark(&i) * pp.getMark(&j) );
          for(ui=0; ui < ndir; ui++){
            // distance from the cross section hitting i.e. projection plane going through 0
            dot = 0;
            for(l=0; l < dim; l++)  {
              dif[l] = pp.getCoord(&j,&l)-pp.getCoord(&i,&l);
              dot += dif[l] * directions(ui, l);
            }
            d = fabs(dot); 
            // distance from the line in direction u
            dist  = 0;
            for(l = 0; l < dim; l++) dist += pow( dif[l] - dot*directions(ui,l), 2);
            dist = sqrt(dist);
            { 
              for(ri=0; ri < nr; ri++){
                if(dist < epsilon(ri) && d < r(ri) && (bdist(i) > r(ri) | border==0) ) out(ri, ui) += w;
              }
            }
          }
        }
      }
    }
  }
  return out;
}
