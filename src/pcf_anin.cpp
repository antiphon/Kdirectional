#include <Rcpp.h>
#include "helpers.h"
#include "Pp.h"

using namespace Rcpp;


double indicator_function(double r, double a) { // less than a?
  if(r < a) return 1;
  return 0;
}



// [[Rcpp::export]]
NumericMatrix pcf_anin_fry_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, 
                             NumericVector r, double bw, NumericMatrix directions,
                             int border=1) {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);
  
  int i,j,l, ui, ri;
  double d, d2, d3, w;
  double rmax = r(nr-1) + bw;
  double bw2 = bw * bw;
  
  for(i=0; i < pp.size()-1 ; i++) {
    for(j= i + 1; j < pp.size(); j++) {
      d = pp.getDist(&i, &j); 
      if(d < rmax) {
        w =  1.0/(pp.getMark(&i) * pp.getMark(&j) ); // intensities
        if(border==1) w /= pp.getWeight(&i, &j) ; // translation correction
        for(ri=0; ri < nr; ri++){
          if(d < r(ri)+bw && r(ri) < d+bw){ //
            for(ui=0; ui < ndir; ui++){
              // fry to target point distance. Use antipodal estimation.
              d2 = 0;
              d3 = 0;
              for(l=0; l < dim; l++) {
                d2 += pow( pp.getCoord(&j, &l)-pp.getCoord(&i, &l) - directions(ui, l) * r(ri) , 2);
                d3 += pow( pp.getCoord(&i, &l)-pp.getCoord(&j, &l) - directions(ui, l) * r(ri) , 2);
              }
              d2 = d2/bw2;
              d3 = d3/bw2;
              if(d2 < 1) out(ri, ui) += (1-d2) * w; // normalise the kernel in R
              if(d3 < 1) out(ri, ui) += (1-d3) * w; // normalise the kernel in R
            }
          }
        }
      }
    }
  }
  
  return out;
}



// [[Rcpp::export]]
NumericMatrix pcf_anin_conical_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, 
                           NumericVector r, double r_h, 
                           NumericMatrix directions,
                           double epsilon, // box kernel width
                           int border=1, 
                           int divisor_i = 0, 
                           int ang_kernel = 0) // 0: indicator function (not kernel), 1:kernel epa. 
  {
  
  Pp pp(coord, lambda, bbox); 
  
  int nr = r.size();
  int dim = directions.ncol();
  int ndir = directions.nrow();
  NumericMatrix out(nr, ndir);

  int i,j,l, ui, ri;
  double d, w, dot, ang;
  double ka, kd;
  double rmax = r(nr-1)+r_h;
  
  double (*ang_kern)(double d, double h) = indicator_function;
  if(ang_kernel>0) ang_kern = &kernel_epa;
  
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
          ang = fmin(ang, M_PI-ang);
          // check cone inclusion or kernel range
          if(ang < epsilon){
            ka = ang_kern(ang, epsilon);
            for(ri=0; ri < nr; ri++){
              div  = pow(divisorf(r(ri), d), dim-1);
              kd =  kernel_epa(fabs(d-r(ri)), r_h);
              out(ri, ui) += ka * kd * w / div;
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
                                    NumericVector epsilon, // half width of cylinder (radius)
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
          { // inside the infinite cylinder
            for(ri=0; ri < nr; ri++){
              if(dist < epsilon(ri)){
                k = kernel_epa(fabs(d-r(ri)), r_h);
                out(ri, ui) += k*w / div;
              }
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
          ang = fmin(ang, M_PI-ang);
          if(ang < epsilon){
            for(ri=0; ri < nr; ri++){
              k = kernel_epa(fabs(d-r(ri)), r_h);
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