#include "helpers.h"
double min(double a, double b){
  if(a < b ) return a;
  return b;
}
double kernel_epa(double r, double h){
  if(r > h) return 0;
  //return 0.5/h;
  // Epanechnikov
  return 0.75 * (1-r*r/(h*h))/h;
}

double abs(double x){
  if(x < 0) return -x;
  return x;
}
