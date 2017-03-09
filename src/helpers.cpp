#include "helpers.h"
double min(double a, double b){
  if(a < b ) return a;
  return b;
}
double max(double a, double b){
  if(a > b ) return a;
  return b;
}
double kernel_epa(double r, double h){
  if(r > h) return 0;
  return 0.75 * (1-r*r/(h*h))/h;
}

// unit bandwidth, no normalising constant, for any dimension then.
double kernel_epa_un(double r){
  if(r > 1) return 0;
  return (1-r*r);
}


double kernel_box(double r, double h){
  if(r > h) return 0;
  return 1/(2*h);
}


double abs(double x){
  if(x < 0) return -x;
  return x;
}


int max(int i, int j){
  if(i > j) return i;
  return j;
}

int min(int i, int j){
  if(i < j) return i;
  return j;
}

double divisor_r(double r, double d){ return r; }
double divisor_d(double r, double d){return d;};// return max(d,r); }
