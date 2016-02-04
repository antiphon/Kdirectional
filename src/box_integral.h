#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

#ifndef BOXI_H_
#define BOXI_H_

NumericVector box_integral(NumericMatrix x,
                              NumericMatrix bbox, double, int);

NumericVector box_integral_grid(NumericMatrix x,
                           NumericMatrix bbox, double, int);

#endif