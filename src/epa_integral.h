#include <Rcpp.h>
#include "box_integral.h"
#include "helpers.h"
#include <vector>

using namespace Rcpp;

#ifndef EPAI_H_
#define EPAI_H_

NumericVector epa_integral(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n);

NumericVector epa_integral_biased(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n);

NumericVector epa_integral_grid(NumericMatrix x,
                           NumericMatrix bbox,
                           double bw,
                           int n);

#endif