#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

#ifndef PP_H_
#define PP_H_

class Pp
{
  int npoints;
  int dim;
  NumericMatrix X;
  NumericMatrix bbox;
  NumericVector boxlen;
  NumericVector marks;
  double (Pp::*dist)(int*, int*);
  std::vector<double> distTriangle;
  std::vector<double> * pdists;

  double distEuclidian(int*, int*);
  double distGreatCircle(int*, int*);
  double distPrecalculated(int*, int*);

  double (Pp::*weight)(int*, int*);
  double weightAll1(int *, int *);
  std::vector<double> weightTriangle;

  std::vector<int> typevec;

public:
  Pp(NumericMatrix, NumericVector, NumericMatrix );
  virtual ~Pp();

  double getCoord(int *i, int *d);
  int    size();
  int    d();
  double getMark(int *);
  double getDist(int *, int *);
  double getAngle(int *, int *);
  double getWeight(int *, int *);
};

#endif /*PP_H_*/
