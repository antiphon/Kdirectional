#include "Pp.h"

/********************************************************************************************/
Pp::~Pp()
{
}
/********************************************************************************************/
Pp::Pp(NumericMatrix coord, NumericVector m, NumericMatrix bb)
{
  X = coord;
  bbox = bb;
  marks = m;
  npoints = X.nrow();
  dim = X.ncol();
  //
  int i;
  boxlen = NumericVector(dim);
  for(i=0; i < dim; i++) boxlen(i) = bbox(1,i) - bbox(0,i);
  // distance
  dist = &Pp::distEuclidian;
}
/********************************************************************************************/
int Pp::size(){
  return npoints;
}
int Pp::d() {
  return dim;
}
/********************************************************************************************/
double Pp::getCoord(int *i, int *d){
  return X(*i, *d);
}
double Pp::getMark(int *i){
  return marks(*i);
}

/********************************************************************************************/
double Pp::distEuclidian(int *i, int *j)
{
  if(*i==*j) return 0.0;
  if(*i>*j) return distEuclidian(j, i);
  double s=0;
  for(int k=0; k < dim; k++) s+=pow(getCoord(i, &k)-getCoord(j, &k) , 2.0);
  return sqrt(s);
}
/********************************************************************************************/
double Pp::getDist(int *i, int *j)
{
  return (this->*dist)(i,j);
}

/********************************************************************************************/
double Pp::getAngle(int *i, int *j)
{
  return 0;
}


/********************************************************************************************/
double Pp::getWeight(int *i, int *j)
{
  double wij = 1;
  for(int k=0; k < dim; k++) wij *= boxlen(k) - fabs(X(*i,k)-X(*j,k));
  return wij;
}
