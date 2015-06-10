
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List c_geom(NumericMatrix x, IntegerVector from, IntegerVector to, double r) {
    int i, j, l, f, t;
    int n = x.nrow();
    int dim = x.ncol();
    double d;
    List nlist(n);
    std::vector<int> neighs;
    for(i=0; i < from.length(); i++) {
      neighs.clear();
      f = from[i]-1;
      for(j=0; j < to.length(); j++) {
        t = to[j]-1;
        if(t!=f) {
          d=0;
          for(l=0; l < dim; l++)  d += pow(x(f,l)-x(t,l), 2);
          d = sqrt(d);
          if(d < r) neighs.push_back(t+1); 
        }
      }
      nlist(f) = wrap( neighs );
    }
    return nlist;
}
