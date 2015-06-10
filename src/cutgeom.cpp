
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List c_cutgeom(NumericMatrix x, List nlist, double r) {
    int i, j, k, l;
    int n = x.nrow();
    int dim = x.ncol();
    double d;
    List nlist2(nlist);
    std::vector<int> new_neighs; 
    IntegerVector neighs;
    GenericVector tt;
    
    for(i=0; i < n; i++) {
      new_neighs.clear();
      tt = nlist(i); 
      if(tt.length()){
        neighs = nlist(i);
        for(j=0; j < neighs.length(); j++) {
          k = neighs(j) - 1; // use 1-> indexing
          d=0;
          for(l=0; l < dim; l++)  d += pow(x(i,l)-x(k,l), 2);
          d = sqrt(d);
          if(d < r) new_neighs.push_back(k+1); 
        }
        nlist2(i) = wrap(new_neighs);
      }
    }
    return nlist2;
}
