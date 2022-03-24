/* Rosenberg's intensities
 * 
 * + helpers
 * 
 */

#include <Rcpp.h>
using namespace Rcpp;

#define TWOPI 6.283185

List bbox2bbquad(NumericMatrix bb) {
  NumericMatrix vert(4,2);
  vert(0,_) = bb(0,_);
  vert(1,0) = bb(1,0); vert(1,1) = bb(0,1);
  vert(2,0) = bb(0,0); vert(2,1) = bb(1,1);
  vert(3,0) = bb(1,0); vert(3,1) = bb(1,1);
  NumericVector idx(4); idx(0)=2-1; idx(1)=4-1; idx(2)=3-1;idx(3)=1-1;
  return List::create(Named("vertices") = vert, Named("idx") = idx);
}

NumericMatrix bbquad_planes(List bq) {
  
  int i;
  double ux, uy,l, nx, ny;
  
  NumericMatrix M(4,4);
  NumericVector idx = bq["idx"];
  NumericMatrix vert = bq["vertices"];
  
  M(0,1) = -1;
  M(1,0) =  1;
  
  idx.push_back(idx(0)); // close the polygon
  
  for(i = 0; i < idx.length()-1; i++) {
    // normal
    ux = vert(idx(i+1),0) - vert(idx(i),0);
    uy = vert(idx(i+1),1) - vert(idx(i),1);
    nx = -uy;
    ny =  ux;
    l = sqrt( nx*nx + ny*ny );
    M(0, i) = nx/l;
    M(1, i) = ny/l;
    // point
    M(2, i) = (vert(idx(i+1),0)+vert(idx(i),0))/2.0;
    M(3, i) = (vert(idx(i+1),1)+vert(idx(i),1))/2.0;
  }
  
  return M;
}

// [[Rcpp::export]]
NumericMatrix line_hit_planes(NumericVector line, NumericMatrix planes) {
  // line = (x1, x2, u1, u2) , x loc, u direction unit
  // planes cols=(n1, n2, p1, p2) n normal unit, p loc on plane (2d line)
  
  NumericMatrix hit(planes.ncol(), 2);
  
  double dx, dy, s, t;
  int i;
  
  for(i = 0; i < planes.ncol(); i++) {
    dx = planes(2,i) - line(0);
    dy = planes(3,i) - line(1);
    s  = line(2) * planes(0,i) + line(3) * planes(1,i);
    t  = 0.9999999 * ( dx * planes(0,i) + dy * planes(1,i) )/s; // rounding errors? Fix properly? How?
    hit(i,0) = line(0) + t * line(2);
    hit(i,1) = line(1) + t * line(3);
  }
  
  return hit;
}


// [[Rcpp::export]]
NumericMatrix c_rosenberg_intensities(NumericMatrix x, NumericMatrix bbox, int steps = 180) {
  NumericMatrix counts(x.nrow(), steps);
  NumericMatrix areas( x.nrow(), steps);
  //NumericMatrix sectorhits(steps, 1);
  // 1st, compute local directional intensities
  int i, j, k, n;
  double ang = 0, dx = 0, dy = 0, r2 = 0;
  
  double width =  M_PI / (double)steps;
  double halfwidth  = width / 2.0;
  
  double base =  atan(M_PI/(2.0 * steps)); // sector triangle base, for approximation
  
  NumericMatrix hits;
  NumericVector line(4);
  NumericMatrix planes = bbquad_planes( bbox2bbquad(bbox) );
  
  for(i = 0; i < x.nrow(); i++) {
    for(j = i + 1; j < x.nrow(); j++) {
      dx = x(i, 0) - x(j, 0);
      dy = x(i, 1) - x(j, 1);
      ang = atan2(dy, dx) + halfwidth;
      if(ang < 0 )  ang = ang + M_PI;
      if(ang >= M_PI)  ang = ang - M_PI;
      k = (int)(ang / width); // slot
      counts(i, k)+= 1;
      counts(j, k)+= 1;
      //sectorhits(k,0)+=1;
    }
    //
    // approximation of the double-cone cap window area
    // use rosenberg's triangle approxiamtion
    line(0) = x(i,0);
    line(1) = x(i,1);
    for(k = 0; k < steps; k++) {
      ang = halfwidth * (2 * k); // central-direction of the cone, upper half
      line(2) = cos(ang);
      line(3) = sin(ang);
      // first get hit points
      n = 0; // check
      hits = line_hit_planes(line, planes);
      areas(i,k) = 0;
      for(j = 0; j < hits.nrow(); j++) {
        // check which actually hit the box
        if( hits(j,0) <  bbox(0,0) || hits(j,0) > bbox(1,0) || hits(j,1) < bbox(0,1) || hits(j,1) > bbox(1,1)) {}
        else{
          // we should have two of these
          n++;
          dx = hits(j,0) - x(i,0);
          dy = hits(j,1) - x(i,1);
          r2 = dx*dx + dy*dy; // distance from point to edge, ^2
          areas(i, k) += base * r2; // area of the isosceles triangle
        }
      }
      if(n != 2) Rprintf("%i,%f:%i\n", i, ang, n); // rounding error or something happenedd...
      // normalise the count
      counts(i,k) /= areas(i,k);
    }
  }
  //return bbquad_planes( bbox2bbquad(bbox) );
  //return  bbox2bbquad(bbox) ;
  //return areas;
  //return sectorhits;
  return counts;
}
