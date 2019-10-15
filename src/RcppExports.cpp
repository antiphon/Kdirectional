// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Kest_anin_c
NumericMatrix Kest_anin_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double epsilon, int border);
RcppExport SEXP _Kdirectional_Kest_anin_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(Kest_anin_c(coord, lambda, bbox, r, directions, epsilon, border));
    return rcpp_result_gen;
END_RCPP
}
// Kest_anin_border_c
NumericMatrix Kest_anin_border_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector bdist, NumericVector r, NumericMatrix directions, double epsilon, int border);
RcppExport SEXP _Kdirectional_Kest_anin_border_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP bdistSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bdist(bdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(Kest_anin_border_c(coord, lambda, bbox, bdist, r, directions, epsilon, border));
    return rcpp_result_gen;
END_RCPP
}
// Kest_anin_cylinder_c
NumericMatrix Kest_anin_cylinder_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double epsilon, int border);
RcppExport SEXP _Kdirectional_Kest_anin_cylinder_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(Kest_anin_cylinder_c(coord, lambda, bbox, r, directions, epsilon, border));
    return rcpp_result_gen;
END_RCPP
}
// Kest_anin_cylinder_border_c
NumericMatrix Kest_anin_cylinder_border_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector bdist, NumericVector r, NumericMatrix directions, double epsilon, int border);
RcppExport SEXP _Kdirectional_Kest_anin_cylinder_border_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP bdistSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bdist(bdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(Kest_anin_cylinder_border_c(coord, lambda, bbox, bdist, r, directions, epsilon, border));
    return rcpp_result_gen;
END_RCPP
}
// Kest_gaussian_c
NumericMatrix Kest_gaussian_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double kappa, int border);
RcppExport SEXP _Kdirectional_Kest_gaussian_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP kappaSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(Kest_gaussian_c(coord, lambda, bbox, r, directions, kappa, border));
    return rcpp_result_gen;
END_RCPP
}
// c_angles
List c_angles(NumericMatrix x, IntegerVector from, IntegerVector to);
RcppExport SEXP _Kdirectional_c_angles(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(c_angles(x, from, to));
    return rcpp_result_gen;
END_RCPP
}
// c_angles_in_a_cone
List c_angles_in_a_cone(NumericMatrix x, NumericVector unit, double theta, IntegerVector from, IntegerVector to, bool antipodal);
RcppExport SEXP _Kdirectional_c_angles_in_a_cone(SEXP xSEXP, SEXP unitSEXP, SEXP thetaSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP antipodalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unit(unitSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    Rcpp::traits::input_parameter< bool >::type antipodal(antipodalSEXP);
    rcpp_result_gen = Rcpp::wrap(c_angles_in_a_cone(x, unit, theta, from, to, antipodal));
    return rcpp_result_gen;
END_RCPP
}
// anisotropic_markcor_c
List anisotropic_markcor_c(NumericMatrix coord, NumericVector marks, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double bw_r, double bw_a, Function f);
RcppExport SEXP _Kdirectional_anisotropic_markcor_c(SEXP coordSEXP, SEXP marksSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP bw_rSEXP, SEXP bw_aSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type marks(marksSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type bw_r(bw_rSEXP);
    Rcpp::traits::input_parameter< double >::type bw_a(bw_aSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(anisotropic_markcor_c(coord, marks, bbox, r, directions, bw_r, bw_a, f));
    return rcpp_result_gen;
END_RCPP
}
// anisotropic_markcor_c_d
List anisotropic_markcor_c_d(NumericMatrix coord, NumericVector marks, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double bw_r, double bw_a, Function f);
RcppExport SEXP _Kdirectional_anisotropic_markcor_c_d(SEXP coordSEXP, SEXP marksSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP bw_rSEXP, SEXP bw_aSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type marks(marksSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type bw_r(bw_rSEXP);
    Rcpp::traits::input_parameter< double >::type bw_a(bw_aSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(anisotropic_markcor_c_d(coord, marks, bbox, r, directions, bw_r, bw_a, f));
    return rcpp_result_gen;
END_RCPP
}
// c_anisotropic_pcf
NumericVector c_anisotropic_pcf(NumericMatrix x, List theta, NumericVector r, NumericVector h, NumericMatrix bbox, int correction);
RcppExport SEXP _Kdirectional_c_anisotropic_pcf(SEXP xSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    rcpp_result_gen = Rcpp::wrap(c_anisotropic_pcf(x, theta, r, h, bbox, correction));
    return rcpp_result_gen;
END_RCPP
}
// c_anisotropic_unit_pcf
List c_anisotropic_unit_pcf(NumericMatrix x, NumericMatrix r, NumericVector h, NumericMatrix bbox, int correction);
RcppExport SEXP _Kdirectional_c_anisotropic_unit_pcf(SEXP xSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    rcpp_result_gen = Rcpp::wrap(c_anisotropic_unit_pcf(x, r, h, bbox, correction));
    return rcpp_result_gen;
END_RCPP
}
// box_integral
NumericVector box_integral(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_box_integral(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(box_integral(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// box_integral_grid
NumericVector box_integral_grid(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_box_integral_grid(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(box_integral_grid(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// c_cutgeom
List c_cutgeom(NumericMatrix x, List nlist, double r);
RcppExport SEXP _Kdirectional_c_cutgeom(SEXP xSEXP, SEXP nlistSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type nlist(nlistSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cutgeom(x, nlist, r));
    return rcpp_result_gen;
END_RCPP
}
// c_directed_geom
List c_directed_geom(NumericMatrix x, NumericVector u, double theta, double r, IntegerVector from, IntegerVector to);
RcppExport SEXP _Kdirectional_c_directed_geom(SEXP xSEXP, SEXP uSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(c_directed_geom(x, u, theta, r, from, to));
    return rcpp_result_gen;
END_RCPP
}
// c_directed_geom_by_cut
List c_directed_geom_by_cut(NumericMatrix x, NumericVector u, List pregraph, double theta, double r, IntegerVector from, IntegerVector to);
RcppExport SEXP _Kdirectional_c_directed_geom_by_cut(SEXP xSEXP, SEXP uSEXP, SEXP pregraphSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< List >::type pregraph(pregraphSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(c_directed_geom_by_cut(x, u, pregraph, theta, r, from, to));
    return rcpp_result_gen;
END_RCPP
}
// epa_integral
NumericVector epa_integral(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_epa_integral(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(epa_integral(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// epa_integral_biased
NumericVector epa_integral_biased(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_epa_integral_biased(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(epa_integral_biased(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// epa_integral_grid
NumericVector epa_integral_grid(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_epa_integral_grid(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(epa_integral_grid(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// epa_integral_2d
NumericVector epa_integral_2d(NumericMatrix x, NumericMatrix bbox, double bw, int n);
RcppExport SEXP _Kdirectional_epa_integral_2d(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(epa_integral_2d(x, bbox, bw, n));
    return rcpp_result_gen;
END_RCPP
}
// c_geom
List c_geom(NumericMatrix x, IntegerVector from, IntegerVector to, double r);
RcppExport SEXP _Kdirectional_c_geom(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(c_geom(x, from, to, r));
    return rcpp_result_gen;
END_RCPP
}
// int_bw_c
NumericMatrix int_bw_c(NumericMatrix x, NumericVector h, int kernel);
RcppExport SEXP _Kdirectional_int_bw_c(SEXP xSEXP, SEXP hSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(int_bw_c(x, h, kernel));
    return rcpp_result_gen;
END_RCPP
}
// int_at_points_c
NumericVector int_at_points_c(NumericMatrix x, NumericMatrix bbox, double bw, int kernel, int border, int loo);
RcppExport SEXP _Kdirectional_int_at_points_c(SEXP xSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP borderSEXP, SEXP looSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    Rcpp::traits::input_parameter< int >::type loo(looSEXP);
    rcpp_result_gen = Rcpp::wrap(int_at_points_c(x, bbox, bw, kernel, border, loo));
    return rcpp_result_gen;
END_RCPP
}
// int_at_anywhere_c
NumericVector int_at_anywhere_c(NumericMatrix x, NumericMatrix loc, NumericMatrix bbox, double bw, int kernel, int border);
RcppExport SEXP _Kdirectional_int_at_anywhere_c(SEXP xSEXP, SEXP locSEXP, SEXP bboxSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loc(locSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(int_at_anywhere_c(x, loc, bbox, bw, kernel, border));
    return rcpp_result_gen;
END_RCPP
}
// c_knnangles
List c_knnangles(NumericMatrix x, int k, IntegerVector from, IntegerVector to);
RcppExport SEXP _Kdirectional_c_knnangles(SEXP xSEXP, SEXP kSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(c_knnangles(x, k, from, to));
    return rcpp_result_gen;
END_RCPP
}
// markcor_anin_c
List markcor_anin_c(NumericMatrix coord, NumericVector marks, NumericVector lambda, NumericMatrix bbox, NumericVector r, NumericMatrix directions, double bw_r, double bw_a, int border, Function f, int divisor_i);
RcppExport SEXP _Kdirectional_markcor_anin_c(SEXP coordSEXP, SEXP marksSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP directionsSEXP, SEXP bw_rSEXP, SEXP bw_aSEXP, SEXP borderSEXP, SEXP fSEXP, SEXP divisor_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type marks(marksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type bw_r(bw_rSEXP);
    Rcpp::traits::input_parameter< double >::type bw_a(bw_aSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type divisor_i(divisor_iSEXP);
    rcpp_result_gen = Rcpp::wrap(markcor_anin_c(coord, marks, lambda, bbox, r, directions, bw_r, bw_a, border, f, divisor_i));
    return rcpp_result_gen;
END_RCPP
}
// c_oh_K
NumericVector c_oh_K(NumericMatrix x, List theta, NumericVector r, NumericMatrix bbox);
RcppExport SEXP _Kdirectional_c_oh_K(SEXP xSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP bboxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    rcpp_result_gen = Rcpp::wrap(c_oh_K(x, theta, r, bbox));
    return rcpp_result_gen;
END_RCPP
}
// c_p_of_KStest
double c_p_of_KStest(int n, double d);
RcppExport SEXP _Kdirectional_c_p_of_KStest(SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(c_p_of_KStest(n, d));
    return rcpp_result_gen;
END_RCPP
}
// c_pairwise_dist_angle
NumericMatrix c_pairwise_dist_angle(NumericMatrix x);
RcppExport SEXP _Kdirectional_c_pairwise_dist_angle(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pairwise_dist_angle(x));
    return rcpp_result_gen;
END_RCPP
}
// c_pairwise_dist_angle_subset
NumericMatrix c_pairwise_dist_angle_subset(NumericMatrix x, IntegerVector from, IntegerVector to);
RcppExport SEXP _Kdirectional_c_pairwise_dist_angle_subset(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pairwise_dist_angle_subset(x, from, to));
    return rcpp_result_gen;
END_RCPP
}
// pcf_anin_fry_c
NumericMatrix pcf_anin_fry_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, double bw, NumericMatrix directions, int border);
RcppExport SEXP _Kdirectional_pcf_anin_fry_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP bwSEXP, SEXP directionsSEXP, SEXP borderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    rcpp_result_gen = Rcpp::wrap(pcf_anin_fry_c(coord, lambda, bbox, r, bw, directions, border));
    return rcpp_result_gen;
END_RCPP
}
// pcf_anin_conical_c
NumericMatrix pcf_anin_conical_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, double r_h, NumericMatrix directions, double epsilon, int border, int divisor_i, int ang_kernel);
RcppExport SEXP _Kdirectional_pcf_anin_conical_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP r_hSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP, SEXP divisor_iSEXP, SEXP ang_kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type r_h(r_hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    Rcpp::traits::input_parameter< int >::type divisor_i(divisor_iSEXP);
    Rcpp::traits::input_parameter< int >::type ang_kernel(ang_kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(pcf_anin_conical_c(coord, lambda, bbox, r, r_h, directions, epsilon, border, divisor_i, ang_kernel));
    return rcpp_result_gen;
END_RCPP
}
// pcf_anin_cylindrical_c
NumericMatrix pcf_anin_cylindrical_c(NumericMatrix coord, NumericVector lambda, NumericMatrix bbox, NumericVector r, double r_h, NumericMatrix directions, double epsilon, int border, int divisor_i);
RcppExport SEXP _Kdirectional_pcf_anin_cylindrical_c(SEXP coordSEXP, SEXP lambdaSEXP, SEXP bboxSEXP, SEXP rSEXP, SEXP r_hSEXP, SEXP directionsSEXP, SEXP epsilonSEXP, SEXP borderSEXP, SEXP divisor_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type r_h(r_hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type border(borderSEXP);
    Rcpp::traits::input_parameter< int >::type divisor_i(divisor_iSEXP);
    rcpp_result_gen = Rcpp::wrap(pcf_anin_cylindrical_c(coord, lambda, bbox, r, r_h, directions, epsilon, border, divisor_i));
    return rcpp_result_gen;
END_RCPP
}
// c_sector_pcf
NumericVector c_sector_pcf(NumericMatrix x, NumericVector u, double theta, NumericVector r, double h, NumericMatrix bbox, int correction);
RcppExport SEXP _Kdirectional_c_sector_pcf(SEXP xSEXP, SEXP uSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    rcpp_result_gen = Rcpp::wrap(c_sector_pcf(x, u, theta, r, h, bbox, correction));
    return rcpp_result_gen;
END_RCPP
}
// c_translation_weights
NumericVector c_translation_weights(NumericMatrix x, NumericMatrix bbox);
RcppExport SEXP _Kdirectional_c_translation_weights(SEXP xSEXP, SEXP bboxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    rcpp_result_gen = Rcpp::wrap(c_translation_weights(x, bbox));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Kdirectional_Kest_anin_c", (DL_FUNC) &_Kdirectional_Kest_anin_c, 7},
    {"_Kdirectional_Kest_anin_border_c", (DL_FUNC) &_Kdirectional_Kest_anin_border_c, 8},
    {"_Kdirectional_Kest_anin_cylinder_c", (DL_FUNC) &_Kdirectional_Kest_anin_cylinder_c, 7},
    {"_Kdirectional_Kest_anin_cylinder_border_c", (DL_FUNC) &_Kdirectional_Kest_anin_cylinder_border_c, 8},
    {"_Kdirectional_Kest_gaussian_c", (DL_FUNC) &_Kdirectional_Kest_gaussian_c, 7},
    {"_Kdirectional_c_angles", (DL_FUNC) &_Kdirectional_c_angles, 3},
    {"_Kdirectional_c_angles_in_a_cone", (DL_FUNC) &_Kdirectional_c_angles_in_a_cone, 6},
    {"_Kdirectional_anisotropic_markcor_c", (DL_FUNC) &_Kdirectional_anisotropic_markcor_c, 8},
    {"_Kdirectional_anisotropic_markcor_c_d", (DL_FUNC) &_Kdirectional_anisotropic_markcor_c_d, 8},
    {"_Kdirectional_c_anisotropic_pcf", (DL_FUNC) &_Kdirectional_c_anisotropic_pcf, 6},
    {"_Kdirectional_c_anisotropic_unit_pcf", (DL_FUNC) &_Kdirectional_c_anisotropic_unit_pcf, 5},
    {"_Kdirectional_box_integral", (DL_FUNC) &_Kdirectional_box_integral, 4},
    {"_Kdirectional_box_integral_grid", (DL_FUNC) &_Kdirectional_box_integral_grid, 4},
    {"_Kdirectional_c_cutgeom", (DL_FUNC) &_Kdirectional_c_cutgeom, 3},
    {"_Kdirectional_c_directed_geom", (DL_FUNC) &_Kdirectional_c_directed_geom, 6},
    {"_Kdirectional_c_directed_geom_by_cut", (DL_FUNC) &_Kdirectional_c_directed_geom_by_cut, 7},
    {"_Kdirectional_epa_integral", (DL_FUNC) &_Kdirectional_epa_integral, 4},
    {"_Kdirectional_epa_integral_biased", (DL_FUNC) &_Kdirectional_epa_integral_biased, 4},
    {"_Kdirectional_epa_integral_grid", (DL_FUNC) &_Kdirectional_epa_integral_grid, 4},
    {"_Kdirectional_epa_integral_2d", (DL_FUNC) &_Kdirectional_epa_integral_2d, 4},
    {"_Kdirectional_c_geom", (DL_FUNC) &_Kdirectional_c_geom, 4},
    {"_Kdirectional_int_bw_c", (DL_FUNC) &_Kdirectional_int_bw_c, 3},
    {"_Kdirectional_int_at_points_c", (DL_FUNC) &_Kdirectional_int_at_points_c, 6},
    {"_Kdirectional_int_at_anywhere_c", (DL_FUNC) &_Kdirectional_int_at_anywhere_c, 6},
    {"_Kdirectional_c_knnangles", (DL_FUNC) &_Kdirectional_c_knnangles, 4},
    {"_Kdirectional_markcor_anin_c", (DL_FUNC) &_Kdirectional_markcor_anin_c, 11},
    {"_Kdirectional_c_oh_K", (DL_FUNC) &_Kdirectional_c_oh_K, 4},
    {"_Kdirectional_c_p_of_KStest", (DL_FUNC) &_Kdirectional_c_p_of_KStest, 2},
    {"_Kdirectional_c_pairwise_dist_angle", (DL_FUNC) &_Kdirectional_c_pairwise_dist_angle, 1},
    {"_Kdirectional_c_pairwise_dist_angle_subset", (DL_FUNC) &_Kdirectional_c_pairwise_dist_angle_subset, 3},
    {"_Kdirectional_pcf_anin_fry_c", (DL_FUNC) &_Kdirectional_pcf_anin_fry_c, 7},
    {"_Kdirectional_pcf_anin_conical_c", (DL_FUNC) &_Kdirectional_pcf_anin_conical_c, 10},
    {"_Kdirectional_pcf_anin_cylindrical_c", (DL_FUNC) &_Kdirectional_pcf_anin_cylindrical_c, 9},
    {"_Kdirectional_c_sector_pcf", (DL_FUNC) &_Kdirectional_c_sector_pcf, 7},
    {"_Kdirectional_c_translation_weights", (DL_FUNC) &_Kdirectional_c_translation_weights, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Kdirectional(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
