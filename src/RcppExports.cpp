// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// c_angles
List c_angles(NumericMatrix x, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_angles(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_angles(x, from, to));
    return __result;
END_RCPP
}
// c_angles_in_a_cone
List c_angles_in_a_cone(NumericMatrix x, NumericVector unit, double theta, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_angles_in_a_cone(SEXP xSEXP, SEXP unitSEXP, SEXP thetaSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unit(unitSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_angles_in_a_cone(x, unit, theta, from, to));
    return __result;
END_RCPP
}
// c_anisotropic_pcf
NumericVector c_anisotropic_pcf(NumericMatrix x, List theta, NumericVector r, NumericVector h, NumericMatrix bbox, int correction);
RcppExport SEXP Kdirectional_c_anisotropic_pcf(SEXP xSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    __result = Rcpp::wrap(c_anisotropic_pcf(x, theta, r, h, bbox, correction));
    return __result;
END_RCPP
}
// c_anisotropic_unit_pcf
List c_anisotropic_unit_pcf(NumericMatrix x, NumericMatrix r, NumericVector h, NumericMatrix bbox, int correction);
RcppExport SEXP Kdirectional_c_anisotropic_unit_pcf(SEXP xSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    __result = Rcpp::wrap(c_anisotropic_unit_pcf(x, r, h, bbox, correction));
    return __result;
END_RCPP
}
// c_cutgeom
List c_cutgeom(NumericMatrix x, List nlist, double r);
RcppExport SEXP Kdirectional_c_cutgeom(SEXP xSEXP, SEXP nlistSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type nlist(nlistSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    __result = Rcpp::wrap(c_cutgeom(x, nlist, r));
    return __result;
END_RCPP
}
// c_directed_geom
List c_directed_geom(NumericMatrix x, NumericVector u, double theta, double r, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_directed_geom(SEXP xSEXP, SEXP uSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_directed_geom(x, u, theta, r, from, to));
    return __result;
END_RCPP
}
// c_directed_geom_by_cut
List c_directed_geom_by_cut(NumericMatrix x, NumericVector u, List pregraph, double theta, double r, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_directed_geom_by_cut(SEXP xSEXP, SEXP uSEXP, SEXP pregraphSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< List >::type pregraph(pregraphSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_directed_geom_by_cut(x, u, pregraph, theta, r, from, to));
    return __result;
END_RCPP
}
// c_geom
List c_geom(NumericMatrix x, IntegerVector from, IntegerVector to, double r);
RcppExport SEXP Kdirectional_c_geom(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    __result = Rcpp::wrap(c_geom(x, from, to, r));
    return __result;
END_RCPP
}
// c_knnangles
List c_knnangles(NumericMatrix x, int k, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_knnangles(SEXP xSEXP, SEXP kSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_knnangles(x, k, from, to));
    return __result;
END_RCPP
}
// c_oh_K
NumericVector c_oh_K(NumericMatrix x, List theta, NumericVector r, NumericMatrix bbox);
RcppExport SEXP Kdirectional_c_oh_K(SEXP xSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP bboxSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    __result = Rcpp::wrap(c_oh_K(x, theta, r, bbox));
    return __result;
END_RCPP
}
// c_p_of_KStest
double c_p_of_KStest(int n, double d);
RcppExport SEXP Kdirectional_c_p_of_KStest(SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    __result = Rcpp::wrap(c_p_of_KStest(n, d));
    return __result;
END_RCPP
}
// c_pairwise_dist_angle
NumericMatrix c_pairwise_dist_angle(NumericMatrix x);
RcppExport SEXP Kdirectional_c_pairwise_dist_angle(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(c_pairwise_dist_angle(x));
    return __result;
END_RCPP
}
// c_pairwise_dist_angle_subset
NumericMatrix c_pairwise_dist_angle_subset(NumericMatrix x, IntegerVector from, IntegerVector to);
RcppExport SEXP Kdirectional_c_pairwise_dist_angle_subset(SEXP xSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP);
    __result = Rcpp::wrap(c_pairwise_dist_angle_subset(x, from, to));
    return __result;
END_RCPP
}
// c_sector_pcf
NumericVector c_sector_pcf(NumericMatrix x, NumericVector u, double theta, NumericVector r, double h, NumericMatrix bbox, int correction);
RcppExport SEXP Kdirectional_c_sector_pcf(SEXP xSEXP, SEXP uSEXP, SEXP thetaSEXP, SEXP rSEXP, SEXP hSEXP, SEXP bboxSEXP, SEXP correctionSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    Rcpp::traits::input_parameter< int >::type correction(correctionSEXP);
    __result = Rcpp::wrap(c_sector_pcf(x, u, theta, r, h, bbox, correction));
    return __result;
END_RCPP
}
// c_translation_weights
NumericVector c_translation_weights(NumericMatrix x, NumericMatrix bbox);
RcppExport SEXP Kdirectional_c_translation_weights(SEXP xSEXP, SEXP bboxSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bbox(bboxSEXP);
    __result = Rcpp::wrap(c_translation_weights(x, bbox));
    return __result;
END_RCPP
}