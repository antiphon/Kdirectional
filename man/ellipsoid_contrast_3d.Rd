% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_ellipsoid_contrasts.R
\name{ellipsoid_contrast_3d}
\alias{ellipsoid_contrast_3d}
\title{Testing contrast for ellipsoids}
\usage{
ellipsoid_contrast_3d(els, ijk = c(3, 1, 2), out = TRUE, z0 = 1.85)
}
\description{
Give this to confint.ellipsoid to compute the contrast for
 testing semi_i - 0.5 * (semi_j + semi_k) = 0
}
\details{
@param els list of ellipsoid estimates
 @param ijk on of c(1,2,3), c(2,1,3) or c(3,1,2). Default: c(3,1,2)
 @param out if >0 Use outlier detection to stabilise (useful with large variances)



The default contrast will compare third semi-axis to the average of first two. This seems to be the best, as the semi-axes are returned in an increasing order by the OLS-algorithm.

The 'out' is used like 'coef' in boxplot.stats, see there for details.
}
