% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_ellipsoid_confint.R
\name{confint.ellipsoid}
\alias{confint.ellipsoid}
\title{Confidence interval for ellipsoid parameters using simulation}
\usage{
\method{confint}{ellipsoid}(
  object,
  parm,
  level,
  fun,
  nsim = 1000,
  probs = c(0.025, 0.975),
  Asolve.eps = 0.01,
  ...
)
}
\arguments{
\item{object}{Fitted ellipsoid}

\item{fun}{Optional contrast function. See ellipse_contrast_2d for an example.}

\item{probs}{The quantiles for the CI. Default is c(0.025, 0.975)}

\item{Asolve.eps}{Inflate the A matrix diagonals by a factor of eps, to avoid numerical problems with (near) singular parameter sets.}

\item{...}{Passed on to sample_ellipse_beta}

\item{tol}{tolerance for absolute  deviation in the ||beta||=1 constraint.}
}
\description{
Assuming normality of errors.
}
