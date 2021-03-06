% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcf_anin_cylinder.R
\name{pcf_anin_cylinder}
\alias{pcf_anin_cylinder}
\title{Inhomogeneous anisotropic pcf function, cylinder version}
\usage{
pcf_anin_cylinder(
  x,
  u,
  epsilon,
  r,
  lambda = NULL,
  lambda_h,
  r_h,
  stoyan = 0.15,
  renormalise = TRUE,
  border = 1,
  aspect = 1/3,
  ...
)
}
\arguments{
\item{x}{pp, list with $x~coordinates $bbox~bounding box}

\item{u}{unit vector(s) of direction, as row vectors. Default: x and y axis.}

\item{epsilon}{The cylinder half-width}

\item{r}{radius vector at which to evaluate the function}

\item{lambda}{optional vector of intensity estimates at points}

\item{lambda_h}{if lambda missing, use this bandwidth in a kernel estimate of lambda(x)}

\item{r_h}{smoothing for range dimension, epanechnikov kernel}

\item{stoyan}{If r_h not given, use r_h=stoyan/lambda^(1/dim). Same as 'stoyan' in spatstat's pcf.}

\item{renormalise}{See details.}

\item{border}{Use translation correction? Default=1, yes. Only for cuboidal windows.}

\item{aspect}{If given, instead of using a fixed halfwidth (epsilon) take the halfwidth to be r/(2*aspect). Default : 1/3}

\item{...}{passed on to e.g. \link{intensity_at_points}}
}
\value{
Returns a dataframe.
}
\description{
Estimate a cylinder-pcf function for second order reweighted ("inhomogeneous") pattern.
}
\details{
Computes a second order reweighted version of the cylinder-pcf, defined as the function to integrate in range over [0,R] to get the cylinder-K(R) function.

Lambda(x) at points can be given, 
or else it will be estimated using Epanechnikov kernel smoothing. See 

If 'renormalise=TRUE', we normalise the lambda estimate so that sum(1/lambda(x))=|W|. This corresponds in \code{spatstat}'s \code{Kinhom} to setting 'normpower=2'.
}
