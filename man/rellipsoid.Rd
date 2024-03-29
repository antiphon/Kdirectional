% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_ellipsoid_rellipsoid.R
\name{rellipsoid}
\alias{rellipsoid}
\title{Sample points Roughly Uniformly on the surface of an ellipsoid, 2D or 3D}
\usage{
rellipsoid(n, axes = c(1, 1, 1), center, noise.sd = 0, R = NULL, pieces = 1000)
}
\arguments{
\item{axes}{vector giving the semi-axis lengths}

\item{center}{center vector, c(x,y) or c(x,y,z)}

\item{noise.sd}{noise sd, 0=no noise}

\item{R}{rotation matrix}

\item{pieces}{Number of pieces to divide the surface into. See details.}
}
\description{
Sample points Roughly Uniformly on the surface of an ellipsoid, 2D or 3D
}
\details{
Sampling is done by first sampling uniformly on the surface of a
sphere and then transforming with the ellipsoid axes. Corrections:

* 2D: The sampling weights are set to piecewise linear approximation of the ellipse surface integral, 
the resolution of which is controlled by 'pieces'.
* 3D: No weighting at the moment, so very biased for elongated ellipsoids.

Additional nD Gaussian noise can be added by setting `noise.sd`>0.
}
