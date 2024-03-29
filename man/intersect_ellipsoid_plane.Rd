% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_ellipsoid_projections.R
\name{intersect_ellipsoid_plane}
\alias{intersect_ellipsoid_plane}
\title{Intersection of a 3d Ellipsoid and a Plane}
\usage{
intersect_ellipsoid_plane(x, n, q, r)
}
\arguments{
\item{x}{ellipsoid, with $semi_axes and $rot components}

\item{n}{normal of the plane}

\item{q}{a point defining the plane, inside the ellipsoid}

\item{r}{optional, perpendicular to n}
}
\description{
Intersection of a 3d Ellipsoid and a Plane
}
\details{
Note that there is no rotation as the result is in an alternative
basis. The basis is stored 
in the returned ellipse's component $basis3d. Third column is n.
}
\references{
P. P. Klein,"On the Ellipsoid and Plane Intersection Equation," Appl. Math., vol. 3, no. November, pp. 1634-1640, 2012.
}
