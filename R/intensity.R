#' Estimate intensity of point pattern at given locations
#' 
#' @param x points, list with $x and $bbox or just a matrix of coordinates.
#' @param loc Matrix of locations to estimate the intensity at
#' @param bw bandwidth, Epanechnikov kernel domain [-bw,bw]
#' @param b Type of border correction estimate.
#' @param ... ignored
#' @details 
#' Uses Epanechnikov kernel smoothing. Border correction is approximative, try "b=0,1,2,3". Default 0 (works atm only in 2D, in 3D use grid with n~30)
#' 0: Exact integral 1: Box kernel approximation, 2: Biased epa integral (not working in 3d) >3: fine grid sum with resolution b, very slow.
#' 
#' @export

intensity_somewhere <- function(x, loc, bw, b=1, ...) {
  x <- check_pp(x)
  bbox <- x$bbox
  coord <- x$x
  if(ncol(loc)!=ncol(coord)) stop("pattern and locations of different dimension.")
  intensity_at_other_points_c(coord, loc, bbox, bw, b)
}

#' Estimate intensity of point pattern at the points of the pattern
#' 
#' @param x points, list with $x and $bbox or just a matrix of coordinates.
#' @param bw bandwidth, Epanechnikov kernel domain [-bw,bw]
#' @param b Type of border correction estimate.
#' @param ... ignored
#' @details 
#' Uses Epanechnikov kernel smoothing. Border correction is approximative, try "b=1,2,3".
#' 1: Box rectangle approximation, 2: Biased box integral (not working in 3d) >3: fine grid sum with resolution b, very slow.
#' 
#' @export

intensity_at_points <- function(x, bw, b=1, ...) {
  x <- check_pp(x)
  bbox <- x$bbox
  coord <- x$x
  intensity_at_points_c(coord, bbox, bw, b)
}

#' Find the optimal smoothing for intensity estimation kernel width
#' 
#' @param x point pattern
#' @param bw_vector Vector of bandwidth values to optimize over
#' @param ... passed on to \code{intensity_at_points}
#' @details 
#' 
#' Optimize bandwidth $h$ of the Epanechnikov kernel using loss function
#' \deqn{latex}{\sum 1/\hat\lambda_h(x_i) - |W|)^2}
#' 
#' @return 
#' 
#' Vector of losses.
#' 
#' @export

intensity_bandwidth_profile <- function(x, bw_vector, ...){
  x <- check_pp(x)
  V <- bbox_volume(x$bbox)
  err <- NULL
  for(bw in bw_vector) {
      o <- intensity_at_points(x, bw=bw, ...)
      # discrepancy
      err <- c(err, (sum(1/o)-V)^2)
  }
  err
}