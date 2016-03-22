#' Estimate intensity of point pattern at given locations
#' 
#' @param x points, list with $x and $bbox or just a matrix of coordinates.
#' @param loc Matrix of locations to estimate the intensity at
#' @param bw bandwidth, Epanechnikov kernel domain [-bw,bw]
#' @param b Type of border correction estimate.
#' @param ... ignored
#' @details 
#' Uses Epanechnikov kernel smoothing. Border corrections controlled by what b is:
#' 0: 2D exact, not available for 3D atm; 
#' 1: Box rectangle approximation; 
#' 2: 2D Biased box integral, 3D none;
#' >3: fine grid sum with resolution b^dim (slow). 
#' 
#' if b < 0 no edge correction is applied.
#' 
#' @export

intensity_somewhere <- function(x, loc, bw, b=0, ...) {
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
#' Uses Epanechnikov kernel smoothing. Border corrections controlled by what b is:
#' 0: 2D exact, not available for 3D atm; 
#' 1: Box rectangle approximation; 
#' 2: 2D Biased box integral, 3D none;
#' >3: fine grid sum with resolution b^dim (slow). 
#' 
#' if b < 0 no edge correction is applied.
#' 
#' 
#' @export

intensity_at_points <- function(x, bw, b=0, ...) {
  x <- check_pp(x)
  bbox <- x$bbox
  coord <- x$x
  intensity_at_points_c(coord, bbox, bw, b)
}

#' Find the optimal smoothing for intensity estimation kernel width
#' 
#' @param x point pattern
#' @param bw_vector Vector of bandwidth values to optimize over
#' @param keep_min if true, keep the intensity_at_points for the best bw
#' @param ... passed on to \code{intensity_at_points}
#' @details 
#' 
#' Optimize bandwidth $h$ of the Epanechnikov kernel using loss function
#' \deqn{latex}{\sum 1/\hat\lambda_h(x_i) - |W|)^2}
#' 
#' @return 
#' 
#' Vector of losses, or if keep_min=TRUE, list with vector of losses and vector of intensity values with the bw that gives the minimum loss.
#' 
#' @export

intensity_bandwidth_profile <- function(x, bw_vector, keep_min=FALSE, ...){
  x <- check_pp(x)
  V <- bbox_volume(x$bbox)
  err <- NULL
  for(bw in bw_vector) {
      o <- intensity_at_points(x, bw=bw, ...)
      # discrepancy
      err1 <- (sum(1/o)-V)^2
      if(keep_min) if(err1 <= min(err)) i <- o
      err <- c(err, err1)
  }
  if(keep_min) list(err=err, intensity_best=i)
  else err
}