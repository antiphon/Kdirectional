#' Rosenberg Wavelet Transform
#' 
#' Estimate the wavelet variance suggested by Rosenberg 2004. 2D only, rectangular box only.
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param theta directions to estimate the wavelet variance at
#' @param scales angular scales to average over, in radians. 
#' @param steps divide the half-circle into this many sectors around each point (default: 180)
#' @param include Only estimate the variance from points in the central region, area fraction given by 'include'. default: 0.5. 
#' 
#' @details 
#' Compute the estimator given in Rosenberg 2004.
#'  
#' @references 
#' Rosenberg M S, Wavelet analysis for detecting anisotropy in point patterns, J of Veg Sci 2004
#' 
#' @export

rosenberg <- function(x, theta = seq(0, pi, l = 50), scales = 1:45 / 180 * pi, 
                      steps = 180, include = .5) {
  x <- check_pp(x)
  if(ncol(x$x) != 2) stop("Use only for 2D data")
  if( include > 1 || include < 0 ) stop("'include' needs to be between 0 and 1.")
  #
  loc <- x$x
  bb <- x$bbox
  # pointwise directional densities
  etas <- c_rosenberg_intensities(loc, bb, steps)
  step <- pi/steps
  thetai <- seq(0, pi-step, by = step) + step/2
  # 
  
  # Wavelet to use
  french_hat <- function(v) {
    o <- 0 * v
    a <- abs(v)
    o[a < 3/2] <- -1
    o[a < 1/2] <- 2
    o
  }
  
  # cyclic distance?
  #arg <- pmin(outer(thetai, theta, "-"), outer(thetai, pi-theta, "-"))
  #browser()
  arg <- outer(thetai, theta, "-")
  
  bk <- scales
  WM <- simplify2array( lapply( scales, function(bk){ etas %*% french_hat( arg / bk ) / bk } )  )
  # pointwise "variance"
  Px <- apply(WM^2, c(1,2), mean)
  # Edge correction
  bb_inc <-  t( (t(bb) - colMeans(bb)) * sqrt(include) + colMeans(bb)  ) 
  inc <- loc[,1] > bb_inc[1,1] & loc[,1] < bb_inc[2,1] & loc[,2] > bb_inc[1,2] & loc[,2] < bb_inc[2,2]
  Px <- Px[inc,, drop = FALSE]
  # overall variance
  data.frame( theta = theta, Phat = colMeans(Px) )
}
