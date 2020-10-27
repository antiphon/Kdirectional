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
                      steps = 180, include = .75) {
  x <- check_pp(x)
  if(ncol(x$x) != 2) stop("Use only for 2D data")
  if( include > 1 || include < 0 ) stop("'include' needs to be between 0 and 1.")
  #
  loc <- x$x
  bb <- x$bbox
  
  # pointwise directional densities
  etas <- c_rosenberg_intensities(loc, bb, steps)
  
  # Wavelet to use
  french_top_hat <- function(v) {
    o <- 0 * v
    a <- abs(v)
    o[a < 1.5] <- -1
    o[a < 0.5] <-  2
    o
  }
  
  # To avoid edge effects:
  step <- pi/steps
  thetai0 <- seq(0, pi-step , by = step) # actual interesting
  thetai <- seq(-max(scales), pi - step + max(scales), by = step)  # buffered
  thetaij <- thetai 
  thetaij[ thetaij < 0]  <- thetaij[thetaij < 0 ]  + pi
  thetaij[ thetaij > pi] <- thetaij[thetaij > pi ] - pi
  # find match
  ij <- match( round(thetaij, 9), round(thetai0, 9)) 
  if(any(is.na(ij))) ij <- sapply(thetaij, function(v) which.min(abs(v-thetai0))) # rounding errors again
  # wrap
  arg <- outer(thetai, theta, "-")
  etass  <-  etas[,ij]
  #arg <- pmin(outer(thetai, theta, "-"), outer(thetai, pi-theta, "-"))
  # 
#  browser()
  
  WM <- simplify2array( lapply( scales, function(bk){ etass %*% french_top_hat( arg / bk ) / bk } )  ) 
  # pointwise "variance"
  #browser()
  Px <- apply(WM^2, c(1,2), mean)
  # Edge correction
  bb_inc <-  t( (t(bb) - colMeans(bb)) * sqrt(include) + colMeans(bb)  ) 
  inc <- loc[,1] > bb_inc[1,1] & loc[,1] < bb_inc[2,1] & loc[,2] > bb_inc[1,2] & loc[,2] < bb_inc[2,2]
  Px <- Px[inc,, drop = FALSE]
  # overall variance
  data.frame( theta = theta, Phat = colMeans(Px) )
}
