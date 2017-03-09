#' Second order orientation density
#' 
#' As in Illian et al 2008, Ch. 4.5.3.
#' 
#' @param x point pattern
#' @param r1 lower limit of the ranges
#' @param r2 upper limit of the ranges
#' @param alpha angles to estimate the density at, on range [0, pi]
#' @param bw Bandwidth (half-width) for the Epanechnicov kernel.
#' 
#' @export

fry_orientation_density <- function(x, r1, r2, alpha, bw = 0.35) {
  x <- check_pp(x)
  loc <- x$x
  bbox <- x$bbox
  n <- nrow(loc)
  dim <- ncol(bbox)
  if(dim != 2 ) stop("only 2D implemented at the moment.")
  
  M <- pi # maximum angle
  if(missing(r1) | missing(r2)) stop("Ranges r1 < r2 required.")
  if(missing(alpha)) {
    alpha <- seq(0, M, l = 150)
  }
  
  # Lets go:
  fry <- fry_points(x, double = FALSE, border = FALSE)
  ok_r <- with(fry, fry_r < r2 & fry_r >= r1)
  unit_ok <- fry$fry_units[ok_r,, drop=FALSE]
  
  # Border correction
  w <- translation_weights(x)
  w_ok <- w[ok_r]
  
  #
  kern <- function(a) {
    v <- 0.75 * (1 - (a/bw)^2)/bw
    v[a < -bw] <- 0
    v[a > bw] <- 0
    v
  }
  
  # angles
  angs <- unit_2_angle(unit_ok)
  # this maps the alpha range automatically to [0,pi]
  angs[angs < 0] <- angs[angs < 0] + pi
  
  # Estimate the density
  est <- vapply(alpha, function(a) {
    d <- abs(angs - a)
    # toroidal wrap over [0,pi]
    d <- pmin(d, pi - d)
    sum( kern(d) / w_ok)
  },  1.0  )
  #
  # Normalising constant:
  fest <- est/sum(1/w_ok)
  # 
  # Done
  out <- data.frame(alpha = alpha, density = fest)
  class(out) <- c("angle_function", is(out))  
  out
}


