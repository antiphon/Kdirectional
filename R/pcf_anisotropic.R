#' Anisotropic pair correlation function
#' 
#' Estimate the anisotropic pair correlation function (2d and 3d), as defined in f 1991, f. 5.2-5.3.
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param r radius vector at which to estimate
#' @param theta vector or list of vectors for for angles, see details.
#' @param h widths of epanechnicov kernels, vector of two values, for ranges and angles.
#' @param f If h not given, use h=c( f/lambda^(1/dim), h=f*pi). Same as 'stoyan' in spatstat's pcf.
#' @param correction "none" or translation. Translation only for rectangle box.
#' @param n_dir Angle grid resolution (in case theta not given)
#' 
#' The antipode symmetry makes it necessary to compute only on one half of the circle/sphere.
#' 
#' @useDynLib Kdirectional
#' @export

pcf_anisotropic <- function(x, r, theta, h, f=0.15, correction="translation", n_dir=7) {
  x <- check_pp(x)
  bbox <- as.matrix(x$bbox)
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  lambda <- nrow(x$x)/bbox_volume(bbox)
  
  if(missing(theta)){
    dirs <- check_directions(n_dir=n_dir ,dim=dim)
    theta <- dirs$theta
    if(dim==3) theta <- append(theta, theta)
  }
  if(!is.list(theta)) stop("Theta should be a list.")
  if(any(sapply(theta, abs)>pi)) stop("Theta should be in range [0, pi]")
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  if(missing(h)) {
    h <- c( f/lambda^(1/dim), f*pi )  
  }
  if(length(h) < 2) h <- rep(h, 2)
  # correction
  correction_i <- pmatch(correction, c("none", "translation"))
  #' start:
  xc <- as.matrix(x$x)
  #'
  units <- theta_2_unit(theta)
  #'
  directions <- NULL 
  for(ri in r) directions <- rbind(directions, ri*units)
  #
  res <- c_anisotropic_unit_pcf(xc, directions, h, bbox, correction_i)
  #'  correction
  grid <- expand.grid(append(list(r=r), theta))
  
  rho <- res[[1]]
  #' pcf
  g <- rho/lambda^2
  g_m <- matrix(g, byrow=T, nrow=length(r))
  #'
  res <- list(est=g_m, r=r, theta=theta, directions=directions, correction=correction, h=h, rho=res[[1]], dim=dim, counts=res[[2]])
  class(res) <- "pcf_anisotropic"
  res
}
