#' Sector pair correlation function
#' 
#' This is the same as "sector pair correlation function" in Stoyan 1991.
#' Direction defined with a unit vector. Estimated in a cone of angle theta and different radii.
#' 
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param r Ranges
#' @param theta Direction sector central angles, list of ncol(x$x)-components.
#' @param unit Alternative to theta, a matrix of direction vectors.
#' @param epsilon Direction sector width angles
#' @param h width of epanechnicov kernel
#' @param stoyan if h not given, use h=stoyan/sqrt(lambda)
#' @param correction "none" or translation. Translation only for rectangle box.
#' @param n_dir If unit missing, compute on a grid with resolution n_dir
#' @details
#' Otherwise same as default pcf in spatstat, but now with an angle. 
#' @useDynLib Kdirectional
#' @export

pcf_sector <- function(x, r, theta, unit, epsilon, h, stoyan=0.15, correction="translation", n_dir=12) {
  x <- check_pp(x)
  bbox <- as.matrix(x$bbox)
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  lambda <- nrow(x$x)/prod(sidelengths)
  
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  if(missing(h)) h <-  h <- stoyan/lambda^(1/dim)
  
  #' 
  dirs <- check_directions(theta, epsilon, unit, dim=dim, antipode = FALSE, n_dir = n_dir)
  epsilon <- dirs$epsilon
  # start:
  xc <- as.matrix(x$x)
  
  correction_i <- pmatch(correction, c("none", "translation"))
  #' go for the direction
  dirgrid <- dirs$grid$unit
  v <- NULL
  for(i in 1:nrow(dirgrid)){
    v <- cbind(v, c_sector_pcf(xc, dirgrid[i,], epsilon, r, h, bbox, correction_i))
  }
  sector <- if(dim==2)  2*epsilon * r else  2*pi * r^2 * (1-cos(epsilon))
  g <- v/(lambda^2 * sector)
  
  res <- list(est=g, r=r, grid=dirs$grid, theta=dirs$theta, unit=dirs$unit, epsilon=epsilon, 
              theo=rep(1, length(r)), correction=correction, h=h, delta=v, dim=dim)
  
  class(res) <- c("pcf_sector", "sector", is(res))
  res
}
