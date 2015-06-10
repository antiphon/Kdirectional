#' Ohser-Stoyan K-function
#' 
#' @param x point pattern, $x and $bbox
#' @param r vector of ranges
#' @param theta vector of angles. list with two components if 3D
#' 
#' @export

osK <- function(x, r, theta, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  if(missing(theta)){
    theta <- list(angle=seq(0, pi, length=20))
    #' add also pi if not already
    theta[[1]] <- union(theta[[1]], pi)
    if(dim==3) theta <- append(theta, list(angle2=seq(0, pi, length=20)))
  }
  
  if(r[1]!=0) stop("r[1] needs to be 0.")
  
  if(!is.list(theta)) 
    stop("Theta should be a list giving 
          increasing azimuth (and in 3d inclination) vector(s), in physical spherical coordinate system.")    
  xc <- x$x
  v <- c_oh_K(xc, theta, r, bbox)
  #'
  lambda <- nrow(xc)/prod(sidelengths)
  K <- v/lambda^2
  #'
  grid <- expand.grid(append(list(r=r), theta))
  #' compile useful data
  ntheta <- sapply(theta, length)
  if(dim==2) {
    z <- matrix(K, nrow=length(r))
    #' isotropic K
    K_iso <- data.frame(r=r, trans=2*z[,ntheta])
  }
  else{
    dd <- c(length(r), ntheta); names(dd) <- c("range", "azimuth", "inclination (polar angle)")
    z <- array(K, dim=dd)
    K_iso <- data.frame(r=r, trans=2*z[ , ntheta[1], ntheta[2]])
  }
  #'
  #'
  res <- list(Z=z, K_iso=K_iso, v=K, total_orientation=tot_ori, theta=theta, r=r, grid=grid)
  class(res) <- "os_K"
  res
}