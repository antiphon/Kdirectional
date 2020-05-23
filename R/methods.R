#' method plot for O-S K-function
#' 
#' @export
#' @exportMethod plot

plot.os_K <- function(x, ...){
  image(z=x$Z, x=x$r, y=x$theta[[1]], xlab="range", ylab="angle", ...)
}

#' method plot for anisotropic pcf-function
#' 
#' @export
#' @exportMethod plot

plot.pcf_anisotropic <- function(x, ab, ...){
  s <- NULL
  # 2d:
  if(x$dim==2) image(z=x$est, x=x$r, y=x$theta[[1]], xlab="range", ylab="angle", ...)
  else{
    stop("Plot not implemented for 3D.")
    # more tricky. We will plot a integral over the range:
#     lon <- azi2lon(x$theta[[1]])
#     lat <- inc2lat(x$theta[[2]])
#     
#     if(missing(ab)) ab <- range(x$r)
#     ok <- x$r <= ab[2] & x$r >= ab[1]
#     v <- apply(x$est[ok,,], 3, apply, 2, sum) - 1
#     ll <- ai2ll(expand.grid(lat, lon))
#     #' add antipode
#     ll <- rbind(ll, antipode(ll))
#     s <- rgltools.smoothsphere(ll, c(c(v), c(v)))
#     rgltools.plotsmoothsphere(s, ...)
    
  }
  invisible(s)
}

#' method plot2d for anisotropic pcf-function
#' 
#' @export
#' @exportMethod plot2d

plot2d.pcf_anisotropic <- function(x, ab, ...){
  s <- NULL
  # 2d:
  if(x$dim==2) image(z=x$est, x=x$r, y=x$theta[[1]], xlab="range", ylab="angle", ...)
  else{
    # more tricky. We will plot a integral over the range:
    if(missing(ab)) ab <- range(x$r)
    ok <- x$r <= ab[2] & x$r >= ab[1]
    s <- apply(x$est[ok,,], 3, apply, 2, sum) * diff(x$r[1:2])
    # then we plot this using image
    image(z=s , x=x$theta[[1]], y=x$theta[[2]], xlab="azimuth", ylab="inclination (polar angle)", ...)
  }
  invisible(s)
}
