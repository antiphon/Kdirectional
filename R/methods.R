#' method plot for O-S K-function
#' 
#' @export
#' @exportMethod plot

plot.os_K <- function(x, ...){
  image(z=x$Z, x=x$r, y=x$theta[[1]], xlab="range", ylab="angle", ...)
}

#' method plot for anisotropic pcf-function
#' 
#' @import sphere
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
#' @import sphere
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

#' Print for nearest neighbour angles
#' 
#' @exportMethod print
#' @export
print.nnangle <- function(x, ...) {
  cat("Nearest neighbour angles and distances\n")
  #'
  k <- attr(x, "k")
  if(!is.null(k)) cat(paste0(k,"th nearest neighbours\n"))
  # check if in cone:
  dim <- ncol(x)
  unit <- attr(x, "unit")
  if(!is.null(unit)){
    type <- ifelse(dim==2, "sector", "cone")
    cat("inside a directed ", type, ":\n", sep="")
    cat(" * central direction unit-vector: (", paste0(unit, collapse=","), ")\n" , sep="")
    cat(" * Angle (width) of the ", type, ": ", attr(x, "theta"), sep="")
  }
  cat("\n\n Data looks like:\n")
  print(head(x))
  if(nrow(x)>6)cat("     ...\n")
  cat("Number of points:", nrow(x))
  if(dim==2){
    cat("\nFirst column: Angle [0, 2pi]")
    cat("\nSecond column: Distance")
  }else{
    cat("\nFirst column: Azimuth [0, 2pi], observed")
    cat("\nSecond column: Inclination (polar angle) [0, pi]")
    cat("\nThird column: Distance")
  }
  
}

#' Summary method for nearest neighbour angles
#' 
#' @exportMethod summary
#' @export

summary.nnangle <- function(x, ...) {
  dim <- ncol(x)
  ranges <- round(apply(x, 2, range), 5)
  if(dim==2) {
    cat("Angle [0, 2pi], observed [", paste0(ranges[,1], collapse=","),"]", sep="")
    cat("\nDistance, observed [", paste0(ranges[,2], collapse=","),"]", sep="")
    cat("\n\n")
    # Kolmogorov-Smirnov test
    p0 <- function(x) punif(x, 0, 2*pi)
    angle <- x[,1]
    K <- ks.test(angle, p0, ...)
    cat("Kolmogorov-Smirnov test against uniformity:\n")
    cat("D=", K$statistic, " p=", K$p.value, "\n")
  }
  else{
    cat("\nAzimuth [0, 2pi], observed [", paste0(ranges[,1], collapse=","),"]", sep="")
    cat("\nInclination (polar angle) [0, pi], observed [", paste0(ranges[,2], collapse=","),"]", sep="")
    cat("\nDistance, observed [", paste0(ranges[,3], collapse=","),"]", sep="")
    # Kolmogorov-Smirnov tests
    p0 <- function(x) punif(x, 0, 2*pi)
    angle1 <- x[,1]
    K1 <- ks.test(angle1, p0, ...)
    p1 <- function(x) 0.5 * (1-cos(x))
    angle2 <- x[,2]
    K2 <- ks.test(angle2, p1, ...)
    cat("\n\nKolmogorov-Smirnov test against isotropy:\n")
    cat("Azimuth: D=", K1$statistic, " p=", K1$p.value, "\n")
    cat("Inclination: D=", K2$statistic, " p=", K2$p.value, "\n")
    K <- list(azimuth=K1, inclination=K2, pazimuth=p0, pinclination=p1)
  }
  invisible(K)
}

#' Plot for nnangle
#' 
#' @exportMethod plot
#' @export
plot.nnangle <- function(x, ...) {
  dim <- ncol(x)
  ang <- x[,1]
  lab <- ifelse(dim==2, "Angle", "Azimuth")
  hist(ang, xlab=lab, main="", freq=FALSE, ...)
  abline(h=1/(2*pi), col=3)
  if(dim==3){
    ang2 <- x[,2]
    h <- hist(ang2, main="", xlab="Inclination (polar angle)", freq=FALSE, ...)
    p1 <- function(x) 0.5*sin(x)
    curve(p1, add=T, col=3, from=0, to=pi)
  }
}



