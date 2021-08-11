#' nearest neighbour angles and distances
#' 
#' @param x Matrix of coordinates, n rows and d columns
#' 
#' @details 
#' 
#' If 2d, returns the angle in [-pi,pi]. In 3d also the angle from xy-plane in
#' [0,pi]
#' 
#' @return 
#' The values in 2d are (azimuth, distance).
#' The values in 3d are physical coordinate system (azimuth, inclination, distance)
#' 
#' @examples 
#' x <- matrix(runif(200), ncol = 2)
#' nn <- nnangle(x)
#' summary(nn)
#' 
#' @useDynLib Kdirectional 
#' @export
nnangle <- function(x, from, to){
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  d<-c_angles(x, from, to)
  v<-do.call(cbind, d)
  colnames(v) <- c("ang", if(ncol(v)==2) "dist" else c("ang2", "dist"))
  class(v) <- c("nnangle", is(v))
  
  v
}


#' Summary method for nearest neighbour angles
#' 
#' @param x object from nnangle-function
#' @param ... passed on to ks.test for testing uniformity of the angles
#' 
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
