#' Bartlett's spectral density
#' 
#' Periodogram estimator
#' 
#' @param x point pattern, list of $x and $bbox
#' @param k The frequencies in integers from -n to n. See details.
#' @param std scale to unit window?
#' @param ... ignored
#' 
#' @details We estimate the spectral density using a periodogram as suggested by Bartlett 1964.
#' The frequencies that the spectrum is estimated are omega = 2*pi*k/n, 
#' where k is a d-dimensional signed integer, (=integer frequency).
#' The k as input therefore should be either 1) a vector such as -16:16 (default) 
#' which will be expanded to integer frequencies by taking all pairs using expand.grid or 
#' 2) column matrix of dimension m x d, each row giving one integer frequency. 
#' 
#' No edge correction is applied. 
#' 
#' If std = TRUE, the pattern is first scaled to unit cube.
#' 
#' @export

bartletts_spectral_density <- function(x, k = -16:16, std = FALSE, ...) {
  x <- check_pp(x)
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box), or something that can be interpreted as such.")
  #
  bbox <- x$bbox
  loc <- as.matrix(x$x)
  dim <- ncol(bbox)
  n <- nrow(loc)
  #
  if(missing(k)) {
    k <- -16:16
  }
  if(is(k,"vector")){
    kv <- k
    k <- as.matrix( expand.grid(k, k) )
    if(dim == 3) k <- as.matrix( expand.grid(kv, k) )
  }
  else{
    if(ncol(k) != dim) stop("k should be a vector or a m x d matrix, d=dimension, m=n. of frequencies to compute.")
  }
  #
  #
  sl <- bbox_sideLengths(bbox)
  omega <- 2 * pi * t(t(k)/sl)
  #
  # shift to origin
  loc <- t(t(loc) - bbox[1,])
  bbox <- t(t(bbox) - bbox[1,])
  #
  # Scaling voodoo
  if(std){
    loc <- n * t(t(loc) / sl) # this weird scaling
    bbox <- t(t(bbox) / sl) #* n^dim
    omega <- t(t(omega)/sl)
  }
  V <- bbox_volume(bbox)
  #
  # Compute:
  J <- colSums(exp( -1i * loc %*% t(omega)))
  sdf <- (Re(J)^2 + Im(J)^2)  / V
  #
  #
  # Drop (0,0,..) value
  zero <- which( apply(k==0, 1, all) )
  zerov <- sdf[zero]
  sdf[ zero ] <- NA
  #zerov <- NULL
  
  # correct for bias?
  #sdf <- sdf - n^2/V
  #
  stops <- lapply(1:dim, function(d) unique(omega[,d]))
  #
  list(sdf_estimate=sdf, sdf_matrix = matrix(sdf, nrow = length(stops[[1]]) ), 
       stops = stops, 
       bbox=bbox, zerov=zerov)
}




