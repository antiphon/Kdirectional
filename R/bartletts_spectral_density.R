#' Bartlett's spectral density
#' 
#' Periodogram estimator
#' 
#' @param x point pattern, list of $x and $bbox
#' @param omega The frequencies
#' @param ... ignored
#' 
#' @details We estimate the spectral density using a periodogram as suggested by Bartlett 1964,
#' 
#' \deqn{\mathcal{F}(\omega) = \lambda + \lambda^2\int_{R^d}[g(z)-1]e^{-i\omega^T z}dz}
#' 
#' where we assume that the process is stationary with intensity $\lambda$ and pair correlation function $g$. 
#' Isotropy is not assumed.
#' 
#' This function deliberately does not scale or assume a form for the frequencies (such as 2kpi/n, k integer), as 
#' there is no consensus on the best form. 
#' 
#' Additionally, no scaling or other transformation of the pattern is conducted. 
#' 
#' The frequencies that the spectrum is estimated are given by omega: 
#' 
#' 1) If omega is a vector, we expand it to d-dimensional frequencies using expand.grid.
#' 
#' 2) If omega is a column matrix of dimension m x d, each row is interpreted as a frequency. 
#' 
#' No edge correction is applied as none is known. The 
#' 
#' @export

bartletts_spectral_density <- function(x, omega, ...) {
  x <- check_pp(x)
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box), or something that can be interpreted as such.")
  #
  bbox <- x$bbox
  loc <- as.matrix(x$x)
  dim <- ncol(bbox)
  n <- nrow(loc)
  #
  if(missing(omega)) {
    stop("omega is missing.")
  }
  o <- omega
  if(is(o, "vector")){
    ov <- o
    o <- as.matrix( expand.grid(ov, ov) )
    if(dim == 3) o <- as.matrix( expand.grid(o, ov) )
  }
  else{
    if(ncol(o) != dim) stop("omega should be a vector or a m x d matrix, d=dimension, m=number of frequencies to compute.")
  }
  #
  V <- bbox_volume(bbox)
  #
  # Compute:
  J <- colSums(exp( -1i * loc %*% t(o)))
  sdf <- (Re(J)^2 + Im(J)^2)  / V
  #
  #
  # Drop (0,0,..) value
  zero <- which( apply(o==0, 1, all) )
  zerov <- sdf[zero]
  sdf[ zero ] <- NA
  #
  stops <- lapply(1:dim, function(d) unique(o[,d]))
  #
  list(sdf_estimate=sdf, 
       sdf_matrix = matrix(sdf, nrow = length(stops[[1]]) ), 
       omega = o,
       stops = stops, 
       bbox=bbox, 
       zerov=zerov)
}




