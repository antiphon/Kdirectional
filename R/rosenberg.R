#' Rosenberg Wavelet Transform
#' 
#' Estimate the wavelet variance suggested by Rosenberg 2004. 2D only.
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' 
#' @export
rosenberg <- function(x) {
  x <- check_pp(x)
  if(ncol(x$x) != 2) stop("Use only for 2D data")
  #
  # sectors:
  
}
