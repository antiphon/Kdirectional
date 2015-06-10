#' check pattern format
#' @param x pattern candidate
#' 
#' @export

check_pp <- function(x){
  if("ppp"%in%is(x)){
    x <- list(x=coords(x))
    x$bbox <- apply(x$x, 2, range)
  }
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  x$x <- as.matrix(x$x)
  if(is.bbquad(x$bbox)){ # ok, we have 3d box 
  }
  else{
    x$bbox <- as.matrix(x$bbox)
    if(ncol(x$x)!=ncol(x$bbox)) stop("x$x and x$bbox have different dimensions.")
  }
  x
}