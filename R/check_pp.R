#' check pattern format
#' @param x pattern candidate
#' 
#' @export

check_pp <- function(x){
  if("ppp"%in%is(x)){
    window <- x$window
    x <- list(x=as.matrix( cbind(x$x, x$y) ))
    x$bbox <- owin_to_bbox(window)
  }
  if(!is(x,"list")){
    if(is(x, "matrix")) x <- list(x=x, bbox = bbox_make(x))
    else stop(paste("Can not parse x of type", is(x) )  )
  }
  else if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  x$x <- as.matrix(x$x)
  if(is.bbquad(x$bbox)){ #hmm
  }
  else{
    x$bbox <- as.matrix(x$bbox)
    if(ncol(x$x)!=ncol(x$bbox)) stop("x$x and x$bbox have different dimensions.")
  }
  x
}
