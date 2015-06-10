#' Translation correction weights for point pattern
#' 
#' Defined in a bounding box.
#' 
#' @param x Pattern, $x coordinate matrix, $bbox bounding box.
#' 
#' @export

translation_weights <- function(x, asMatrix=FALSE){
  x <- check_pp(x)
  d <- ncol(x$x)
  n <- nrow(x$x)
  
  w <- c_translation_weights(x$x, x$bbox)
  
  if(asMatrix){
    B <- diag(0, n)
    B[lower.tri(B)] <- w
    B <- B+t(B)
    diag(B) <- prod( apply(x$bbox, 2, diff)  )
    w <- B
  }
  
  w
}