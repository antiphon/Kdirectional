#' Translation correction weights for point pattern
#' 
#' Defined in a bounding box.
#' 
#' @param x Pattern, $x coordinate matrix, $bbox bounding box.
#' @param asMatrix Return a square symmetric matrix? Otherwise only the upper triangle row by row.
#' 
#' @return if asMatrix=FALSE, the vector gives the weights in order 
#' (1,2), (1,3), ..., (1,n), (2,3), (2,4), ... (n-2, n-1), (n-2,n), (n-1, n)
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