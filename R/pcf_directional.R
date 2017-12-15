#' Directed pcf function
#' 
#' We use translation edge correction if we have an axis-oriented bounding box. Otherwise minus-correction.
#'
#' @param x pp, spatstat's ppp-object, or a coordinate matrix, or a list with $x~coordinates $bbox~bounding box
#' @param u unit direction(s), if many then one direction per row.
#' @param epsilon Central angle for the directed cone (total angle is 2 epsilon)
#' @param r radius vector at which to evaluate K
#' @param ... Passed on to \code{\link{pcf_anin}}
#' @param cylindrical If TRUE, compute the cylindrical version using \code{\link{pcf_anin_cylinder}}
#' 
#' @details 
#' 
#' Compute the sector/cone/cylindrical pcf function. This version uses the more general anisotropic-inhomonogeneous \link{pcf_anin_conical} and \code{\link{pcf_anin_cylinder}}, by setting the intensity = constant. See there for further parameters, especially kernel smoothing options.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @examples
#' 
#' x <- matrix(runif(300), ncol=2)
#' k <- pcf_directional(x)
#' plot(k, rmax = 0.1)
#' @import Matrix
#' @useDynLib Kdirectional
#' @export

pcf_directional <- function(x, u, epsilon, r, ..., cylindrical = FALSE) {
  x <- check_pp(x)
  bbox <- x$bbox
  n <- nrow(x$x)
  lambda1 <- n/bbox_volume(bbox)
  lambda <- rep(lambda1, n)
  gest_f <- if(cylindrical) pcf_anin_cylinder else pcf_anin_conical
  gest <- gest_f(x, u, epsilon, r, lambda, ...)
  if(cylindrical) attr(gest, "fun_name") <- "Cylindrical pcf"
  else attr(gest, "fun_name") <- "Conical pcf"
  gest
}

