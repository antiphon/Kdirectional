#' geometric graph
#'
#' @param x matrix of location coordinates
#' @param from default 1:nrow(x), indices for which to compute edges
#' @param to default 1:nrow(x), indices which account as potent. neighbours
#' @param r radius for neighbourhood
#' @export

geom <- function(x, from, to, r) {
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  if(missing(r)) stop("r needed")
  c_geom(x, from, to, r) 
}

