#' Directed geometric graph
#'
#' Compute the geometric graph \deqn{x_i ~ x_j \Leftrightarrow ||x_i-x_j|| < r \& -\theta < angle(x_j-x_i) < \theta}{x1 ~ x2 <=> ||x1-x2||< r & -theta < angle(x2-x1, u) < theta} 
#'
#' @param x matrix of location coordinates
#' @param u unit vector of direction
#' @param theta angle for the directed cone
#' @param double use the double cone, i.e. look also in direction -u?
#' @param r radius for neighbourhood
#' @param from default 1:nrow(x), indices of points from which to compute edges
#' @param to default 1:nrow(x), indices of points which account as potential neighbours
#' @param pregraph A list-of-neighbour indices where to look for potential neighbours
#' 
#' @details
#' If pregraph is given, "to" is ignored. 
#' 
#' @return 
#' Neighbour index vectors as a list.
#' @export

directed_geom <- function(x, u, theta, r, from, to, double=TRUE, pregraph) {
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  if(missing(r)) stop("r needed")
  if(missing(u)) stop("u needed")
  if(missing(theta)) stop("theta needed")
  if(!double) stop("Not implemented.")
  # make sure unit vector
  u <- u/sqrt(t(u)%*%u)
  if(missing(pregraph)) c_directed_geom(x, u, theta, r, from, to) 
  else c_directed_geom_by_cut(x, u, pregraph, theta, r, from, to) 
}
