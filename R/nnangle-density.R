#' Nearest neighbour angle density
#'
#' 2D only.
#' 
#' @param x point pattern
#' @param bw Epanechnikov kernel bandwidth
#' @param antipodal Work with angles only on the 0 - pi range? If false, work on 0 - 2*pi.
#' @param avec vector of angles to estimate the density at
#' @param ... ignored
#' @param data alternative to x (expert use only)
#' @param justData Just return the components for calculation?
#' 
#' 
#' @export
nn_orientation_density <- 
  nnangle.density <- function(x, bw = 0.35, antipodal = TRUE, avec, ..., data, justData = FALSE) {
  
  M <- ifelse(antipodal, 1, 2)*pi # maximum
  
  # compute nesessary components
  if(!missing(data)) {
    ang <- data$ang
    V <- data$V
    l <- data$lambda
  }
  else{
    # assume rectagle
    x <- check_pp(x)
    loc <- x$x
    bbox <- x$bbox
    n <- nrow(loc)
    # nn's
    nd <- nnangle(loc)
    ang <- nd[,1]
    nd <- nd[,2]
    if(antipodal) ang[ ang > pi] <- ang[ ang > pi] - pi
    # window erosions
    sl <- bbox_sideLengths(bbox)
    V <- sapply(nd, function(d) prod(sl-2*d) )
    # edge corrections
    ed <- bbox_distance(loc, bbox)
    #
    ok <- nd < ed
    ang <- ang[ok]
    V <- V[ok]
    # adapted intensity estimate
    l <- sum(1 / V)
    #
    #
    if(justData) return(list(ang=ang, V=V, lambda=l, bw=bw))
  }
  #
  # smoothing kernel
  kern <- function(d) {
    v <- 0.75 * (1-d^2)
    v[abs(d)>1] <- 0
    v
  }
  # estimation grid
  if(missing(avec)) avec <- seq(0, M, l = 50)
  
  # estimate
  val <- sapply(avec, function(av) {
    b <- abs(ang - av) 
    b <- pmin(b, M - b) # toroidal distance [0,M]
    v <- kern(b/bw)/bw
    sum(v/V)
  }  )
  
  est <- val/l
  data.frame(angle = avec, density = est)
}

