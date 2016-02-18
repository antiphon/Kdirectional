#' Directed K function
#' 
#' We use translation edge correction.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector of direction
#' @param epsilon Central angle for the directed cone (total angle is 2 epsilon)
#' @param r radius vector at which to evaluate K
#' @param pregraph A super graph for determining neighbourhoods
#' 
#' @details 
#' 
#' Compute the sector K function. Try \code{\link{Kest_anin}} for a bit better version that does also second order inhomogeneous patterns.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @import Matrix
#' @useDynLib Kdirectional
#' @export

Kest_directional <- function(x, u, epsilon, r, pregraph) {
  x <- check_pp(x)
  bbox <- x$bbox
  if(is.bbquad(bbox)) stop("bbquad window not yet supported.")
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  if(missing(u)) stop("u needed")
  if(missing(epsilon)) stop("epsilon needed")
  if(abs(epsilon)>pi/2) stop("epsilon should be in range [0, pi/2]")
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  # make sure unit vector
  u <- u/sqrt(t(u)%*%u)
  # start:
  xc <- as.matrix(x$x)
  from <- 1:nrow(xc)
  to <- 1:nrow(xc)
  Nr <- length(r)
  N <- nrow(xc)
  #' lambda:
  lambda <- N/prod(apply(bbox, 2, diff))
  #' the largest graph:
  edges <- directed_geom(xc, u, epsilon, r[Nr], from, to, TRUE, pregraph) 
  #'
  #' compute weights for translation correction
  W <- translation_weights(x, asMatrix=TRUE)
  #' the K value
  Kv <- function(edges) {
    has_neighs <- which(sapply(edges, length)>0)
    S <- if(length(has_neighs))
            sum(sapply(has_neighs, 
                 function(i) sum( 1/W[i, edges[[i]] ] )
                )
      ) else 0
    S/lambda^2
  }
  #'
  res <- c(Kv(edges), r[Nr])
  for(i in (Nr-1):1) {
    edges <- c_cutgeom(xc, edges, r[i]) 
    res <- rbind(c(Kv(edges), r[i]), res)
  }
  #
  bd <- pi * if(dim==3) 4/3  else 1
  #'
  list(trans=res[,1], r=r, theo=bd*r^dim)
}


#' Compute all directed K functions along axis
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param ... parameters except u for Kest_directional
#' @export

Kest_along_axis <- function(x, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(x$x)
  Kvals <- NULL
  vecs <- rbind(c(1,0,0), c(0,1,0), c(0,0,1))
  vecs <- vecs[1:dim, 1:dim]
  # compute Kdir in different directions along axis
  for(i in 1:dim){
    kest <- Kest_directional(x, u=vecs[i,], ...)
    Kvals <- rbind(Kvals, kest$trans)
  }
  r <- kest$r
  list(r=r, vals=Kvals, dim=dim, theo=kest$theo)
}




