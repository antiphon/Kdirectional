#' Directed F, empty space function
#' 
#' We use toroidal
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector of direction
#' @param theta angle for the directed cone
#' @param r radius vector at which to evaluate K
#' @param n number of uniform dummies
#'
#' @useDynLib Kdirectional
#' @export

Fest_directional <- function(x, u, theta, r, n) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  N <- nrow(x$x)
  sidelengths <- apply(bbox, 2, diff)
  if(missing(u)) stop("u needed")
  if(missing(theta)) stop("theta needed")
  if(abs(theta)>pi/2) stop("Theta should be in range [0, pi/2]")
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  #'
  nDummies <- if(missing(n)) 4*N else n
  rcen <- max(r)
  dummies <- sapply(1:dim, function(i) runif(nDummies, bbox[1,i]+r, bbox[2,i]-r))
  #'
  
  # make sure unit vector
  u <- u/sqrt(t(u)%*%u)
  # start:
  xc <- as.matrix( rbind(as.matrix(x$x), dummies) )
  from <- 1:nDummies + N
  to <- 1:nrow(x$x)
  Nr <- length(r)
  
  #' the F value
  Fv <- function(edges) {
    # edge correction
    sum( sapply(edges[from], length ) > 0)/nDummies
  }
  #' the largest graph:
  edges <- as.list( c_directed_geom(xc, u, theta, r[Nr], from, to) )
  #'
  res <- c(Fv(edges), r[Nr])
  for(i in (Nr-1):1) {
    edges <- c_cutgeom(xc, edges, r[i]) 
    res <- rbind(c(Fv(edges), r[i]), res)
  }
  #'
  lambda_ang <- N/prod(apply(bbox, 2, diff)) * 2*theta/pi
  bd <- if(dim==3) 4/3 * pi else pi
  #'
  warning("Experimental.")
  list(rs=res[,1], theo=1-exp(-lambda_ang*bd*r^dim), r=r)
}


#' Compute all directed F functions along axis
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param ... parameters except u for Kest_directional
#' @export

Fest_along_axis <- function(x, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(x$x)
  vals <- NULL
  vecs <- rbind(c(1,0,0), c(0,1,0), c(0,0,1))
  if(dim==2) vecs <- vecs[1:2, 1:2]
  # compute Kdir in different directions along axis
  for(i in 1:dim){
    est <- Fest_directional(x, u=vecs[i,], ...)
    vals <- rbind(vals, est$rs)
  }
  r <- est$r
  list(r=r, vals=vals, dim=dim, theo=est$theo)
}




