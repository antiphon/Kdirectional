#' Construct a bbpoly from column range matrix
#' 
#' @param bb bbox
#' 
#' @export
bbox2bbpoly  <- function(bb){
  if(ncol(bb)==3) bbpoly_default(bb[,1],bb[,2],bb[,3])
  else bbpoly_default(bb[,1], bb[,2], NA)
}

#' Return the bounding box of a bbpoly
#' 
#' @param x bbpoly
#' 
#' @export
bbpoly2bbox <- function(x) {
  apply(x$vertices, 2, range)
}

#' Simplex object
#' 
#' create a simplex (triangle)
#' 
#' @param d dimension, 2 or 3
#' 
#' @export
bbpoly_simplex <- function(d=3){
  if(d==3) {
    vertices <- cbind(x=c(0,.5, 1, .5), y=c(0, 1, 0, 0.5), z=c(0,0,0,1))
    e0 <- cbind(c(1,1,1,2,2,3),
              c(2,3,4,3,4,4))
    f0 <- list(vertices = list(c(1,2,4), c(2,3,4), c(3,1,4), c(1,2,3)),
             edges = c(1,5,3), c(4,6,5), c(2,3,6),c(1,4,2))
  
  as_bbpoly(list(vertices = vertices,
                 edges = e0,
                  faces = f0))
  }
  else{
    as_bbpoly(list(vertices = cbind(x=c(0,0,0.86), y=c(0,1,.5)),
                   edges = cbind(c(1,1,2),c(2,3,2)) ) )
  }
}


#' Affine transformation of a polytope
#' 
#' @param bbpoly bbpoly object
#' @param A linear transformation
#' @param s shift vector
#' @param center_before_A Center the polygon to mass centrum before applying A? Reversed after A and before shift.
#' 
#' 
#' 
#' @export
bbpoly_affine <- function(bbpoly, A, s=c(0,0,0), center_before_A=FALSE){
  if(is.bbpoly(bbpoly)){
    v <- bbpoly$vertices
    d <- ncol(v)
    c0 <- (if(center_before_A) c0 <- apply(v, 2, mean)  else c(0,0,0))[1:d]
    v <- t( t(v) - c0 )
    v1 <- coord_affine(v, A)
    v1 <- t( t(v1) + c0)
    bbpoly$vertices <- coord_affine(v1, diag(1,d), s[1:d])
    colnames(bbpoly$vertices) <- c("x","y","z")[1:d]
    bbpoly$volume <- bbpoly$volume * abs( det(A) )
  }
  else{
    if(is.matrix(bbpoly)){
      d <- ncol(bbpoly)
      bbpoly <- t(t(bbpoly %*% t(A)) + s[1:d])
    }
    else stop("Dont recognise input.")
  }
  bbpoly
}


#' bbpoly from 2d polygon given by x and y coordinates
#' 
#' @param x x coordinates
#' @param y y coordinates
#' 
#' @details 
#' Assuming an open polygon, meaning the last vertex is not the 
#' same as the first. Will check and adapt.
#' 
#' @export
poly2bbpoly <- function(x, y) {
  n <- length(x)
  if(!(x[n]!=x[1] | y[n]!=y[1])) { # closed, unclose.
    x <- x[-n]
    y <- y[-n]
  }
  n <- length(x)
  if(n<3) stop("at least 3 unique points required.")
  bbp <- list(vertices=cbind(x,y), 
              edges = cbind(1:n, c((1:(n-1)+1), 1) ))
  # store the volume
  px <- bbp$vertices
  n <- nrow(px)
  px <- bbp$vertices
  # Shoelace method
  bbp$volume <- 0.5 * abs( sum(px[-n,1]*px[-1,2]) - sum(px[-1,1]*px[-n,2]) )
  class(bbp) <- "bbpoly"
  bbp
}

#' Turn something into a bbpoloy
#'
#' @param x input
#' 
#' @details understands: matrix (bbox); bbquad; list with vertices, egdes and faces in 3d
#' 
#' @export
as_bbpoly <- function(x) {
  if(is.bbquad(x)) stop("Not yet implemented")
  if(is.matrix(x)) bbox2bbpoly(x)
  else if(is.list(x)) {
    if(is.null(x$vertices)) stop("vertices element missing in the input")
    d <- ncol(x$vertices)
    if(is.null(x$edges)) stop("Edges element missing in the input")
    if(d==3)if(is.null(x$faces)) stop("Faces element missing in the input")
    class(x) <- c("bbpoly", is(x))
    x
  }
  else stop("Don't know what to do with the input.")
}


#' Regular bounding box as a polytope
#' 
#' @param xl x range Default -.5,5
#' @param yl y range. Default = xl
#' @param zl z range. Default = xl. Set it to NA for 2D.
#' 
#' @details 
#' Default output is a 3D cube with vertices in the corners of [-.5,.5]^3. 
#' 
#' See also \code{\link{poly2bbpoly}} for more flexible windows in 2D.
#' 
#' @export
bbpoly_default <- function(xl=c(-.5,.5), yl=xl, zl=xl){
  # 2D
  if(any(is.na(zl))){
    vert <- as.matrix(expand.grid(x=xl, y=yl))
    # quadrilateral polytope, or rectangle
    bbp   <- list(vertices=vert, edges=cbind(c(1,2,4,3), 
                                             c(2,4,3,1)))
    # store the volume
    bbp$volume <- abs( diff(xl)*diff(yl) )
    #    
  }
  else{
    verts <- as.matrix(expand.grid(x=xl, y=yl, z=zl))
    edges <- cbind(c(1,1,3,4,5,5,7,8,1,3,4,2),
                   c(2,3,4,2,6,7,8,6,5,7,8,6))
    # quads, 6 of them. List as we can have,
    # eg triangles and rectangles mixed
    faces_v <- list(c(1,3,4,2),
                    c(1,5,7,3),
                    c(3,7,8,4),
                    c(2,4,8,6),
                    c(1,2,6,5),
                    c(5,6,8,7))
    # edge indices of each face 
    faces_e <- list(c(1,2,3,4),
                    c(2,9,6,10), 
                    c(10,7,11,3), 
                    c(4,12,8,11),
                    c(1,9,5,12),
                    c(5,6,7,8))
    # reverse the orders to get invards pointing normals
    faces_v <- lapply(faces_v, rev)
    bbp <- list(vertices = verts, 
                edges = edges,
                faces = list(vertices = faces_v, 
                                edges = faces_e))
    # store the volume
    bbp$volume <-abs( diff(xl)*diff(yl)*diff(zl) )
    #
  }
  class(bbp) <- "bbpoly"
  bbp
}

#' check class bbpoly
#' @export
is.bbpoly <- function(x, ...){
  is(x, "bbpoly")
}

#' Check with error
#' 
#' @param x object to check for bbpoly class
#' 
#' @export
check_bbpoly <- function(x) if(!is.bbpoly(x)) stop("x not a bbpoly-object.")


#' Surface planes of bbpoly
#'
#' @param x bbpoly object
#' 
#' @details Returns the normals of the edges (2d) or faces (3d). Assumes  that all vertices of each face are on the same plane. 
#' 
#'@export
bbpoly_planes <- function(x) {
  d <- ncol(x$vertices)
  if(d == 3) {
    sapply(x$faces$vertices, function(id) {
      v <- x$vertices[id[1:3],] # Here we assume co-planar face
      n <- cross(v[2,]-v[1,], v[3,]-v[1,])
      n <- n/sqrt(sum(n^2))
      names(v) <- c("x", "y", "z")
      # location point, take geometric mean
      p <- colMeans(x$vertices[id,])
      c(normal=n, p=p)
    } )
  }
  else{  # 2D
    M <- cbind(c(0,1),c(-1,0))
    idx <- unique( c(t(x$edges)) )
    idx <- c(idx, idx[1])
    v <- x$vertices
    sapply(1:(length(idx)-1), function(i){
      u <- v[idx[i+1],]-v[idx[i],]
      n <- M%*%u
      n <- n/sqrt(sum(n^2))
      p <- (v[idx[i+1],]+v[idx[i],])/2
      c(normal=n, p=p)
    })
  }
} 


#' Plot bbpoly object
#' 
#' @param x bbpoly object
#' @param normals draw normals
#' @param faces draw faces
#' @param ecol edge color
#' @param ncol normal color
#' @param add add or not
#' @param ... passsed to shade3d or lines (2d)
#' 
#' @export
#' @import rgl

plot.bbpoly <- function(x, normals=FALSE, faces = FALSE,
                        ecol=2, ncol=4, add=TRUE, ...){
  d <- ncol(x$vertices)
  if(d==3){
    plot3d(x$vertices, aspect = FALSE, add=add, ..., col=ecol)
    if(faces) {
      null <- lapply(x$faces$vertices, function(f)
      shade3d(qmesh3d(t(x$vert), f, homogeneous = FALSE ), ... ))
    }
    edgs <- cbind(x$vertices[x$edges[,1],],x$vertices[x$edges[,2],])
    null<-apply(edgs, 1, function(z)
                lines3d(rbind(z[1:3],z[4:6]), col=ecol)  )
  }
  else{ # 2D
    if(!add) plot(x$vertices, pch=NA, asp=TRUE)
    lines(x$vertices[t(x$edges),], col=ecol, ...)
  }
  #
  if(normals){
    pl <- bbpoly_planes(x)
    n <- ncol(pl)
    bi <- 1:d
    ai <- bi+d
    linesf <- list(lines, if(exists("lines3d")) lines3d else NULL) [[d-1]]
    for(i in 1:n) linesf((rbind(pl[ai,i], pl[bi,i]+pl[ai,i])), col=ncol)
  }
}


#' print method for bbpoly, a cuboid made of quads
#'
#' @param x object
#' @param ... ignored
#'
#' @export

print.bbpoly <- function(x, ...){
  d <- ncol(x$vertices)
  cat(d, "D polytopy, ", sep="")
  cat(nrow(x$vertices), "vertices")
  if(!is.null(x$edges))
    cat(",", nrow(x$edges), "edges")
  if(!is.null(x$faces)) 
    cat(",", length(x$faces$vertices), "faces")
  
  cat("\nVolume:", x$volume)
  cat("\n")
}
