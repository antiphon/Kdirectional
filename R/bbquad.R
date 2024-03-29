#' Construct a bbquad from column range matrix
#' 
#' @export
bbox2bbquad  <- function(bb){
  if(ncol(bb)==3) bbquad_default(bb[,1],bb[,2],bb[,3])
  else bbquad_default(bb[,1],bb[,2],NA)
}

#' Return the bounding box of a bbquad
#' 
#' @param x bbquad
#' 
#' @export
bbquad2bbox <- function(x) {
  apply(x$vertices, 2, range)
}

#' bbquad from 2d polygon
#' 
#' @param x,y polygon coordinates
#' 
#' @export
poly_to_bbquad <- function(x, y) {
  quads <- rbind(1:length(x))
  bbq<-list(vertices=cbind(x,y), idx=quads)
  # store the volume
  px <- bbq$vertices
  n <- nrow(px)
  px <- rbind(px,px[1,])
  # Shoelace method
  bbq$volume <- 0.5 * abs( sum(px[-n,1]*px[-1,2]) - sum(px[-1,1]*px[-n,2]) )
  class(bbq) <- "bbquad"
  bbq
}

#' Regular bounding box as bounding box of quadrilaterals
#' 
#' @param xl x range Default -.5,5
#' @param yl y range. Default = xl
#' @param zl z range. Default = xl. Set it to NA for 2D.
#' 
#' @details 
#' Default output is a 3D cube with vertices in the corners of [-.5,.5]^3. 
#' 
#' See also \code{\link{poly_to_bbquad}} for more flexible windows in 2D.
#' 
#' @export
bbquad_default <- function(xl=c(-.5,.5), yl=xl, zl=xl){
  # special case of 2D
  if(any(is.na(zl))){
    vert <- as.matrix(expand.grid(x=xl, y=yl))
    # quads, just 1 of them
    quads <- rbind(2,4,3,1)
    bbq<-list(vertices=vert, idx=quads)
    # store the volume
    bbq$volume <-abs( diff(xl)*diff(yl) )
    #    
  }
  else{
    vert <- as.matrix(expand.grid(x=xl, y=yl, z=zl))
    # quads, 6 of them
    quads <- cbind(c(1,3,4,2),
                   c(1,5,7,3),
                   c(3,7,8,4),
                   c(2,4,8,6),
                   c(1,2,6,5),
                   c(5,6,8,7))
    quads <- quads[4:1,]
    #m<-qmesh3d(t(vert), quads, homog=F) 
    #text3d(vert, texts = 1:8)
    #plot3d(vert, size=20, col=3)
    #shade3d(m, col=rgb(1,0,0,0.1))
    bbq<-list(vertices=vert, idx=quads)
    # store the volume
    bbq$volume <-abs( diff(xl)*diff(yl)*diff(zl) )
    #
  }
  class(bbq) <- "bbquad"
  bbq
}


#' print method for bbquad, a cuboid made of quads
#' 
#' @param x bbquad
#' @param ... ignored
#' 
#' @export

print.bbquad <- function(x, ...){
  d <- ncol(x$vertices)
  cat(d, "D quadrilateral bounding box", sep="")
  cat("\n", nrow(x$vert), "vertices")
  cat("\n", ncol(x$idx), "faces")
  cat("\nVolume:", x$volume)
  cat("\n")
}

#' check class bbquad
#'
#' @param x object to check
#' @param ... ignored
#' 
#' @export
is.bbquad <- function(x, ...){
  "bbquad"%in%is(x)
}

#' Plot bbquad object
#'
#' @param x bbquad
#' @param edges draw edges
#' @param normals draw normals
#' @param faces draw faces
#' @param ecol edge color
#' @param ncol normal color
#' @param add add or not
#' @param ... passsed to shade3d or lines (2d)
#' 
#' @importFrom rgl qmesh3d lines3d shade3d
#' @export

plot.bbquad <- function(x, edges = FALSE, 
                        normals = FALSE, 
                        ecol=2, ncol=4, add=TRUE, ...){
  d <- ncol(x$vertices)
  if(d==3){
    m<-qmesh3d(t(x$vert), x$idx, homogeneous = FALSE)
    shade3d(m, ...)
    if(edges){
      edgs <- bbquad_edges(x)
      null<-apply(edgs, 2, function(z)  lines3d(rbind(z[1:3],z[4:6]), col=ecol)  )
    }
  }else{
    if(!add) plot(x$vertices, pch=NA, asp=TRUE)
    lines(x$vertices[c(x$idx,x$idx[1]),], col=ecol, ...)
  }
  if(normals){
    pl <- bbquad_planes(x)
    n <- ncol(pl)
    bi <- 1:d
    ai <- bi+d
    linesf <- list(lines, if(exists("lines3d"))lines3d else NULL) [[d-1]]
    for(i in 1:n) linesf((rbind(pl[ai,i], pl[bi,i]+pl[ai,i])), col=ncol)
  }
}

#' Surface planes
#' 
#' normal vector and a point -format. In 2D returns the surface lines.
#' 
#' @param b bbquad
#' 
#' 
#' @export
bbquad_planes <- function(b){
  d <- ncol(b$vertices)
  if(d==3){
    apply(b$idx, 2, function(id){
      x<-b$vertices[id[1:3],] 
      n <- cross(x[2,]-x[1,], x[3,]-x[1,])
      n <- n/sqrt(sum(n^2))
      # location point, take geometric mean
      p <- colMeans(b$vertices[id,])
      c(normal=n, p=p)
    })
  }
  else{ # 2d
    M <- cbind(c(0,1),c(-1,0))
    idx <- c(b$idx,b$idx[1])
    v <- b$vertices
    sapply(1:(length(idx)-1), function(i){
      u <- v[idx[i+1],]-v[idx[i],]
      n <- M%*%u
      n <- n/sqrt(sum(n^2))
      p <- (v[idx[i+1],]+v[idx[i],])/2
      c(normal=n, p=p)
    })
  }
}



#' Distance from points to bbquad walls
#' 
#' @param x n x 3 matrix of locations
#' @param b quadrilateral bounding box 
#' @param warn Emit a warning for points outside box? (default: TRUE)
#' 
#' Note: Returned distance will be negative if point outside the box.
#' 
#' @export

bbquad_distance <- function(x, b, warn = TRUE){
  # form all planes/lines
  planes <- bbquad_planes(b)
  # compute distance from each point to each plane
  dists <- apply(planes, 2, dist_point_to_plane, x=x, absit = FALSE)
  # check every point inside
  if(warn){
    neg <- which(apply(dists, 1, function(d) any(d<0)))
    nn <- length(neg)
    if(nn) warning(paste0(nn, " points ", ifelse(nn<5, paste0("(idx: ", paste0(neg, collapse=","), ") "), ""), "outside the box."))
  }
  apply(dists, 1, min)
}


#' Cross-product of two 3D vectors
#' 
#' @export
cross <- function(u, v){
  c(u[2]*v[3]-u[3]*v[2], 
    u[3]*v[1]-u[1]*v[3],
    u[1]*v[2]-u[2]*v[1])
}

#' Three points to a plane
#' 
#' @param x vector of 3 coordinates
#' 
#' @export
three_points_to_plane <- function(x){
  x <- rbind(x)
  if(nrow(x)!=3) stop("x should be 3 non-collinear points.")
  # normal
  n <- cross(x[2,]-x[1,], x[3,]-x[1,])
  n <- n/sqrt(sum(n^2))
  # location point, take geometric mean
  p <- colMeans(x)
  c(normal=n, p=p)
}

#' Point to plane distance
#' 
#' @param x vector of 3 coordinates
#' @param pl plane in norm-point format
#' @param absit take abs? Otherwise sign gives side
#' 
#' @export
dist_point_to_plane <- function(x, pl, absit = TRUE){
  x <- rbind(x)
  d <- ncol(x)
  n <- pl[1:d]
  p0 <- pl[1:d+d]
  d <- apply(x, 1, function(p) t(p-p0)%*%n )
  if(absit) abs(d) else d
}

#' Get the edges of the quad box
#' 
#' @param x bbquad
#' 
#' @export
bbquad_edges <- function(x){
  d <- ncol(x$vertices)
  ei <- if(d==3)rbind(c(1,2), c(1,3), c(1,5), 
                      c(2,4), c(2,6), c(3,4), 
                      c(3,7), c(4,8), c(5,7),
                      c(5,6), c(6,8), c(7,8)) 
         else rbind(c(1,2), c(1,3),c(2,4),c(3,4))
  apply(ei, 1, function(id) {
      z <- t(x$vertices[id,])
      c(start=z[1:d], end=z[1:d+d]) 
    })
}

#' Triangulate a bbquad 
#' 
#' Due to assuming quadrilateral polytope, we get 2 simplices (triangles/tetrahedrons)
#' 
#' @param x bbquad object
#'
#' @export

bbquad_simplices <- function(x) {
  if(!is.bbquad(x)) stop("x not a bbquad-object.")
  d <- ncol(x$vertices)
  if(d==2){
    list(x$vertices[1:3,], x$vertices[c(1,3,4),])
  }
  else if(d==3){
    if(nrow(x$vertices)!=8) stop("3D bbquad object not in quadrilateral form.")
    s <- list()
    # the central vertex
    c0 <- apply(x$vertices, 2, mean)
    # split each face to two tetrahedrons
    k <- 0
    for(i in 1:6){
      j <- x$idx[,i]
      s[[k<-k+1]] <- rbind(x$vertices[j[1:3],], c0)
      s[[k<-k+1]] <- rbind(x$vertices[j[c(1,3,4)],], c0)
    }
    s # done
  }
  else stop("bbquad not 2d or 3d.")
}

#' Volume of a bbquad
#' 
#' @param x bbquad-object
#' 
#' @details 
#' Will split the bbquad into simplices and return the sum of simplex volumes.
#' 
#' @export

bbquad_volume <- function(x) {
  s <- bbquad_simplices(x)
  sum(simplex_volume(s))
}


#' Affine transformation of a bbquad
#' 
#' @param x bbquad object
#' @param A transformation matrix
#' @param s shift vector
#' 
#' @export
bbquad_affine <- function(x, A, s=c(0,0,0)){
  d <- ncol(x$vertices)
  x$vertices <- coord_affine(x$vertices, A, s[1:d])
  x$volume <- bbquad_volume(x)
  x
}
