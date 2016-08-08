#' Intersection of two polytopes
#' 
#' Compute "a intesection b"
#' 
#' @param a Polytope a
#' @param b Polytope b
#' 
#' @export
bbpoly_intersection <- function(a, b) {
  if(!is.bbpoly(a)|!is.bbpoly(b)) stop("a and b should be bbpoly-objects.")
  # start with poly a and cut away using the planes of b
  ab <- a
  d <- ncol(a$vertices) 
  if(d!=3) stop("Intersection at the moment only for 3D polytopes.")
  # 
  planes <- bbpoly_planes(b)
  #
  
  for(iii in 1:ncol(planes)) {
    nn <- -planes[1:d, iii]
    np <- planes[1:d+d, iii]
    dist <- c(nn %*% np)
    skp <- ab$ve %*% nn
    # check which vertices are on which side of the plane
    th <- skp <= dist
    smallerDist <- which(th)
    largerDist <- which(!th)
    # does the plane hit the polytope
    hit <- length(smallerDist) & length(largerDist)
    #
    if(hit) {
      vert0 <- ab$ve[ ab$edges[,1], ]
      vert1 <- ab$ve[ ab$edges[,2], ]
      skp0 <- vert0%*%nn
      skp1 <- vert1%*%nn
      t0 <- skp0 - dist
      t1 <- skp0 - skp1
      s <- t0/t1
      # in each cut edge
      # replace vertex outside intersection by new vertex
      out0 <- which(skp0 > dist & s >= 0 & s <= 1)
      out1 <- which(skp1 > dist & s >= 0 & s <= 1)
      # edges where both vertices are out -> remove whole edge
      outEdges <- which( skp0 > dist & skp1 > dist)
      # indices of cut edges
      id <- c(out0, out1)
      #browser()
      nnewVertices <- length(id)
      nVertices <- nrow(ab$ve)
      # Fill position of old vertices by new vertices
      newVertices <- matrix(NA, nVertices + nnewVertices, d)
      newVertices[smallerDist,] <- ab$ve[smallerDist,]
      newVert <- as.vector(1-s)*vert0 + as.vector(s) * vert1
      #indices of new vertices (attached to the end)
      ind<-(nVertices+1):(nVertices+nnewVertices)
      #attach remaining new vertices to end
      newVertices[ind,] <- newVert[id,]
      # delete faces that are no longer needed
      outFaces <- which(sapply(ab$faces$vertices, 
                               function(face, outVertices)
                                 all(face %in% largerDist)
                               )
                        )
      if(length(outFaces>0)) {
        #delete faces
        ab$faces$vertices <- ab$faces$vertices[-outFaces]
        ab$faces$edges <- ab$faces$edges[-outFaces]
      }
      newEdges<-matrix(nrow=nrow(ab$edges), ncol=2)
      newEdges[out0,1]<-ind[1:length(out0)]
      newEdges[out1,2]<-ind[(length(out0)+1):length(ind)]
      if(length(out0)>0){
        newEdges[-out0, 1]<- ab$edges[-out0, 1]
      }
      else{
        newEdges[,1] <- ab$edges[, 1]  #set all edges
      }
      if(length(out1)>0){
        newEdges[-out1, 2] <- ab$edges[-out1,2]
      }
      else{
        newEdges[,2] <- ab$edges[,2]  #set all edges
      }
      newEdges[outEdges,]<-c(NA,NA)
      newFacesEdges <- ab$faces$edges
      newFacesVertices <- ab$faces$vertices
      #
      # set new vertices to faces
      for (j in 1:length(ab$faces$vertices)) {
        # edges in this face where vertex 1 is out
        j1<- which(ab$faces$edges[[j]]%in%out0)
        j2<- which(ab$faces$edges[[j]]%in%out1)
        if (length(j1)>0 | length(j2)>0) {
          # which labels do these vertices have?
          v1 <- ab$edges[ab$faces$edges[[j]][j1], 1]
          v2 <- ab$edges[ab$faces$edges[[j]][j2], 2]
          # these vertices should be replaced by these numbers
          k1<-newEdges[ab$faces$edges[[j]][j1], 1]
          k2<-newEdges[ab$faces$edges[[j]][j2], 2]
          # indices to be removed          
          n<-which( ab$faces$vertices[[j]] %in% c(v1,v2) )
          newFacesVertices[[j]] <- c(ab$faces$vertices[[j]][-n],k1,k2)
        }
      }
      if (nnewVertices == 3) { #new face is a triangle (easiest case)
        new<-t(combn(ind,2)) #any combination of new vertices forms an edge
        newEdges<-rbind(newEdges, new) #add new edges to edge list
        # include a new face
        # edges are the last three in the list
        n<- nrow(ab$edges)
        newFacesEdges[[length(ab$faces$edges)+1]] <- c(n+1, n+2, n+3)
        newFacesVertices[[length(ab$faces$vertices)+1]] <- ind
      }
      else { #here we need to test which vertices must be joined
        # all combinations of new vertices
        nEd <- t(combn(1:nnewVertices,2))
        n <- nrow(newEdges)
        newFaceVertices<-numeric()
        for (i in (1:nrow(nEd))) {
          # vector joining two vertices
          diff<-  newVertices[ind[nEd[i,1]],] - newVertices[ind[nEd[i,2]],]
          # compute plane spanned by normal vector of the cutting plane and diff 
          norm <- cross(nn, diff)
          # take any other new vertex 
          scalarProd <- apply(newVertices[ind[-c(nEd[i,1],nEd[i,2])],], 1, function(u) u%*%norm)
          dist<-(newVertices[ind[nEd[i,1]],]%*%norm)[1,1]
          # are all other vertices on the same side of the plane
          sameSide <- length(which(scalarProd>dist))== 0||(length(which(scalarProd<dist))== 0)
          # if all points are on one side -> generate edge  
          if(sameSide) {
            newEdges <- rbind(newEdges,c(ind[nEd[i,1]],  ind[nEd[i,2]]))
            if (! ind[nEd[i,1]] %in% newFaceVertices )
              newFaceVertices<-c(newFaceVertices, ind[nEd[i,1]] )
            if (! ind[nEd[i,2]] %in% newFaceVertices )
              newFaceVertices<-c(newFaceVertices, ind[nEd[i,2]] )
          }
        } 
        #generate a face
        newFacesEdges[[length(ab$faces$vertices)+1]]<- (n+1):nrow(newEdges)
        newFacesVertices[[length(ab$faces$vertices)+1]]<-newFaceVertices 
      } 
      
      # adjust faces
      #the edges of the new face
      edgesNewFace <- newEdges[newFacesEdges[[length(ab$faces$vertices)+1]], ]
      indNewEdges <- newFacesEdges[[ length(ab$faces$vertices) + 1 ]]
      # note: faces completely outside the intersection have already been removed
      for(j in 1:length(ab$faces$vertices)) {
        # variable id contains indices of cut edges is any of these contained in the face?
        cut <- which(id %in% ab$faces$edges[[j]]) 
        if (length(cut)>0)  {
          # cut edges have already been shortened and new vertices have been assigned
          # remove edges outside 
          indoutEdges<-which(newFacesEdges[[j]]%in%outEdges)
          if ( length(indoutEdges) > 0 )
            newFacesEdges[[j]] <- newFacesEdges[[j]][-indoutEdges]
          
          # assign edges of the intersection face
          # these should link two vertices of the face
          #add an edge to the current face i if both of its vertices are in face i
          index <- which(apply(edgesNewFace, 1, 
                             function(edge){ 
                               all(edge %in% newFacesVertices[[j]])
                               }))
          newFacesEdges[[j]] <- c(newFacesEdges[[j]], indNewEdges[index])
        }
      } #else face is not cut
      # faces which are completely outside have already been deleted the others stay as they are
      ab$vertices <- newVertices 
      ab$edges <- newEdges
      ab$faces$edges <- newFacesEdges
      ab$faces$vertices <- newFacesVertices
    } # hit
  } # for( faces )
  # re-index
  
  map <- ab$vertices[,1]
  mape <- ab$edges[,1]
  na <- is.na(ab$vertices[,1])
  map[!na] <- 1:sum(!na)
  nae <- is.na(ab$edges[,1])
  mape[!nae] <- 1:sum(!nae)
  ab$vertices <- ab$vertices[!na,]
  # mapping to new index
  mapf <- function(v, map)  apply(v, 2, function(j) map[na.omit(j)])
  mapfv <- function(v, map)  map[na.omit(v)]
  ab$edges <- mapf(ab$edges, map)
  ab$faces$vertices <- lapply(ab$faces$vertices, mapfv, map=map)
  ab$faces$edges <- lapply(ab$faces$edges, mapfv, map=mape)
  # fix the face facing
  out <- bbpoly_correct_faces(ab)
  out$volume <- NULL #bbpoly_volume(out)
  out
}

#' Split a polytope to simplices.
#' 
#' @param x bbpoly object
#' @param as.bbpoly Return the simplices as bbpoly-objects?
#' @details Assumes convexity, will not check.
#' 
#' @export
bbpoly_simplices <- function(x, as.bbpoly = FALSE) {
  check_bbpoly(x)
  v <- x$vertices
  d <- ncol(v)
  n <- nrow(v)
  f <- if(as.bbpoly) as_bbpoly else identity 
  s0 <- bbpoly_simplex(d = d)
  if(d==2){
    # use the first vertex as central
    edg <- x$edges
    tri <- triangulate_indices_cohen_hickey(edg)
    slist <- apply(tri, 1, function(i) 
              f(list(vertices=v[i,],
                     edges = s0$edges))
    )
  }
  else if(d==3){
    s <- list()
    # the central vertex
    c0 <- apply(x$vertices, 2, mean)
    e0 <- s0$edges
    f0 <- s0$faces
    # triangulate each face
    # split each face into tetrahedrons
    k <- 0
    for(e in x$faces$edges){
      edg <- x$edges[e,]
      tri <- triangulate_indices_cohen_hickey(edg)
      for(i in 1:nrow(tri)){
        j <- tri[i,]
        sp <- f(list(vertices = rbind(c0, v[j,]),
                   edges = e0,
                   faces = f0))
        s[[k <- k + 1]] <- sp
      }
    }
    s # done
  }
  else stop("bbpoly not 2d or 3d.")
}

#' Check and correct for orientation of faces
#' 
#' Assumes faces with vertices on the same plane
#' 
#' @param x bbpoly object in 3D
#' 
#' @export
bbpoly_correct_faces <- function(x) {
  check_bbpoly(x)
  d <- ncol(x$vertices)
  if(d == 2) stop("x not in 3d")
  c0 <- apply(x$vertices, 2, mean)
  pl <- bbpoly_planes(x)
  flip <- NULL
  for(i in 1:ncol(pl) ) {
    dist <- (c0-pl[1:d+d,i]) %*% pl[1:d, i]
    flip <- c(flip, dist < 0)
    if(dist < 0) x$faces$vertices[[i]][1:3] <- rev(x$faces$vertices[[i]][1:3])
  }
  x
}

#' Volume of a bbpoly
#' 
#' @param x bbpoly object
#' 
#' @details splits the bbpoly into simplices and sums  their volume.
#'   
#' @export

bbpoly_volume <- volume.bbpoly <- function(x) {
  s <- bbpoly_simplices(x)
  # get only the vertices
  s <- lapply(s, getElement, "vertices")
  sum(simplex_volume(s))
}

#' Volume of a simplex
#' 
#' @param x simplex, as a (d+1,d) matrix of vertices where d is dimension
#' 
#' @details 
#' If x is a list, we apply the function to each element.
#' 
#' @export
simplex_volume <- function(x) {
  if(is.list(x) & ! is.bbpoly(x)) sapply(x, simplex_volume)
  else{
    if(is.bbpoly(x)) return( simplex_volume(x$vertices) )
    d <- nrow(x)-1
    if(d != ncol(x)) stop("Can't interpret the simplex.")
    A <- t(x[-1,])-x[1,]
    abs(det(A)/factorial(d))
  }
}

