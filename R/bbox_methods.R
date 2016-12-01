#' compute the bounding box for coordinates
#'
#' @param x coordinate matrix
#' @param expand expand a bit? (ripras, see details)
#' 
#' @details if expand = TRUE, the bbox will be expanded using an estimate of the true window area
#' as described by Ripley & Rasson (1977). The convex hull is replaced with a cuboid aligned with the axes.
#' 
#' @export
bbox_make <- function(x, expand=TRUE) {
  bb <- apply(as.matrix(x), 2, range)
  # expand a bit
  if(expand){
    n <- nrow(x)
    d <- ncol(x)
    m <- 2^d # number of vertices
    # Estimate the true volume
    #Vbb <- bbox_volume(bb)
    #Vtrue <-  Vbb / (1 - m/n)
    # need to extend the volume
    #ex <- Vtrue/Vbb
    ###
    # This if for convex hull
    ex <- 1 / min(1, max(0.01, (1-m/n)))
    ex_d <- ex^(1/d)
    # per dimension
    apply(bb, 2, function(ab) {
        r <- diff(ab) * (ex_d-1)
        ab + ex_d * c(-1,1) * r/2
      }  )
  } 
  else bb
}

#' Turn spatstat window object to bbox
#' 
#' @export
owin_to_bbox <- function(x) {
  if(x$type == "rectangle") cbind(x$xrange, x$yrange)
  else if(x$type == "polygonal")
    with(x$bdry[[1]], poly_to_bbquad(x,y)) # just one pieace of the window supported
  else stop("Unable to interpret the window.")
}


#' distance from points to bbox border
#'
#' @param x coordinate matrix
#' @param bbox bounding box 
#' @export

bbox_distance <- function(x, bbox){
  if(is.bbquad(bbox)) {
    bbquad_distance(x, bbox)
  }
  else{
    dim <- ncol(bbox)
    di <- function(loc){
      d <- sapply(1:dim, function(i) min(loc[i]-bbox[1,i], abs(loc[i]-bbox[2,i])))
      min(d)
    }
    c(apply(x, 1, di))
  }
}


#' Affine transformation of a bbox or bbquad
#' 
#' @param bbox bounding box
#' @param A transformation matrix
#' @param s shift vector
#' 
#' @details 
#' If A is not diagonal, bbox will be transformed into a bbquad.
#' 
#' @export
bbox_affine <- function(bbox, A, s=c(0,0,0)){
  if(is.bbquad(bbox)){
    bbox <- bbquad_affine(bbox, A, s)
  }
  else{ # 
    if(any(c(A[upper.tri(A)],A[lower.tri(A)])!=0)) {
      bbox <- bbox2bbquad(bbox)
      bbox <- bbquad_affine(bbox, A, s)
    }
    else {
      bbox <- coord_affine(bbox, A, s)
    }
  }
  bbox
}

#' Affine transfomration of coordinates
#' 
#' @param x coordinates, one row per coordinate
#' @param A matrix
#' @param s shift

#' @export
coord_affine <- function(x, A, s=c(0,0,0)){
  x <- rbind(x)
  d <- ncol(x)
  t(t(x %*% t(A)) + s[1:d])
}



#' Bounding box volume
#' 
#' @export
bbox_volume <- function(bbox){
  if(is.bbquad(bbox)){
    if(is.null(bbox$volume)){
      bbquad_volume(bbox)  
    }
    else bbox$volume
  }
  else{
    prod(apply(bbox, 2, diff))
  }
}

#' Bounding box side lengths
#' 
#' @export
bbox_sideLengths <- function(bbox){
  if(is.bbquad(bbox)){
    d <- ncol(bbox$vertices)
    edges <- bbquad_edges(bbox)
    e1 <- edges[1:d,]
    e2 <- edges[1:d+d,]
    sqrt(colSums((e2-e1)^2 ))
  }
  else{
    apply(bbox, 2, diff)
  }
}

#' Bounding box dimension
#' 
#' @param x Bounding box or bbquad
#' @export

bbox_dim <- function(x){
  if(is.bbquad(x)) ncol(x$vertices)
  else ncol(x)
}