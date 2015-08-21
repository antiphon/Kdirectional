#' compute the bounding box for coordinates
#'
#' @param x coordinate matrix
#' @export
bbox_make <- function(x) {
  apply(as.matrix(x), 2, diff)
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


#' Affine transformation of a bbox
#' 
#' @export
bbox_affine <- function(bbox, A, s=c(0,0,0)){
  if(is.bbquad(bbox)){
    d <- ncol(bbox$vertices)
    bbox$vertices <- bbox_affine(bbox$vertices, A, s[1:d])
    bbox$volume <- bbox$volume * abs( det(A) )
  }
  else{
    d <- ncol(bbox)
    bbox <- t(t(bbox %*% t(A)) + s[1:d])
  }
  bbox
}


#' Bounding box volume
#' 
#' @export
bbox_volume <- function(bbox){
  if(is.bbquad(bbox)){
    if(is.null(bbox$volume)) stop("Can't get the volume of the box.")
    bbox$volume
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


