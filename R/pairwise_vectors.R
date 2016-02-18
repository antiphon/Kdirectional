#' Compute pairwise distances and angles
#'
#' Computes pairwise vectors of a 2d or 3d point pattern. Returns 
#' circular- or spherical coordinates, depending on the dimension.
#'
#' @param x matrix of coordinates.
#' @param from Indices from which to compute
#' @param to Indices to which to compute
#' @param asMatrix return n x n matrices
#' @param as.xyz return as xyz-coordinates instead of polar/spherical
#' @details
#' 
#' Default is to return the upper triangle (without diagonal) of 
#' pairwise distance and angle matrices.
#' 
#' Result will have 2 or 3 columns, where first is the distance
#' and 2(:3) are angle(s).
#'   
#' @export

pairwise_vectors <- function(x, from, to, asMatrix=FALSE, as.xyz=FALSE){
  x <- as.matrix(x)
  d <- ncol(x)
  #
  # check for subsetting
  if(!missing(from) | !missing(to)){
    if(missing(from)) from <- 1:nrow(x)
    if(missing(to)) to <- 1:nrow(x)
    U <- c_pairwise_dist_angle_subset(x, from, to)
    # drop empty ones
    nonempty <- U[,1]>0
    U <- U[nonempty,]
  }
  else{ # all pairs once
    U <- c_pairwise_dist_angle(x)
  }
  colnames(U) <- c("distance", if(d==2) "angle" else c("azimuth", "inclination"))
  
  # if as matrix
  if(asMatrix){
    if(!missing(from)|!missing(to)) stop("Not implemented for subsets.")
    R <- list()
    M0 <- A <- D <- diag(0, nrow = nrow(x))
    ltri <- lower.tri(D)
    # distance matrix
    D[ltri] <- U[,1]
    D <- D + t(D)
    # angles:
    if(d==3) A2 <- A
    # i -> j, transposing
    A[ltri] <- U[,2]
    # j->i = antipode
    B <- M0
    Aa <- U[,2]
    i <- U[,2]>pi
    Aa[i] <- Aa[i]-pi
    Aa[!i] <- Aa[!i]+pi
    B[ltri] <- Aa
    A <- B+t(A)
    diag(A) <-NA
    if(d==3) {
      A2[ltri] <- U[,3]
      B2 <- M0
      B2[ltri] <- pi-U[,3]
      A2 <- B2+t(A2)
      diag(A2) <-NA
    }
    
    diag(D) <- 0
    R$distance <- D
    if(d==2)R$angle <- A
    else{R$azimuth <- A; R$inclination <- A2}
  }
  else R <- U
  if(!asMatrix & as.xyz){
    if(d==2) R <- R[,1] * cbind(cos(R[,2]), sin(R[,2]))
    else R <- R[,1] * ai2xyz(R[,-1])
  }
  R
}

