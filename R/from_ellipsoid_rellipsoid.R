#'Sample points Roughly Uniformly on the surface of an ellipsoid, 2D or 3D
#'
#'@param axes vector giving the semi-axis lengths
#'@param R rotation matrix
#'@param center center vector, c(x,y) or c(x,y,z)
#'@param noise.sd noise sd, 0=no noise
#'@param pieces Number of pieces to divide the surface into. See details.
#'  
#'@details 
#'
#' Sampling is done by first sampling uniformly on the surface of a
#' sphere and then transforming with the ellipsoid axes. Corrections:
#'
#' * 2D: The sampling weights are set to piecewise linear approximation of the ellipse surface integral, 
#' the resolution of which is controlled by 'pieces'.
#' * 3D: No weighting at the moment, so very biased for elongated ellipsoids.
#'
#' Additional nD Gaussian noise can be added by setting `noise.sd`>0.
#'@export 

rellipsoid <- function(n, axes=c(1,1,1), center, noise.sd=0, R=NULL, pieces=1000){
  d <- length(axes) 
  if(d > 3 | d < 1) stop("Only 2D or 3D ellipsoids are supported.") 
  if(missing(center)) center <- rep(0, d)
  if(length(center) != d) stop(paste0("'center' should be a vector of length ", d))
  M <- diag(axes)
  # 2D
  if(d==2){
    # linear approximation to the surface integral
    theta <- seq(0, 2 * pi, length=pieces)
    xy <- cbind(cos(theta), sin(theta))%*%M
    # the surface integrals ~ distance between consecutive points
    integ <- rowSums((xy[-1,]-xy[-pieces,])^2)/(diff(theta[1:2]))
    x <- NULL
    # for some reason this works better than just sample(1:pieces, n, prob=integ)
    while(length(x) < 2*n){
      i <- sample(1:(pieces-1), prob = integ)
      u <- sapply(i, function(j) runif(1, theta[j], theta[j+1]) )
      x <- rbind(x, cbind(cos(u), sin(u)))
    }
    x <- x[1:n,]
    # transform
    x <- x%*%M
  }
  # 3D is still biased
  else if(d==3){
    # stretch sphere sample, will be biased
    u <- runif(n)
    v <- runif(n)
    a <- u * 2 * pi
    i <- acos(2*v - 1)
    ai<- cbind(azi=a, inc=i)
    x <- ai2xyz(ai)
    x <- x%*%M
    warning("3D rellipsoid is truely uniform sample only for a sphere.")
  }
  
  #
  # Rotate if needed
  y <- if(!is.null(R)) t(R%*%t(x)) else x
  #
  # Add noise if needed 
  if(noise.sd) y <- y + rnorm(nrow(y)*d, 0, noise.sd)
  # shift
  y <- t( t(y) + center)
  y
}









#' Sample points uniformly on the surface of an ellipsoid, 2D or 3D
#' 
#' @param axes vector giving the semi-axis lengths
#' @param R rotation matrix
#' @param noise.sd noise
#' 

rellipsoid_dev <- function(n, axes=c(1,1,1), noise.sd=0, R=NULL, method=1, rej=FALSE){
  # transform points on the surface of an circle
  # points on a sphere/circle
  d <- length(axes)
  if(d > 3 | d < 1) stop("Only 2D or 3D ellipsoids are supported.") 
  M <- diag(axes)
  
  ## Normal sample method (wrong)
  #if(0){
  #     rn <- rmvnorm(n, rep(0, d), M)
  #     # project onto surface
  #     u <- rn/sqrt(diag(rn%*%t(rn)))
  #     r <- sqrt( 1/diag(u%*%M%*%t(u))  )
  #     x <- r*u
  #}
  # 
  # 2D
  if(d==2){
    if(method == 1){
      # linear approximation to the surface integral
      stops <- 1000 + 1
      theta <- seq(0, 2 * pi, length=stops)
      xy <- cbind(axes[1] * cos(theta), axes[2] * sin(theta))
      # the surface integrals ~ distance between consecutive points
      integ <- rowSums((xy[-1,]-xy[-stops,])^2)
      #
      # this gives us a sampling acceptance rate for the angle
      den <- integ/(diff(theta[1:2]))
      integf <- function(u) den[cut(u, theta, labels=FALSE)]
      #
      x <- NULL
      while(length(x) < 2*n){
        m <- floor((n - length(x)/2)*2)
        if(rej){
          u <- runif(m, 0, 2*pi)
          pro <- integf(u)
          keep <-  runif(m) < pro
          u <- u[keep]
        }
        else{
          i <- sample(1:(stops-1), prob=den)
          u <- sapply(i, function(j) runif(1, theta[j], theta[j+1]) )
        }
        x <- rbind(x, cbind(cos(u), sin(u)))
      }
      x <- x[1:n,]
      # transform
      x <- x%*%M
    }
    else{ # the naive wrong approach, for comparison
      u <- runif(n, 0, 2*pi)
      x <- cbind(cos(u), sin(u))%*%M
      warning("This is not true uniform sample except for a sphere.")
    }
  }
  #
  # 3D
  else if(d==3){
    # piecewise linear approximation to the surface integral
    tri <- ellipsoid_shape(2)
    # stretch sphere sample, rejection sampling
    weights <- triangulation_areas(tri)
    probs <- weights/max(weights)
    # drop the homog. 
    vb <- (t(tri$vb[1:3,])/tri$vb[4,])
    nvertex <- nrow(vb)
    #
    x <- NULL
    nx <- n #??#min(n, 5000) # to not overload memory 
    while(length(x) < 3*n){
      xn <- runifsphere(nx)
      # culling
      gu <- as.matrix(dist(rbind(vb,xn)))[1:nvertex,1:nx+nvertex]<0.3 #geom(rbind(xn, vb), from=1:nvertex+nx, to=1:nx, r=0.3)
      # Thin according to weights
      xkept <- sapply(1:ncol(tri$it), function(i){
        ijk <- tri$it[,i]
        # cull
        ok <- apply(gu[ijk,], 2, any) #unique(unlist(gu[ijk+n]))
        inside <-ok[ which(check_vector_passes_triangle(xn[ok,], vb[ijk,])) ]
        kept <- inside[ runif(length(inside)) < probs[i] ]
        if(length(kept)) xn[kept,]
        else NULL
      })
      x <- rbind(x, do.call(rbind, xkept))
    }  
    x <- x[sample(nrow(x), n),] %*% M
  }
  
  #
  # Rotate if needed
  if(is.null(R)) R <- diag(1, d)
  y <- t(R%*%t(x))
  #
  # Add noise if needed 
  if(noise.sd) y <- y + rnorm(nrow(y)*d, 0, noise.sd)
  y
}




#' Check if a vector (or vectors), multiplied by inf, go(es) through a triangle in 3D
#' 
#' @param u unit vector
#' @param triangle triangle coordinates
#' 
#' @import sp
#' @export

check_vector_passes_triangle <- function(u, triangle){
  # project
  pl <- three_points_to_plane(triangle)
  
  proj <- project_to_plane(u, pl[1:3])
  tproj <- project_to_plane(triangle, pl[1:3])
  # check
  ( point.in.polygon(proj[,1], proj[,2], tproj[,1], tproj[,2]) == 1 )
}



