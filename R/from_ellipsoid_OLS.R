# From ellipsoid-package
#
#############################################################
#' Ordinary Least Squares Ellipsoid
#' 
#' OLS fitting of an ellipsoid to a d-dimensional point cloud
#'
#' @param x point coordinate matrix
#' @param origin Fix center of the ellipsoid to be the origin?
#' 
#' @details Assumes the normal error model.
#' 
#' @export
ellipsoid_OLS  <- function(x, origin=FALSE) {
  n <- nrow(x)
  d <- ncol(x)
  if(n < 3+d) stop("Trying to fit an ellipsoid to too little amount of points.")
  D  <- matrix(2, ncol=d, nrow=d)
  diag(D) <- 1
  D <- D[upper.tri(D,T)]
  ## Y:
  Y<-apply(x, 1, function(x){
    X  <- x%*%t(x)
    vX <- X[upper.tri(X, T)]
    y <- c(D*vX, if(origin) NULL else c(x), 1)
    y
  })
  Y <- t(Y)
  #
  nd <- (d*(d+1)/2)
  nb <- nd + (d * !origin) + 1 
  H <- diag(1, nb)
  Hi <- solve(H)
  # svd
  USV <- svd(Y%*%Hi)
  # the estimate of parameters:
  beta <- c( Hi%*%USV$v[,nb] )
  if(origin){ # add center parameters to comply with other functions
    beta <- c(beta[1:nd], rep(0, d), beta[(nd+1):nb])
  }
  
  # convert to center and matrix -formulation
  res <- ellipsoid_from_beta(beta, d)
  res$ndata <- n
  ####### add the error variance
  # the error variance. Need to center the coords, note how this is not taken into account.
  
  # center
  xc <- t( t(x) - res$center)
  
  r_data <- sqrt(rowSums(xc^2))
  u_data <- xc/r_data
  r_pred <- predict(res, u_data)
  resi <- r_data - r_pred
  n <- sum(!is.na(resi))
  if(d==2){
    s2 <- sum(resi^2, na.rm=TRUE)/(n-1) # this is approximation with 0 curvature TODO better
  }
  else{
    s2 <- sum(resi^2, na.rm=TRUE)/(n-1) # this is approximation to 0 curvature TODO BETTER
  }
  # the variance covariance
  S0 <- solve_S0(Y)
  S <- s2 * S0
  #
  if(origin){ # add dummy for center
    S0 <- diag(1e-9, nb+d)
    S0[1:(nb-1), 1:(nb-1)] <- S[-nb,-nb]
    aa <- c(1:(nb-1), nb+d)
    S0[aa, nb+d] <- S0[nb+d, aa] <- S[nb,]
    S0 <- 0.5 * (S0 + t(S0))
    S <- S0
  }
  # make sure symmetric varcov
  m <- max(S-t(S))>0.01
  
  if(is.na(m)) browser()
  if(m) 
    warning("varcov of quadratic equation terms asymmetric, forcing symmetric.")
  S <- 0.5 * (S + t(S))
  #
  res$ols_fit <- list(varcov=S, beta_est=beta, s2=s2)
  res$origin <- origin
  #
  # done
  res
}


#' Solve covariance matrix from unstable X
#' 
#' @param Y The matrix, nrow < ncol
#' 
#' @details Solves the inverse of t(Y)%*%Y by slowly increasing diagonal (hedging).
solve_S0 <- function(Y){
  warn <- FALSE
  m <- 1e-9 # hedging
  ss <- try(S0 <- solve(t(Y)%*%Y + diag(m, ncol(Y))), TRUE)
  while("try-error"%in% is(ss)){
    ss <- try(S0 <- solve(t(Y)%*%Y + diag(m <- m * 2, ncol(Y))), TRUE)
    warn <- TRUE
  }
  if(warn) warning(paste("Hedging needed for solving covariance matrix. (diag ", m, ")"))
  S0
}
