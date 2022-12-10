#' ALS fitting of ellipsoid
#' 
#' Adjusted LS fitting of ellipsoids, by Kukush et al. 2003.
#' 
#' @param x Point locations, n x d matrix
#' @param s2 Error variance, or a vector over which to optimize.
#' @param ... ignored
#' @param just_U For internal use: Return just smallest eigenvector and bhat.
#'
#' @details 
#' Not really optimised yet. Good grid for s2 optimisation is something like
#' seq(s2 * 0.1, s2, length=20)
#' where s2 is the OLS estimate (access via $ols_fit$s2).
#'
#' @references 
#' Kukush et al (2003): Consistent estimation in an implicit quadratic measurement error model,
#' Computational Stat & Data Analysis
#' 
#' @import parallel
#' @export

ellipsoid_ALS <- function(x, s2, ..., just_U=FALSE) {
  t0 <- function(v) v*0+1
  t1 <- function(v) v
  t2 <- function(v) v^2-s2
  t3 <- function(v) v^3 - 3*v*s2
  t4 <- function(v) v^4-6*v^2*s2+3*s2^2
  # 
  if(missing(s2)) stop("Provide error variance or a vector over which to optimize for it.")
  #
  X <- t(x)
  d <- nrow(X)
  n <- ncol(X)
  fl <- list(t0,t1,t2,t3,t4)
  # vectorisation operator
  vec <- function(A) A[upper.tri(A,TRUE)]
  #
  # In case s2 not given
  if(length(s2)>1){
    # For estimating the s2 we need to iterate
    U <- mclapply(s2, function(s2) ellipsoid_ALS(x, s2=s2, just_U=TRUE))
    u <- sapply(U, getElement, "U")
    i <- which.min(u^2)
    hatb <- U[[i]]$hatb
    U <- cbind(s2=s2, U=u)
    s2 <- s2[i]
  }
  #1 tensors
  Tkil <- lapply(fl, function(f) f(X))
  #2 1,i,M
  v1 <- rep(1, d+1)
  vi <- c(1:d,0)
  nb <- (d+1)*d/2+d+1
  M <- cbind(vec(vi%*%t(v1)), vec(v1%*%t(vi)))
  # skip 3
  # 4 - 6 Psi
  D1 <- 1:( (d+1)*d/2 )
  D2 <- (1:d)*(1:d+1)/2
  D <- setdiff(D1, D2)
  eta <- diag(0, nb)
  A <- matrix(0, nrow=d, ncol=n)
  Psi <- diag(0, nb)
  for(p in 1:nb)
    for(q in p:nb) {
      for(l in 1:n)for(i in 1:d){
        R <- 1+(M[p,1]==i) + (M[p,2]==i) + (M[q,1]==i) + (M[q,2]==i)
        A[i,l] <- Tkil[[R]][i,l]
      }
      eta[p,q] <- sum(apply(A, 2, prod))
      k <- if(p%in%D & q%in% D) 4 else if((!p%in%D)&(!q%in%D) ) 1 else 2
      Psi[p,q] <- Psi[q,p] <- eta[p,q] * k
    }
  # 7 eigen
  ee <- eigen(Psi)
  hatb <- ee$vectors[, which.min(ee$values)]
  # 8 normalize
  hatb <- hatb/sqrt(sum(hatb^2))
  U <- min(ee$values)/n
  if(just_U) return(list(U=U, hatb=hatb))
  #
  # Construct ellipsoid object
  als <- ellipsoid_from_beta(hatb, d)
  # add info
  als$ndata <- n
  # covariance matrix
  S <- s2 * solve_S0(Psi)
  als$ols_fit <- list(s2=s2, U=U, beta_est=hatb, varcov=S)
  # done
  als
}

