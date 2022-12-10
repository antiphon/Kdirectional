#' Testing contrast for ellipses
#' 
#' Give this to confint.ellipsoid to compute the contrast for
#'  testing semi_a - semi_b = 0
#' @param out do outlier removal, the quadratic simulations is unstable.
#' @param range range of outlier removal. See boxplot for details of the formula.
#' 
#' 
#' @export
ellipse_contrast_2d <- function(els, out=TRUE, range=1.85){
  
  if(!ncol(els[[1]]$A)==2) stop("Only for 2D.")
  
  R <- lapply(lapply(els, getElement, "A"), ellipse_solve_rota)
  axes <- t(sapply(R, function(b) b$axes ))
  angs <- sapply(lapply(R, getElement, "R"), function(R) {f <- R %*% c(1,0)
  atan2(f[2],f[1])})
  angs[angs <0 ] <- angs[ angs<0 ] + pi
  ab <- sapply(1:length(angs), function(i) 
    if(angs[i] < 3*pi/4 & angs[i] > pi/4) axes[i,1:2] else axes[i,2:1]  )
  
  e<- ab[1,]-ab[2,]
  e <- e[!is.na(e) & !is.infinite(e)]
  # check outliers
  if(out){
    m <- median(e, na.rm=TRUE)
    out <- abs(e-m) > range * abs(quantile(e-m, prob=3/4))
    e<-e[!out]
  }
  e
  
}

#' Testing contrast for ellipsoids
#' 
#' Give this to confint.ellipsoid to compute the contrast for
#'  testing semi_i - 0.5 * (semi_j + semi_k) = 0
#' 
#'  @param els list of ellipsoid estimates
#'  @param ijk on of c(1,2,3), c(2,1,3) or c(3,1,2). Default: c(3,1,2)
#'  @param out if >0 Use outlier detection to stabilise (useful with large variances)
#' 
#' @details 
#' 
#' The default contrast will compare third semi-axis to the average of first two. This seems to be the best, as the semi-axes are returned in an increasing order by the OLS-algorithm.
#' 
#' The 'out' is used like 'coef' in boxplot.stats, see there for details.
#' @export
ellipsoid_contrast_3d <- function(els, ijk = c(3,1,2), out=TRUE, z0=1.85){ # this is best
  #
  if(!all(1:3%in%ijk)) stop("ijk should be one of c(1,2,3), c(2,1,3) or c(3,1,2)")
  
  # determine what happends to x,y and z axis
  axes <- sapply(els, function(e){r<- diag(e$A);r/prod(r)} )
  #
  # the contrast
  e <- axes[ijk[1],]-0.5*(axes[ijk[2],]+axes[ijk[3],])
  
  if(out){
    m <- median(e, na.rm=TRUE)
    out <- abs(e-m) > z0 * abs(quantile(e-m, prob=3/4, na.rm=TRUE))
    e<-e[!out]
  }
  na.omit(e)
}
