#' Plot an Ellipsoid Object
#' 
#' Internal ellipsoid plotter.
#' 
#' @param x ellipsoid
#' @param add logical
#' @param i Normal 1,2 or 3. Plot intersection of 3d ellipsoid.
#' @param levels If i>0, draw this many levels from origin to tip of the axis
#' @param res 2d resolution
#' @param scale 1. use to enlarge or shrink
#' @param N 2, iteration of refinement in 3D (resolution)
#' @param xlab axis-label 
#' @param ylab axis-label 
#' @param ... passed on to lines or rgl::shade3d
#' 
#' @details 
#' #' For 3D ellipsoid, if i=1,2 or 3, plot 2D intersection of 
#' plane:  i = 1 yz-plane; i = 2 xz-plane; i = 3 xy-plane. 
#' 
#' @importFrom rgl shade3d
#' 
#'@export
plot.ellipsoid <- function(x,
                           add=TRUE, i=0, levels=1, 
                           res=128, scale=1, N=2, 
                           xlab="x",
                           ylab="y", ...){
  if(x$dim==2){
    a <- c(seq(0, 2*pi, length=res))
    y <- cbind(cos(a),sin(a))
    z <- y * predict(x, y) * scale
    z <- t(t(z) + x$center)
    if(!add) plot(NA, xlim=range(z), ylim=range(z),
                  asp=1, xlab=xlab, ylab=ylab)
    lines(z, ...)
  }
  else{
    if(i!=0){
      if(xlab!="") xlab <- c("y", "x", "x")[i]
      if(ylab!="") ylab <- c("z", "z", "y")[i]
      # max
      n <- c(0,0,0); n[i] <- 1
      rr <- c(0,0,0); rr[(i+1)%%3+1] <- 1
      maxq <- predict(x, rbind(n))
      qseq <- seq(0, maxq, l=levels)
      # max reach
      for(j in 1:levels){
        q <-  x$center
        q[i] <- q[i] + qseq[j]
        p <- intersect_ellipsoid_plane(x, n, q)
        basis <- p$basis
        r <- basis[-i,1]
        s <- basis[-i,2]
        cent <- p$center3d[-i]
        # rotation
        ang <- atan2(r[2],r[1])
        R <- rotationMatrix3(az=ang)[-3,-3]
        # here we have the ellipse
        el2 <- as_ellipsoid(p$semi_axes, R, center=cent)
        ### plot it
        plot.ellipsoid(el2, add=add | j>1, 
                       res=res, scale=scale, 
                       xlab = xlab, ylab = ylab,
                       ...)
      }
    }
    else{
      y <- ellipsoid_shape(N=N, axes=x$semi_axes, R=x$rot, center=x$center)
      rgl::shade3d(y, ...)
    }
  }
}

#' ##################################################################
#' Predict i.e. give the length of a direction to be on the ellipsoid
#'
#' @param object ellipsoid 
#' @param u direction into which predict the distance to the surface of ellipsoid
#'
#' @returns  Returns the distance from ellipsoid center to the ellipsoid surface in
#' the given directions.
#'
#'@export
predict.ellipsoid <- function(object, u, ...){
  x <- object
  if(missing(u)) stop("direction(s) u needed")
  d <- 1 / diag(u %*% x$A %*% t(u) )
  r <- sqrt(abs(d))
  r
}


#' Ellipsoid shape for 3d plotting
#'
#' Refine an icosahedron, then transform
#' 
#' @param R Rotation matrix.
#'
#' @importFrom rgl icosahedron3d subdivision3d translate3d asHomogeneous
#' @export
ellipsoid_shape <- function(N=2, axes=c(1,1,1), R=NULL, center=c(0,0,0)){
  ico <- rgl::icosahedron3d()
  for(i in 1:N) ico <- rgl::subdivision3d(ico)
  xy <- t(ico$vb[-4,])
  # units
  xy <- xy/sqrt(rowSums(xy^2))
  D <- diag(axes)
  
  xy <- xy %*% D
  exy <- if(is.null(R)) xy else t(R%*%t(xy))
  ico$vb <-t(  asHomogeneous(exy)  )
  translate3d(ico, center[1], center[2], center[3])
}


#' ###########################################################
#' Print ellipsoid
#' 
#' @export
print.ellipsoid <- function(x, ...){
  type <- ifelse(x$dim==2, "2D ellipse", "3D ellipsoid")
  if(!is.null(x$ave))
    cat(paste0("Average ", type, ", computed from ", x$nellipses, " ", type, "s.\n"))
  else if(!is.null(x$n)) cat(type, "fitted to", x$n, "points.\n")
  else cat(type, "\n")
}

#' ###################################################################
#' Ellipse center and matrix from general parameter form
#' 
#' @param beta OLS estimates
#' @param d dimension
#' @param check Check for definiteness?
#' 
#' @export
ellipse_form <- function(beta, d, check=FALSE){
  nd <- (d*(d+1)/2)
  nb <- nd + d + 1 
  # ellipse parameters
  A <- diag(0, d)
  A[upper.tri(A,T)] <- beta[1:nd]  
  A[lower.tri(A)] <- A[upper.tri(A)]
  b <- beta[(nd+1):(nd+d)]
  dhat <- beta[nb]
  m <- 0
  v <- try(Sa <- solve(A + diag(m, ncol(A))))
  while("try-error"%in% is(v))  v <- try(Sa <- solve(A + diag(m<-m+5e-8, ncol(A))))
  chat <- -0.5 * Sa%*%b
  Ahat <- A / (t(chat) %*% A %*% chat - dhat)[1]
  # check definiteness
  if(check){
    e <- eigen(Ahat)
    if(any(e$values < 0)){
      i <- which(e$values > 0)
      Av <- diag(0, ncol(Ahat))
      for(j in i) Av <- Av + e$val[j] * e$vec[,j]%*%t(e$vec[,j])
      Ahat <- Av
    }
  }
  list(c=chat, A=Ahat)
}

#' ########################################################################
#' Solve rotation and semi-axes from general transform
#' 
#' @param A the trasformation matrix in an ellipsoid equation
#' @param eps Inflate diagonal by this factor, to avoid numerical problems.
#' 
#' @details does not work with improper A, will return inf long axes.
#'   
#' @export
ellipse_solve_rota <- function(A, eps = 0){
  ev <- eigen(A + diag(diag(A)*eps))
  # make sure working with a definite. But if improper, dont go as breaks:
  if(any(ev$value<0) && !all(ev$value <= 0)){
    i <- which(ev$value > 0)
    S <- diag(0, ncol(A))
    for(j in i) S <- S + ev$value[j] * ev$vec[,i]%*%t(ev$vec[,i])
    ev <- eigen(S)
  }
  # in case negative eigenvalues persist, they are
  # extremely close to -0 so we might as well take abs here.
  axes_len <- 1/sqrt(abs(ev$value)) # the semi-axes lengths
  R <- ev$vector # rotation
  list(axes=axes_len, R=R)
}


#' #########################################################################
#' Sample from the OLS estimate of the beta parameters
#' 
#' Simulate the beta parameters, approximately normal conditional on ||beta||^2=1
#' 
#' @param x ellipsoid object
#' @param nsim number of simulations
#' @param tol tolerance when comparing to vector length 1
#' @param maxiter maximum number of iterations to try to get enough samples within tolerance
#' 
#' @importFrom mvtnorm rmvnorm 
#'
#' @export
sample_ellipse_beta <- function(x, nsim=100, tol=0, maxiter=500){
  d <- x$dim
  nb  <- (d*(d+1)/2) + d + 1
  
  b <- rmvnorm(nsim, x$ols_fit$beta_est, x$ols_fit$varcov)
  if(maxiter==0){ 
    tol <- 0
    b <- b/sqrt(rowSums(b^2))
  }
  if(tol>0){
    dev <-  abs(sqrt(rowSums(b^2))-1)
    ok <- dev < tol
    b <- b[ok,]
    it <- 0
    failed <- FALSE
    while(length(b)/nb < nsim){
      b <- rbind(b, rmvnorm( 2*(nsim-length(b)/nb), 
                             x$ols_fit$beta_est, 
                             x$ols_fit$varcov))
      ok <- abs(sqrt(rowSums(b^2))-1) < tol
      b <- b[ok,]
      it<-it+1
      if(it>maxiter) {it<-0; tol <- tol * 10; failed<-TRUE}
    }
    if(failed) warning(paste("beta sampling tolerance was increased to", tol))
    b<-b[1:nsim,]
  }
  b
}


#' ####################################################################################
#' convert beta vector to ellipsoid object
#' 
#' @param beta beta vector, the coefficients in quadratic form
#' @param d dimension
#' @param ... passed on to 'ellipse_solve_rota'
#' 
#' @export
ellipsoid_from_beta <- function(beta, d, ...){
  elform <- ellipse_form(beta, d)
  chat <- c( elform$c )
  Ahat <- elform$A
  # solve rotation and axes:
  rota <- ellipse_solve_rota(Ahat, ...)
  R <- rota$R
  axes_len <- rota$axes
  M <- R%*%diag(axes_len)
  angles <- NULL
  if(d==2) {
    f <- R %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- rotationMatrix2EulerAngles(R)
  }
  
  # check if we got a valid fit
  valid <- all(!is.infinite(c(axes_len, angles)) & !is.na(c(axes_len, angles)))
  # compile
  
  res <- list(center=chat, A=Ahat, 
              semi_axes=axes_len, 
              rot=R, 
              M=M, 
              rot_angle=angles, 
              valid=valid, dim=d)
  class(res) <- "ellipsoid"
  res
}

#' Create an ellipsoid object
#' 
#' @param R rotation matrix
#' @param semi_axes semi_axes lengths
#' @param center center coordinates
#' 
#' @export
as_ellipsoid <- function(semi_axes=c(1,1,1), 
                         R=diag(0,3), 
                         center=c(0,0,0)) {
  # solve rotation and axes:
  if(ncol(R)!=length(semi_axes))stop("dimension mismatch")
  M <- R%*%diag(semi_axes)
  angles <- NULL
  d <- ncol(R)
  if(d==2) {
    f <- R %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- rotationMatrix2EulerAngles(R)
  }
  
  # check if we got a valid fit
  valid <- all(!is.infinite(c(semi_axes, angles)) & !is.na(c(semi_axes, angles)))
  # compile
  # make the symmetric A matrix by
  # reversing eigenvalue decomposition
  A <- R%*%diag(1/semi_axes^2)%*%t(R)
  #
  res <- list(center=c(center), Ahat=A, 
              semi_axes=semi_axes, 
              rot=R, 
              M=M, 
              rot_angle=angles, 
              valid=valid, dim=d)
  class(res) <- "ellipsoid"
  res
}


#' ############################################################
#' Summarise an ellipsoid
#' 
#' @export
summary.ellipsoid <- function(object, ...){
  x <- object
  print(x)
  ang <- round(x$rot_angle, 3)
  angd <- round(x$rot_angle * 180/pi, 1)
  semi_axes <- x$semi_axes
  semi_axes_rel <- round(semi_axes/semi_axes[1], 3)
  angtxt <- ifelse(x$dim==2, format(ang), paste(c("heading", "attitude", "bank"), format(ang), collapse=" "))
  angtxtd <- ifelse(x$dim==2, format(angd), paste(c("heading", "attitude", "bank"), format(angd), collapse=" "))
  cat("\nEstimates:\n Center:\t \t ", paste0("(", paste0(format(x$center), collapse=", "), ")\n"))
  cat(" Semi-axes lengths (absolute):\t ", paste0(format(semi_axes), collapse=" : "), "\n")
  cat(" Semi-axes lengths (relative):\t ", paste0(format(semi_axes_rel), collapse=" : "), "\n")
  
  cat(" Rotation angles (rad):\t ", angtxt,"\n")
  cat(" Rotation angles (deg):\t ", angtxtd,"\n")
  cat(" Error variance: \t ", x$ols_fit$s2, "\n")
  invisible(NULL)
}


#' Perspective plot of an ellipsoid in 3d
#' 
#' Try to make a pretty perspective plot of an 3D ellipsoid.
#' 
#' @param x ellipsoid-object
#' @param add Add to a scene? 
#' @param theta Camera rotation, see persp.default
#' @param phi Camera rotation, see persp.default
#' @param expand zoom 
#' @param triptych plot from all three main axes?
#' @param pmat optional projection matrix (e.g. from 2d persp)
#'
#' @import graphics
#' @export
persp.ellipsoid <- function(x, add=FALSE, theta=25, phi=30, 
                            expand=.9, triptych=FALSE, 
                            pmat, ...){
  if(x$dim != 3) stop("perspective plot only for 3D ellipsoids.")
  if(triptych){
    a <- c(90, 180, 0)
    b <- c(0, 0,90)
    for(i in 1:3) pmat <- persp.ellipsoid(x, add, theta=a[i], 
                                          phi=b[i], expand, triptych=FALSE, pmat, ...)
  }
  else{
    L <- c(-1,1)*max(x$semi_axes) * 1
    if(!add){
      pmat <- persp(L, L, matrix(NA, 2, 2), 
                    xlim=L, ylim=L, zlim=L,
                    theta=theta, phi=phi, expand=expand , 
                    xlab="X", ylab="Y", zlab="Z", 
                    box=TRUE, ticktype="simple", shade=FALSE)
    }else 
      if(missing(pmat)) 
        stop("Need to provide 'pmat'-projection matrix for adding to persp-plot.")
    #'
    add_ellipsoid2persp(x=x, pmat=pmat, ...)
    #
  }
  invisible(pmat)
}

#' add ellipse to persp plot
#' 
#' @param x ellipsoid
#' @param N refining iterations
#' @param colmap function to generate colors from values
#' @param pmat camera projection matrix
#' @param ... passed on to polygon
#' 
#' @import grDevices
add_ellipsoid2persp <- function(x, N=2, colmap=values2colors, pmat, ...){
  s <- ellipsoid_shape(N, x$semi_axes, x$rot)
  w <- s$vb[4,]
  xc <- s$vb[1,]/w
  yc <- s$vb[2,]/w
  zc <- s$vb[3,]/w
  coords <- cbind(xc,yc,zc)
  Lm <- range(zc)
  zm <- apply(s$it, 2, function(n) mean(zc[n]) )
  cols <- colmap(zm)
  # ah the problem is plotting of behind after front:
  # order according to camera position:
  camera <- c(unlist(trans3d(0, 0, 0, pmat)),100)
  # distance surface
  mini <- apply(s$it, 2, function(n) 
  {d<-t(t(coords[n,]) - camera); n[which.min(diag(d%*%t(d)))]}  )
  o <- order(coords[mini,3])
  for (j in o) {
    idx <- s$it[,j]
    co <- coords[idx,]
    col <- cols[j]
    polygon(trans3d(co[,1], co[,2], co[,3], pmat), col=col, ...)
  }
}


