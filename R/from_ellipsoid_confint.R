#################################################################
#' Confidence interval for ellipsoid parameters using simulation
#' 
#' Assuming normality of errors.
#' 
#' @param x Fitted ellipsoid
#' @param fun Optional contrast function. See ellipse_contrast_2d for an example.
#' @param probs The quantiles for the CI. Default is c(0.025, 0.975)
#' @param tol tolerance for absolute  deviation in the ||beta||=1 constraint.
#' @param Asolve.eps Inflate the A matrix diagonals by a factor of eps, to avoid numerical problems with (near) singular parameter sets.
#' @param ... Passed on to sample_ellipse_beta
#' 
#' @import mvtnorm
#' @export

confint.ellipsoid <- function(x, fun, nsim=1000, probs=c(0.025, 0.975),  
                              Asolve.eps=0.01, ...){
  S <- x$ols_fit$varcov
  if(is.null(S)) stop("Ellipsoid object does not contain varcov matrix. (ALS not supported)")
  m <- x$ols_fit$beta_est
  f <- qt(probs, x$ndata-2)
  s <- sqrt(diag(S))
  L <- m + s*f[1]
  U <- m + s*f[2]
  names(m) <- paste0("beta", 1:length(m))
  pval <- round(2*pt(-abs(m/s) , x$ndata - 2) , 5)
  df_orig <- data.frame(mean=m, median=m, sd=s, L,U, p=pval)
  # simulate in the quadratic form
  betas <- sample_ellipse_beta(x, nsim, ...)
  #
  efs <- apply(betas, 1, function(b) ellipse_form(b, x$dim) )
  rots <- lapply(efs, function(e) ellipse_solve_rota(e$A, Asolve.eps) )
  axes <- t(sapply(rots, function(b) b$axes ))
  #
  # summary for one vector
  sum.p <- function(v, h0=0) {
    # drop bad values
    v <- v[!is.infinite(v) & !is.na(v)]
    z <- (h0-mean(v))/sd(v)
    p <- round(2 * pt(-abs(z), x$ndata-2), 5)
    c(mean=mean(v), median=median(v), sd=sd(v), 
      CI=quantile(v, probs=probs[1]), 
      CI=quantile(v, probs=probs[2]), 
      p=p)
  }
  # basics for axes
  df_axes <- t(apply(axes, 2, sum.p))
  #
  df_fun <- df_ex <- df_diff <- NULL
  if(x$dim==2){
    # eccentricity
    df_ex <- sum.p(sqrt(1-(axes[,1]/axes[,2])^2 )) # eccentricity
    df_ex[6] <- NA # these p-values are not correct
    # custom function
    if(missing(fun))
      df_diff <- sum.p(ellipse_contrast_2d(efs))
  }
  
  if(!missing(fun)) df_fun <- sum.p(fun(efs))
  
  names(df_orig) <- colnames(df_axes)
  rbind(beta=df_orig, axes=df_axes, eccentricity=df_ex, isotropy=df_diff, custom=df_fun)
}




#################################################################
#' Confidence interval for ellipsoid parameters using simulation
#' 
#' Assuming normality of errors.
#' 
#' @param x Fitted ellipsoid
#' @param fun Optional contrast function. See ellipse_contrast_2d for an example.
#' @param probs The quantiles for the CI. Default is c(0.025, 0.975)
#' @param tol tolerance for absolute  deviation in the ||beta||=1 constraint.
#' 
#' @import mvtnorm
# #' @export

confint.ellipsoid_old <- function(x, fun, nsim=1000, probs=c(0.025, 0.975),  ...){
  S <- x$ols_fit$varcov
  if(is.null(S)) stop("Ellipsoid object does not contain varcov matrix. (ALS not supported)")
  m <- x$ols_fit$beta_est
  f <- qt(0.975, x$ndata-2)
  s <- sqrt(diag(S))
  L <- m - s*f
  U <- m + s*f
  names(m) <- paste0("beta", 1:length(m))
  pval <- round(2*pt(-abs(m/s) , x$ndata - 2) , 5)
  df_orig <- data.frame(mean=m, median=m, sd=s, L,U, p=pval)
  # simulate in the quadratic form
  betas <- sample_ellipse_beta(x, nsim, ...)
  #
  efs <- apply(betas, 1, function(b) ellipse_form(b, x$dim) )
  rots <- lapply(efs, function(e) ellipse_solve_rota(e$A) )
  axes <- t(sapply(rots, function(b) b$axes ))
  #
  # summary for one vector
  e <- function(v, h0=0) {
    # drop bad values
    v <- v[!is.infinite(v) & !is.na(v)]
    z <- (h0-mean(v))/sd(v)
    p <- round(2*pt(-abs(z), x$ndata-2), 5)
    c(mean=mean(v), median=median(v), sd=sd(v), 
      CI=quantile(v, probs=probs[1]), 
      CI=quantile(v,probs=probs[2]), 
      p=p)
  }
  # basics for axes
  df_axes <- t(apply(axes, 2, e))
  #
  df_fun <- df_ex <- df_diff <- NULL
  if(x$dim==2){
    # eccentricity
    df_ex <- e(sqrt(1-(axes[,1]/axes[,2])^2 )) # eccentricity
    df_ex[6] <- NA # these p-values are not correct
    # custom function
    if(missing(fun))
      df_diff <- e(ellipse_contrast_2d(efs))
  }
  
  if(!missing(fun)) df_fun <- e(fun(efs))
  
  names(df_orig) <- colnames(df_axes)
  rbind(beta=df_orig, axes=df_axes, eccentricity=df_ex, isotropy=df_diff, custom=df_fun)
}


