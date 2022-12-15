
#' bin angles
#' 
#' @param x vector of angles, from -pi to pi
#' @param k number of bins
#' 
#' @export
angles.bin <- function(x, k=25){
  #if(ncol(x)>1) stop("Only for 2D.")
  left <- seq(-pi, pi, length=k)[-k]
  w <- diff(left[1:2])
  i<-apply(cbind(x), 1, function(z) sum(left<z) )
  ang <- left
  v<-table(factor(i, levels=1:(k-1)))
  data.frame(ang=ang, count=c(v))
}

#' Plot a rose of directions for a set of angles
#' 
#' Useful for checking e.g. nnangles- or fry_points- outputs.
#' 
#' @param angles observed angles , -pi to pi
#' @param binned binning of angles. this or angles needed.
#' @param span if >0, span is passed to loess for smoothing and the smoothed values are added to the plot
#' @param k number of bins if binning is done for angles
#' @param ... passed to plot (e.g. main)
#' @param nums plot the values
#' @param ci Plot approximate pointwise 95\% confidence interval
#' 
#' 
#' 
#' @export
angles.flower <- function(angles, binned, span=0, k=25, nums=FALSE, ci=FALSE, ...){
  if(missing(angles) & missing(binned)) stop("angles or binned data needed.")
  if(!missing(angles) & missing(binned)) binned <- angles.bin(angles, k) 
  b <- binned
  w <- diff(b[1:2,1])
  rs <- sqrt(b[,2])
  if(ci) rs <- rs^2/sum(rs^2)
  slice <- function(ang, r, n=nums){
    x <- r*cos(ang+c(0,1)*w)
    y <- r*sin(ang+c(0,1)*w)
    polygon(xx<-c(0, x, 0), yy<-c(0, y, 0), border="white", 
            col=rgb(0.4,0.4,0.4,0.5))
    mx<-mean(c(x))
    my<-mean(c(y))
    if(n)text(mx, my, round(r^2), cex=0.9)
    c(mx,my)
  }
  
  l <- max(rs)
  plot(NA, xlim=c(-l,l), ylim=c(-l,l), xlab="", ylab="", 
       asp=1, ..., xaxt="n", yaxt="n")
  xy <- apply(cbind(b[,1], rs), 1, function(ab)slice(ab[1], r=ab[2]))
  points(0,0,pch=19, cex=0.2)
  
  # add smooth ring
  if(span>0) {
    e <- angles.smooth(b, span=span)
    angs <- e[,1]+w/2; 
    angs <- c(angs, angs[1]) # close the loop
    r <- sqrt(e[,2]); r <- c(r, r[1])
    x <- r*cos(angs)
    y <- r*sin(angs)
    lines(x, y, lwd=4, col="goldenrod")
  }
  # if ci
  if(ci) {
    p <- w/(2*pi)
    s <- sqrt(p*(1-p)/length(angles))
    rl <- p - 1.96*s
    ru <- p + 1.96*s
    symbols(c(0,0,0),c(0,0,0), circles=c(rl, p, ru), fg=c(4,3,4), add=T, inches=F)  }
  invisible(rs)
}



#' smoothen an histogram of angles
#' 
#' @param binned binning of angles
#' @param ... passed to loess
#' @export 

angles.smooth <- function(binned, ...){
  ang <- binned[,1]
  co <- binned[,2]
  n <- nrow(binned)
  x <- c(ang-2*pi, ang, 2*pi+ang)
  v <- sqrt(c(co,co,co))
  f <- predict(loess(v~x, ...))
  b <- pmax(f[1+1:n], 0)
  data.frame(ang=ang, count=b)
}

#' Cumulative distribution function and pointwise CI
#' 
#' @param angles angles vector, from 2D pattern
#' @param n grid resolution
#' @param ... ignored
#' @export
angles.cdf <- function(angles, n=50, ...){
  ang <- angles
  x <- seq(0, 2*pi, length=n)
  if(min(ang)<0) ang <- ang + pi # in case [-pi,pi]
  n <- length(ang)
  f <- ecdf(ang)
  cdf <- f(x)
  cdf0 <- x/(2*pi)
  s <- sqrt( cdf0*(1-cdf0)/n )
  u <- cdf0 + 1.96*s
  l <- cdf0 - 1.96*s
  data.frame(angle=x, cdf=cdf, CI5=l, CI95=u)
}

