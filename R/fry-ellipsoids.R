#' Fit second order intensity ellipses or ellipsoids
#' 
#' @param x point pattern
#' @param nvec which jumps to consider
#' @param r_adjust consider only fry points with length r_adjust*max(window side) 
#' @param eps add to the sector angle. Results in smoothing if >0 (default=0).
#' @param nangles If 3d, number of refinements (def=2), in 2d def=20 sectors.
#' @param cylindric Use a cylinder instead of a sector
#' @param border Should border effects be cancelled by minus sampling, if TRUE, "double" ignored.
#' @param double Should we use both directions of a pair
#' @param origin Constrain the ellipses' centers to origo.
#' @param keep_data Keep the data for each ellipsoid, inside each ellipsoids? Helps with averaging.
#' @param data The vectors to which fit ellipses. Default is missing, and we compute Fry-points.
#'
#' @import sphere
#' @export

fry_ellipsoids <- function(x, nvec=1:5, r_adjust=1, nangles, eps=0, 
                            cylindric=FALSE, double=FALSE, border=TRUE, origin=TRUE, verbose=FALSE,
                            keep_data=TRUE,
                            data, ...) {
  x <- check_pp(x)
  dim <- ncol(x$x)
  #
  cat2 <- if(verbose) cat else function(...)NULL
  #
  ## Border correction
  if(missing(data)){
    fry <- fry_points(x, double=double, border=border)
    data_r <- fry$fry_r
    data_units <- fry$fry_units
  }
  else{ # data given
    data_r <- sqrt(rowSums(data^2))
    data_units <- data/data_r
  }
  # drop long range for easier computations
  r_ok <- data_r < r_adjust * max(bbox_sideLengths(x$bbox))
  data_r <- data_r[r_ok]
  data_units <- data_units[r_ok,]
  data <- data_units * data_r
  
  cat2("Fry points gathered, n=", nrow(data),"\n")
  #
  ##### The directions -grid
  if(missing(nangles)) nangles <- ifelse(dim==3, 2, 20)
  if(dim==2){
    angles <- seq(0, 2*pi, length=nangles+1)[-1]
    delta <- diff(angles[1:2])/2
    angles <- angles -  delta
    grid_unit <- cbind(cos(angles), sin(angles))
  }
  else{
    grid_unit <- directions <- ll2xyz( sphere.grid(N=nangles, ico=TRUE) )
    nng <-  apply(as.matrix(dist(directions)),1, function(v) sort(v)[2])
    mdi <- mean(nng)
    delta <- asin(mdi/2)
  }
  if(cylindric){
    # go for the width of the cylinder grown up to nearest neighbour of origin
    # in any direction
    m <- sort(data_r)[1]
    delta <- sin(delta) * m 
  }
  cat2("Direction grid generated, n=", nrow(grid_unit),"\n")
  #
  #### sector width. eps >0 results in overlaps
  delta <- eps + delta
  #
  #### Compute the inclusions per fry point in each cone:
  if(!cylindric){
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      a <- acos(data_units%*%u)
                      a <- pmin(pi-a, a)
                      a < delta
                    })
  }
  else {
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      s <- c(data%*%u)
                      p <- data - t(sapply(s, "*", u))
                      d <- sqrt(rowSums(p^2))
                      d < delta
                    })
  }
    ## get the ranges of each fry point, sorted
  Kd     <- apply(inside, 2, function(i) cbind(r=sort(data_r[i])) )
  #
  ns <- sapply(Kd, length)
  cat2("Sector inclusions resolved, range of point count varies ", min(ns), "-", max(ns),"\n")
  if(max(nvec) > min(ns)) {
    txt <- paste("Requested neighbour count not available for all directions (min is", min(ns),")\n")
    cat2(txt)
  }
  
  ## directional isocurves, the ranges jumps in each direction
  iso_ranges <- sapply(Kd, function(k) {
      ok <- nvec[nvec <= length(k)]
      e <- k[ok, 1]  # need nvec aout of this
      n <- length(e)
      if(n < length(nvec)) e <- c(e, rep(NA, length(nvec)-n))
      e
    } ) 
  # Fit ellipses:
  iso_ranges <- rbind(iso_ranges)
  counts <- apply(iso_ranges,2,function(v)sum(!is.na(v)))
  cat2("Range vectors per direction jumps computed, count varies ", min(counts), "-",max(counts),"\n")
  perc <- apply(iso_ranges, 1, function(v)sum(!is.na(v)))
  cat2("Data points per contour varies ", min(perc), "-",max(perc),"\n")
  
  good_k <- perc > (3+dim)
  if(any(perc< (3+dim))) cat2(sum(!good_k),"contours not computed, not enough data.\n")
  
  els <- apply(iso_ranges[good_k, ], 1, function(r){
    p <- r * grid_unit
    p <- p[!is.na(r),]
    el <- ellipsoid_OLS(p, origin=origin) # in 'sphere' package 
    if(keep_data) el$data <- p
    el
  } )
  nvec <- nvec[good_k]
  # done
  res <- list(ellipsoids=els, dim=dim, param=list(cylindric=cylindric, delta=delta, eps=eps),
              jumps=iso_ranges, 
              fry=data,
              n=nvec, grid_unit=grid_unit)
  #
  class(res) <- "fryellipsoids"
  res
}

#############################################################
#' plot method for fryellipsoids
#' 
#' @export
#' @exportMethod plot


plot.fryellipsoids <- function(x, ellipsoids=TRUE, used_points=TRUE, sectors=NULL, xlim=NULL, ylim=NULL, ...) {
  if(x$dim==2){
    plot(x$fry, asp=1, xlab="x", ylab="y", main="Fry points, estim points(x) and fitted ellipses", xlim=xlim, ylim=ylim)
    for(i in 1:length(x$ellipsoids)){
      if(used_points) points(x$grid_unit * x$jumps[i,], pch=20, col=i)
      if(ellipsoids) plot(x$ellipsoids[[i]], col=i, ...)
    }
    abline(h=0, col="gray50")
    abline(v=0, col="gray50")
    
    # 
    if(any(sectors)){
      delta <- x$param$delta
      r <- max(x$fry)
      for(i in sectors){
      u <- x$grid_unit[i,]
      if(x$param$cylindric){
        ut <- u[2:1]*c(-1,1)
        uv <- sapply(seq(0, r, length=2), "*", u )
        lines(rbind(delta*ut, -delta*ut))
        lines(t(delta*ut+uv))
        lines(t(-delta*ut+uv))  
      }
      else{
        m <- cbind(c(cos(delta), sin(delta)),c(-sin(delta), cos(delta)))
        u1 <- c(m%*%u)
        u2 <- c(solve(m)%*%u)
        lines(rbind(c(0,0), u1))
        lines(rbind(c(0,0), u2))
      }
      }
    }
    
  }
  else warning("Plot not implemented in 3d")
}

#############################################################
#' Print method for fryellipsoids
#' 
#' @exportMethod print
#' @export
print.fryellipsoids <- function(x, ...) {
  cat("Ellipsoid fit to fry points.\n")
}









