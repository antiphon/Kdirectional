#' Fit second order intensity ellipses or ellipsoids
#' 
#' Difference to v1: Don't use the pseudo-Fry but true Fry.
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
#' @import ellipsoid
#' @export

fry_ellipsoids2 <- function(x, nvec=1:5, r_adjust=1, nangles, eps=0, 
                           cylindric=FALSE, double=FALSE, 
                           border=TRUE, 
                           origin=TRUE, 
                           verbose=FALSE,
                           keep_data=FALSE,
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
  if(missing(nangles)) {
    n <- nrow(x$x)
    # use poisson approximation and requiremnt P(fry's in sector) ~= n(x)/const
    nangles <- if(dim==2) round( n/6 )
    else round( n * pi / 20 )
    if(dim == 3) nangles <- (1:5)[ which.min(abs(nangles- 20*4^(1:5)/2+2) )]
  }
  if(dim==2){
    angles <- seq(0, 2*pi, length=nangles+1)[-1]
    delta <- diff(angles[1:2])/2
    angles <- angles -  delta
    grid_unit <- cbind(cos(angles), sin(angles))
  }
  else{
    grid_unit <- directions <- ll2xyz( sphere.grid(N=nangles, ico=TRUE) )
    nng <-  apply(as.matrix(dist(directions)), 1, function(v) sort(v)[2])
    mdi <- mean(nng) # mean distance between directions in the grid
    # angle between two unit directions with distance mdi
    delta <- asin(mdi) * 0.57 # >.5 factor to get the full coverage, .5 too small. 1/sqrt(3) should do it
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
  # conical
  if(!cylindric){
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      s <- data_units%*%u
                      a <- acos(s)
                      a <- pmin(pi-a, a)
                      a< delta & s>=0
                    })
    if(double) {
      inside_antipode <- apply(-grid_unit, 1, 
                      function(u) {
                        s <- data_units%*%u
                        a <- acos(s)
                        a <- pmin(pi-a, a)
                        a< delta & s>=0
                      })
    }
  }
  else { # cylindrical
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      s <- c(data%*%u)
                      p <- data - t(sapply(s, "*", u))
                      d <- sqrt(rowSums(p^2))
                      d < delta & s >= 0
                    })
  }
  ## get the ranges of each fry point, sorted
  Kd     <- apply(inside, 2, function(i) cbind(idx=which(i)[order(data_r[i])] ) )
  if(double){
    Kd_antipode  <- apply(inside_antipode, 2, function(i) cbind(idx=which(i)[order(data_r[i])] ) )
  }
  #
  ns <- sapply(Kd, length) # counts per sector
  #
  cat2("Sector inclusions resolved, range of point count varies ", min(ns), "-", max(ns),"\n")
  if(max(nvec) > min(ns)) {
    txt <- paste("Requested neighbour count not available for all directions (min is", min(ns),")\n")
    cat2(txt)
  }
  
#   ## directional isocurves, the ranges jumps in each direction
#   iso_ranges <- sapply(Kd, function(k) {
#     ok <- nvec[nvec <= length(k)]
#     e <- k[ok, 1]  # need nvec out of this
#     n <- length(e)
#     if(n < length(nvec)) e <- c(e, rep(NA, length(nvec)-n))
#     e
#   } )
#  iso_ranges <- rbind(iso_ranges)
  # compile the idx table
  
  iso_idx <- sapply(Kd, function(k) {
    ok <- nvec[ nvec <= length(k)]
    ko <- k[ok,1]
    n <- length(ko)
    if(n < length(nvec)) ko <- c(ko, rep(NA, length(nvec)-n))
    ko
  } ) 
  iso_idx <- rbind(iso_idx)
  
  if(double){
    iso_idxa <- sapply(Kd_antipode, function(k) {
      ok <- nvec[ nvec <= length(k)]
      ko <- k[ok,1]
      n <- length(ko)
      if(n < length(nvec)) ko <- c(ko, rep(NA, length(nvec)-n))
      ko
    } ) 
    iso_idxa <- rbind(iso_idxa)
  }
  
  # Fit ellipses:
  
  # check counts
  counts <- apply(iso_idx, 2, function(v) sum(!is.na(v)) )
  cat2("Range vectors per direction jumps computed, count varies ", min(counts), "-",max(counts),"\n")
  perc <- apply(iso_idx, 1, function(v)sum(!is.na(v)))
  cat2("Data points per contour varies ", min(perc), "-",max(perc),"\n")
  
  good_k <- perc > (3+dim)
  if(any(perc< (3+dim))) cat2(sum(!good_k),"contours not computed, not enough data.\n")
  
#   els <- apply(iso_idx[good_k, ], 1, function(idxs){
#     p <- data[idxs,]
#     p <- p[!is.na(idxs),]
#     if(double) p <- rbind(p, -p)
#     el <- ellipsoid_OLS(p, origin=origin) # in 'ellipsoid' package 
#     if(keep_data) el$data <- p
#     el
#   } )
  
  els <- lapply(which(good_k), function(k) {
    p <- data[iso_idx[k,],]
    if(double) p <- rbind(p, data[iso_idxa[k,],])
    el <- ellipsoid_OLS(p, origin=origin) # in 'ellipsoid' package 
    if(keep_data) el$data <- p
    el
  })
  nvec <- nvec[good_k]
  # done
  res <- list(ellipsoids=els, dim=dim, param=list(cylindric=cylindric, delta=delta, eps=eps),
              jump_idx=iso_idx, 
              fry=data,
              n=nvec, grid_unit=grid_unit, v2=TRUE)
  #
  class(res) <- "fryellipsoids"
  res
}





