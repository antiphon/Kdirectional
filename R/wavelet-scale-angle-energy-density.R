#' Scale-angle-energy-density
#' 
#' 2D at the moment. Morlet anisotropic wavelet
#' 
#' @param x point pattern, list of $x and $bbox
#' @param theta angles
#' @param scales scales
#' @param shift_res resolution of the shift-grid
#' @param ... ignored
#' @export
wavelet_saed <- function(x, theta = seq(0, pi, l = 30), 
                            scales, 
                            shift_res = 20,
                            ...) {
  x <- check_pp(x)
  loc <- t( x$x )
  bb <- x$bbox
  sl <- apply(bb, 2, diff)
  
  if(length(shift_res)<2) shift_res <- shift_res * c(1, sl[2]/sl[1])
  if(missing(scales)) scales <- seq(0, max(sl), length.out = 20)
  # wavelet
  # mother wavelet, fixed
  D <- 0.1
  k0 <- cbind(0, 5.5)
  con <- sqrt(D/pi) 
  
  wave <- function(a, b, ang) {
    R <- rotationMatrix3(az = -ang)[-3,-3]
    zx <- R%*%(loc - b) /a  
    #e <- exp(1i * k0 %*% zx - 0.5 * colSums((zx)^2) )
    #psia2 <- abs(e * con)^2
    d <- k0%*%zx
    psia2 <- con^2 * exp(-(zx[1,]^2*D^2 + zx[2,]^2)) * ( sin(d)^2 + cos(d)^2 )
    sum(psia2)/a
  }
  ## 
  bgrid <- as.matrix(expand.grid( bx <-  seq(bb[1,1], bb[2,1], l = shift_res[1]),
                        by <-  seq(bb[1,2], bb[2,2], l = shift_res[2])) )
  db <- diff(bx[1:2])*diff(by[1:2])
  S <- sapply( theta, function(ang) {
    R <- rotationMatrix3(az = ang)[-3,-3]
    v <- sapply(scales, function(a){
      sa <- apply(bgrid, 1, function(b){
        zx <- R%*%(loc - b) /a  
        d <- k0 %*% zx
        psia2 <- con^2 * exp(-(zx[1,]^2*D^2 + zx[2,]^2)) * ( sin(d)^2 + cos(d)^2 )
        sum(psia2)/a
      })
      sum(sa^2) * db
    })
    v
  })
  S
}


