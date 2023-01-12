# #' Kolmogorov-Smirnov test using ecdf's
# #' 
# #' @param x empirical cdf function, see details
# #' @param y theoretical cdf function
# #' @param n Number of 'datapoints'
# #' 
# #' @useDynLib Kdirectional
# #' @export
# 
# ks.test2 <- function(x, y, n) {
#   if(ncol(x)!=2)stop("x should be 2-col matrix with first col the domain and second column the cdf values.")
#   if(missing(n))stop("Need the number of datapoints from which x is computed.")
#   d <- max(  abs(x[,2]-y(x[,1])) )
#   p <- c_p_of_KStest(n, d)
#   c(statistic=d, p.value=p)
# }
## Obsolete.