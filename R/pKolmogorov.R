#' Limiting distribution of the Kolmogorov-Smirnov test statistic
#' @details
#' Will return approximately lim P(sqrt(n)Dn < x) .
#' @export

pKolmogorov <- function(x, upto=999){
  i <- 1:upto
  sapply(x, function(x) sqrt(2*pi)/x * sum( exp(-(2*i-1)^2*pi^2/(8*x^2))  ))
}

