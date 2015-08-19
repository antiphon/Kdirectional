# tes tprofile with one-vs-many

#' test
library(devtools)
load_all(".")

pp <- list(x=matrix(runif(200), ncol=2), bbox=cbind(c(0,1),c(0,1)))

anisotropy_profile(pp, i=1, target=1, cvec=(2:11)/10, epsilon=pi/8)

#'

