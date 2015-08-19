#' pcf in 3D. This is my best tool?

library(devtools)
load_all(".")

library(rstrauss)

set.seed(1)
pp <- rstrauss(100, .1, R<-0.03, perfect=TRUE)
#'

g <- pcf_directions(pp, u=c(1, 0, 0), theta=pi/10, correction = "trans", r=seq(0,0.1, length=20))

plot(g$r, g$trans, type="l")
abline(h=1)
abline(v=R, col=4)
