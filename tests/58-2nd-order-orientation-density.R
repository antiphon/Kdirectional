# Test the 2nd order orientation density estimator

library(devtools)
load_all(".")



p <- x <- matrix(runif(300), nc=2)


x <- out <- fry_orientation_density(x, 0.01, 0.1)

plot(x, ylim=c(0,1))
abline(h = 1/pi, col=3)
