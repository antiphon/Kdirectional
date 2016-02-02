# make test patterns
library(rstrauss)

z <- rstrauss(500, 0.1, 0.03, bbox=cbind(0:1, 0:1)-.5, perfect=T, iter=1e5)
