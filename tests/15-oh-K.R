#' Ohser-Stoyan

library(devtools)
load_all(".")

library(rstrauss)
x <- rstrauss(gamma=0.01, range=0.07, bbox=cbind(0:1, 0:1), perfect=T)
x0 <- rstrauss(gamma=1, range=0.07, bbox=cbind(0:1, 0:1), perfect=T)
est <- osK(x)
est0 <- osK(x0)


par(mfrow=c(2,1))
plot(est0)
plot(est)
