#' Test K-S distribution

library(devtools)
load_all(".")

X<-runif(1000)
v <- ecdf(X)

r <- seq(0,1, length=50)
x <- cbind(r, v(r))
y <- function(r) r

D <- # TODO: Problem is that what is n?

p <- pKolmogorov()


plot(v)
curve(y, add=T, col=4)
print(p)
