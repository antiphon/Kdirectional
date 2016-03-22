# test intersercion

library(devtools)
load_all(".")

# gotta be able to compute intersections really fast

n <- 100^2

x <- bbpoly_default()

shifts <- matrix( runif(n*3, -.2,.2), nc=3)

t0 <- system.time(v <-  apply(shifts, 1, function(s) {
  y <- x
  y$vertices <- t(t(y$vertices)+s)
  w <- bbpoly_intersection(x, y)
  bbpoly_volume(w)
}  )  )
print(t0)
