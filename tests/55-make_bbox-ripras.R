#' Test RipRas bounding box expansion in bbox_make
#' 

library(devtools)
load_all(".")

bb0 <- cbind(c(0, 500), c(0, 500), c(2800, 4500))

set.seed(1)
x <- apply(bb0, 2, function(a) runif(200, a[1], a[2]))
x2 <- x[,-3]

print(bbox_make(x, expand = TRUE))
print(bbox_make(x, expand = FALSE))

# compare to spatstat

library(spatstat)

print( W0 <-as.owin(c(bbox_make(x2, F))) )
print( W<-ripras(x2, shape="rectangle") )
print( w<-as.owin(c(bbox_make(x2, expand = TRUE) )))

print(data.frame(noexp=area(W0), ripras=area(W), oma=area(w)))

plot(w, border=4)
points(x2)
plot(W, border=3, add=T)
#plot(W0, border=1, add=T)
plot(as.owin(c(bb0[,-3])), add=T, border=9)
