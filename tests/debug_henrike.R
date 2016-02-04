# Solved. Was division with area even though we had the translation correction.


library(devtools)
load_all(".")
x <- read.csv("~/Dropbox/work/students/henrike/data/cor_kbranchp.txt", header=F, skip=1, sep=" ")[,-1]
x <- as.matrix(x)
x1 <- x/(s1<-408.7311)
x2 <- x/(s2<-1000)
b1 <- bbox_make(x1)
b2 <- bbox_make(x2)
#z <- epa_integral(x1, as.matrix(bbox_make(x1, F)), bw=.7, n=31)

eps <- pi/2
ll1 <- intensity_at_points(x1, b=31, bw=0.7)
ll2 <- intensity_at_points(x2, b=31, bw=0.5)
#ll1 <- intensity_at_points(x1, b=2, bw=0.5)
#ll2 <- intensity_at_points(x2, b=2, bw=0.5)

# ll1 <- 1/rep(bbox_volume(b1)/nrow(x1), nrow(x1))
# ll2 <- 1/rep(bbox_volume(b2)/nrow(x2), nrow(x2))

p1 <- list(x=x1, bbox=b1)
p2 <- list(x=x2, bbox=b2)

i1 <- Kest_anin(p1, epsilon = eps, 
                renormalise=T, 
                lambda=ll1)
i2 <- Kest_anin(p2, epsilon = eps, 
                renormalise=T,
                lambda=ll2)

par(mfrow=c(2,1))
plot(i1, r_scale=s1)
plot(i2, r_scale=s2)
# 
# 
# #' poisson perkele
# x <- matrix(runif(100*3),ncol=3)
# y <- x/3
# x <- list(x=x, bbox=cbind(0:1,0:1,0:1))
# y<- list(x=y, bbox=cbind(0:1,0:1,0:1)/3)
# #a <- Kest_anin(x, lambda = rep(100, 100), u=u<-cbind(1,0,0), border=0)
# #b <- Kest_anin(y, u=u, lambda= rep(100, 100)/bbox_volume(bbox_make(y,F)), border=0)
# a <- Kest_anin(x, lambda = l<-1, u=u<-cbind(1,0,0), border=1)
# b <- Kest_anin(y, u=u, lambda= l, border=1)
# plot(a)
# lines(b$r, b$`(1,0,0)`)
# plot(b)
# 
