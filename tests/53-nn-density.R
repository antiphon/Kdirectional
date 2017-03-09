# test nn angle density, IPSS p253

library(devtools)
load_all(".")

if(!exists("x")) x <- matrix(runif(200), nc=2)



par(mfrow=c(1,3))

# 0-2pi
f <- nnangle.density(x, antipodal = F)
plot(f, type ="l", ylim=yl<-c(0,.6))
abline(h=1/(2*pi))

# wrapping working?
plot(f, type ="l", xlim=c(-pi,3*pi), ylim=yl)
 lines(f$angle+pi*2, f$density, col=2)
 lines(f$angle-pi*2, f$density, col=2)

# 0-pi
fa <- nnangle.density(x, antipodal = T)
plot(fa, type ="l", ylim=yl)
abline(h=1/(pi))

# data input:
D <- nnangle.density(x, justData = T , antipodal=T)
fb <- nnangle.density(data = D, antipodal = T)

lines(fb, col=2, lty=3)
