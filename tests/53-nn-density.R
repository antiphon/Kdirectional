# test nn angle density, IPSS p253

library(devtools)
load_all(".")

if(!exists("x")) x <- matrix(runif(200), nc=2)



par(mfrow=c(1,2))

f <- nnangle.density(x)
plot(f, type ="l", ylim=c(0,.5))
abline(h=1/pi)



# toroidal ok ?
#plot(f, type ="l", xlim=c(-pi,3*pi), ylim=c(0, .3))
# lines(f$angle+pi*2, f$density, col=2)
# lines(f$angle-pi*2, f$density, col=2)

fa <- nnangle.density(x, antipodal = F)
plot(fa, type ="l", ylim=c(0,.5))
abline(h=1/(2*pi))

# data input
D <- nnangle.density(x, justData = T , antipodal=F)

fb <- nnangle.density(data = D, antipodal = F)

lines(fb, col=2, lty=3)
