# debug henrike 3d - markcor 
library(devtools)
load_all(".")
#library(Kdirectional)
# load data
load("~/Downloads/mcf_data.Rdata")
xb0 <- data$xb0
m0 <- as.numeric( data$m0 )
scale <- 208.7311
# call function
mcor30 <- markcorr_anisotropic(xb0, marks=c(m0), r = seq(0, 70/scale, l=50))
# plot
plot(mcor30$r*scale, mcor30$mcor[,1], type="l", ylim=c(0,2), main="MCF: Sample 1", col=2, lwd=3, ylab="MCF(r)", xlab="r")
for(i in 2:3) lines(mcor30$r*scale, mcor30$mcor[,i], col=i+1, lwd=3)
abline(h=1,lty=2)
legend("topright", legend = c('(1,0,0)','(0,1,0)','(0,0,1)'), cex=1, col=c(2,3,4), lwd=2, bty="n")
