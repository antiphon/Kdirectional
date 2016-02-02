# test epanechnicov kernel integral calculation

library(devtools)
load_all(".")


ss <- seq(0,1, length=20)
x <- as.matrix(expand.grid(ss, ss, ss))
bbox <- cbind(0:1,0:1, 0:1)

t0 <- system.time( v <- epa_integral(x, bbox, bw<-0.2) )

# check it integrates to 1
print(c(epa_integral(cbind(.5,.5), bbox[,-3], 0.1) , 
        epa_integral(cbind(.5,.5,.5), bbox, 0.1)))

#print(rbind(t0,t1))

# library(rgl)
# library(scales)
# co <- gray(rev(i<-scales::rescale(v)))
# 
# inc <-  abs(x[,3]-.5) < 0.1 
# 
# plot3d(x, col=co, aspect=F, alpha=rev(i))
