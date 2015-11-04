#' the sector for sphere grid

library(devtools)
load_all(".")

library(rgl)

g <- ll2xyz( sphere.grid(2, ico = TRUE) ) 

plot3d(g, aspect=F
       )


library(spatgraphs)

G <- spatgraph(p<-list(x=g[,1], y=g[,2], z=g[,3]))

# angles between i and nn(i)

n<-sapply(1:G$N, function(i) acos(t(g[i,])%*%g[G$edges[[i]][1],]) )

