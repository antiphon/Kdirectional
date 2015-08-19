#' Test the pairwise computering
library(devtools)
load_all(".")

d<-2
x <- matrix(runif(10*d), ncol=d)

t0 <- system.time( p <- pairwise_vectors(x, asMatrix=T) )

#' C
t1 <- system.time( pd <- as.matrix(dist(cbind(x))) )

print(all.equal(pd, p$dist, check.attributes = F))

#' the from-to formulation
t2 <- system.time( p3 <- pairwise_vectors(x, from=1:nrow(x)) )




#' check by plotting the arrows

# plot(x, asp=1)
# i <- 2
# for(j in setdiff(1:nrow(x),i)){
#   a <- p$angle[i,j]
#   aa <- p$angle[j,i]
#   d <- p$d[i,j]
#   dd <- p$d[j,i]
#   uf <- function(a,d)d*c(cos(a), sin(a))
#   u <- uf(a,d)
#   uu <- uf(aa,dd)
#   arrows(x[i,1],x[i,2],x[i,1]+u[1], x[i,2]+u[2], col="green", length = .2)
#   arrows(x[j,1],x[j,2],x[j,1]+uu[1], x[j,2]+uu[2], col="red", lty=4, length=.2)
# }