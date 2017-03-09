# test the redenbach erosion terms
# | W\ominus S(u,e,r)| 
# with double cone S with direction u, half-angle e and range r.



W <- cbind(0:1, c(0,2))

u <- c(1,0)
eps <- pi/10
r <- 0.1

plot(W, asp=1)
rect(0, 0, 1, 2)

#xy <- #