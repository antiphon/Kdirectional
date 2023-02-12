# dev and test wavelet transform of intensity

library(devtools)
load_all(".")
source("tests/test-data.R")

x <- pp
#source("~/Sync/Aka/julkaisut/anisotropy-review-overleaf/examples/synthetic-examples/generate-stripe.R")
#x <- stripes

#nn <- 100
#x <-  cbind(runif(nn) -.5, runif(nn) - .5)  * 60
#x[,2] <- round(x[,2])

t0 <- system.time(S <- wavelet_saed(x, k0 = .1, # k0 does not matter
                  scales = scales <- seq(0, 1, l = 50)[-1], 
                  theta =  theta  <- seq(0, pi, l = 50)[-1], shift_res = 10) )
print(t0)

image(y = scales, x = theta, matrix(S, byrow=T, ncol = length(scales)))#, x=theta, y=scales)

