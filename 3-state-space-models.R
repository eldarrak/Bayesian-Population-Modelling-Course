# state-space models


n.years=25
N1 <-30
mean.lambda<-1.02

sigma2.lambda<-0.02 # process varitation

sigma2.y<-20 #observation error

y<-rep(NA, n.years)
N<-rep(NA, n.years)
N[1]<-N1
lambda<-rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))

for (t in 1:(n.years-1)) {
   N[t+1]<-N[t]*lambda[t]
}

for (t in 1:n.years) {
  y[t]<-rnorm(1, N[t], sqrt(sigma2.y))
}

y



