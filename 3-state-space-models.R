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

sink('ssm.jags')
cat("
model {
# Priors
N.est[1] ~ dunif(0, 500) # intital population size
mean.lambda ~ dunif(0, 5)
sigma.proc ~ dunif(0, 10)
sigma2.proc<-sigma.proc*sigma.proc
tau.proc<-pow(sigma.proc, -2)

sigma.obs ~ dunif(0,100)
sigma2.obs<-sigma.obs*sigma.obs
tau.obs<-pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)) {
  lambda[t] ~ dnorm(mean.lambda, tau.proc)
  N.est[t+1]<- N.est[t]*lambda[t]
} #t
# Observation process
for (t in 1:T) {
  y[t] ~ dnorm(N.est[t], tau.obs)
} #t
}
", fill=TRUE)
sink()

jags.data<-list(y=y, T=n.years)

inits<-function() list(sigma.proc=runif(1, 0, 1), mean.lambda=runif(1, 0.1, 2),
sigma.obs=runif(1,0,10), N.est=c(runif(1, 20, 40), rep(NA, n.years-1)))

parameters<-c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)








