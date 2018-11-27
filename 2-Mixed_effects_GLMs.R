# mixed effects GLMs

data.fn<-function(nsite=5, nyear = 40, alpha = 4.18456, beta1 =1.90672, beta2 = 0.10852, beta3 = -1.17121, sd.site=0.5, sd.year=0.2 ) {

# Generate data structure to hold counts and log(lambda)
C<-log.expected.count<-array(NA, dim=c(nyear, nsite))

year<-1:nyear
yr<-(year-mean(year))/sd(year)
site<-1:nsite

alpha.site<- alpha+ rnorm(n=nsite, mean=0, sd=sd.site)
eps.year<- rnorm(n=nyear, mean=0, sd=sd.year)

for (j in 1:nsite) {
  log.expected.count[,j]<-alpha.site[j] + 
  beta1*yr+ 
  beta2*yr^2+
  beta3*yr^3 + 
  eps.year
  
  expected.count<-exp(log.expected.count)
  
  C[,j] <-rpois(n=nyear, lambda=expected.count)
 
} #j
 return(list(nsite=nsite, nyear=nyear, 
             alpha.site=alpha.site,
             beta1=beta1,
             beta2=beta2,
             beta3=beta3,
             year=year,
             yr=yr,
             sd.site=sd.site,
             sd.year=sd.year,
             expected.count=expected.count,
             C=C))
}

data<-data.fn(nsite=10, nyear=40, sd.site=0.3, sd.year=0.2)

data

sink('GLMM_Poisson.jags')
cat("
model {
# priors
for (j in 1:nsite) {
   alpha[j]~dnorm(mu, tau.alpha)
} #j
mu~dnorm(0, 0.01)
tau.alpha<-1/(sd.alpha*sd.alpha)
sd.alpha~dunif(0, 1)

for (p in 1:3){
   beta[p] ~ dnorm(0, 0.01)
   } #p

tau.year <- 1 / (sd.year*sd.year)
sd.year ~ dunif(0, 1)				# Hyperparameter 3

#Likelihood
for (i in 1:nyear) {
  eps[i] ~ dnorm(0, tau.year)
  for (j in 1:nsite) {
     C[i,j] ~ dpois(lambda[i,j]) # random part
     lambda[i,j] <-exp(log.lambda[i,j])
     log.lambda[i,j] <- alpha[j] + beta[1]*year[i] +
     beta[2]*pow(year[i],2)+
     beta[3]*pow(year[i],3)+
     eps[i] 
} #j
} #i
}
", fill=TRUE)
sink()
jags.data<-list(C=data$C, nsite=data$nsite, nyear=data$nyear, year=data$yr)

# Initial values
inits <- function() list(mu = runif(1, 0, 2), alpha = runif(data$nsite, -1, 1), beta = runif(3, -1, 1), sd.alpha = runif(1, 0, 0.1), sd.year = runif(1, 0, 0.1))

# Parameters monitored (may want to add "lambda")
params <- c("mu", "alpha", "beta", "sd.alpha", "sd.year")

# MCMC settings (may have to adapt)
ni <- 200000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (BRT ~5 min)
out <- jags(jags.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())












