## Capture-recapture 1

#Define parameters
n.occasions<-6
marked<-rep(50, n.occasions-1) # Annual number of marked individuals
phi<-rep(0.88, n.occasions-1)

p<-rep(0.75, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break		# If dead, move to next individual 
         # Bernoulli trial: is individual recaptured? 
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }
   
CH<-simul.cjs(PHI=PHI, P=P, marked=marked)

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

sink('cjs-c-c.jags')
cat("
model {
# priors
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
  phi[i,t]<-mean.phi
  p[i,t]<-mean.p
  } #t
} #i
mean.phi ~ dunif(0,1)
mean.p ~ dunif(0,1)
# Likelihood
for (i in 1:nind) {
  #define the latent state 
  z[i,f[i]]<-1
  for (t in (f[i]+1):n.occasions) {
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t]<-phi[i, t-1]*z[i, t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <-p[i, t-1]*z[i,t]
    } #t 
} #i
}
", fill=TRUE)
sink()

jags.data<-list(y=CH, f=f, nind = dim(CH)[1], n.occasions=dim(CH)[2])

known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

   known.state.cjs(CH)
# initials 
inits<-function() {list(mean.phi=runif(1, 0,1), mean.p=runif(1, 0,1), z=   known.state.cjs(CH)
)}


# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)

##########################################
# better estimation
jags.data<-list(y=CH, f=f, nind = dim(CH)[1], n.occasions=dim(CH)[2], z=known.state.cjs(CH))


# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT <1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


#########################################################
############ random effects models

#Define parameters
n.occasions<-6
marked<-rep(30, n.occasions-1) # Annual number of marked individuals
mean.phi<-0.88
var.phi<-1
p<-rep(0.75, n.occasions-1)

logit.phi<-rnorm(n.occasions-1, logit(mean.phi), var.phi^2)
phi<-inv.logit(logit.phi)
phi
# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-temp-raneff.jags")
cat("
model {
# priors
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
  logit(phi[i,t])<- mu + epsilon[t]
  p[i,t]<-mean.p
  } #t
} #i
 for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
  }  

mean.phi ~ dunif(0,1)
mean.p ~ dunif(0,1)
mu <- log(mean.phi/(1-mean.phi)) # logit transform
sigma ~ dunif(0,3)
tau <- pow(sigma, -2)
sigma2<-pow(sigma,2)


# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
", fill=TRUE)
sink()

##########################################
# better estimation
jags.data<-list(y=CH, f=f, nind = dim(CH)[1], n.occasions=dim(CH)[2], z=known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}  
# Parameters monitored
parameters <- c("mean.phi", "mean.p", 'sigma2')

# MCMC settings
ni <- 5000
nt <- 6
nb <- 4000
nc <- 7

# Call JAGS from R (BRT <1 min)
cjs.c.c <- jags.parallel(jags.data, inits, parameters, "cjs-temp-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd(), export_obj_names=c('cjs.init.z', 'nb', 'ni', 'nt', 'CH'))

download.file('https://git.io/fp2ER', destfile='cjs-temp-raneff_cjs.ran.RData', mode='wb')
load('cjs-temp-raneff_cjs.ran.RData')
cjs.c.c <-cjs.ran
cjs.c.c

################################################
# covariate and random effects

#Define parameters
n.occasions<-20
marked<-rep(15, n.occasions-1) # Annual number of marked individuals
mean.phi<-0.88
r.var<-0.2
beta <- -0.3                       # Slope of survival-winter relationship	
p<-rep(0.75, n.occasions-1)

winter <- rnorm(n.occasions-1, 0, 1)
logit.phi<- logit(mean.phi) + beta*winter + rnorm(n.occasions-1, 0 , r.var^0.5)

phi<-inv.logit(logit.phi)
phi
# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-cov-raneff.jags")
cat("
model {
# priors
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
  logit(phi[i,t])<- mu+ beta*x[t]+ epsilon[t]
  p[i,t]<-mean.p
  } #t
} #i
 for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    phi.est[t]<- 1/(1+exp(-mu- beta*x[t]- epsilon[t]))
  }  

mean.phi ~ dunif(0,1)
mean.p ~ dunif(0,1)
mu <- log(mean.phi/(1-mean.phi)) # logit transform
beta ~ dnorm(0, 0.001)I(-10,10)
sigma ~ dunif(0,5)
tau <- pow(sigma, -2)
sigma2<-pow(sigma,2)


# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i

}
", fill=TRUE)
sink()

##########################################
# better estimation
jags.data<-list(y=CH, f=f, nind = dim(CH)[1], n.occasions=dim(CH)[2], z=known.state.cjs(CH), x=winter)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma=runif(1, 0, 5), beta=runif(1, -5,5))}  
# Parameters monitored

parameters <- c("mean.phi", "mean.p", 'sigma2', 'beta', 'phi.est')

# MCMC settings
ni <- 1000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 12 min)
cjs.cov <- jags(jags.data, inits, parameters, "cjs-cov-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


ni <- 12000
nt <- 6
nb <- 8000
nc <- 7

# Call JAGS from R (BRT 7 min)
cjs.cov <- jags.parallel(jags.data, inits, parameters, "cjs-cov-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd(), , export_obj_names=c('cjs.init.z', 'nb', 'ni', 'nt', 'CH'))

download.file('https://git.io/fp2E2', destfile='cjs-cov-raneff_cjs.cov.RData', mode='wb')
load('cjs-cov-raneff_cjs.cov.RData')


plot(cjs.cov$BUGSoutput$mean$phi.est, ylim=c(0,1), pch=19)
points(phi, col='red', pch=19)



###############################################
#  sex specific survival

# Models with individual variation
# Fixed group effects - sex specific survival
# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi.f <- rep(0.92, n.occasions-1)  # Survival of females
p.f <- rep(0.7, n.occasions-1)     # Recapture of females
phi.m <- rep(0.87, n.occasions-1)   # Survival of males
p.m <- rep(0.7, n.occasions-1)     # Reacpture of males

# Define matrices with survival and recapture probabilities
PHI.m <- matrix(phi.m, ncol = n.occasions-1, nrow = sum(marked))
PHI.f<- matrix(phi.f, ncol = n.occasions-1, nrow = sum(marked))
P.m <- matrix(p.m, ncol = n.occasions-1, nrow = sum(marked))
P.f <- matrix(p.f, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH.m <- simul.cjs(PHI.m, P.m, marked)
CH.f <- simul.cjs(PHI.f, P.f, marked)

# Merge capture-histories by row
CH <- rbind(CH.f, CH.m)

group<-c(rep(1, dim(CH.f)[1]), rep(2, dim(CH.m)[1]))

f <- apply(CH, 1, get.first)

sink('cjs-group.jags')
cat("
model {
# priors
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
  phi[i,t]<-phi.g[group[i]]
  p[i,t]<-p.g[group[i]]
  } #t
} #i

for(u in 1:g) {
 phi.g[u] ~ dunif(0,1)
 p.g[u] ~ dunif(0,1)
}

# Likelihood
for (i in 1:nind) {
  #define the latent state 
  z[i,f[i]]<-1
  for (t in (f[i]+1):n.occasions) {
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t]<-phi[i, t-1]*z[i, t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <-p[i, t-1]*z[i,t]
    } #t 
} #i
}
", fill=TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length(unique(group)), 0, 1), p.g = runif(length(unique(group)), 0, 1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 2 min)
cjs.group <- jags(jags.data, inits, parameters, "cjs-group.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

cjs.group

#################################################### group difference

group<-c(rep(0, dim(CH.f)[1]), rep(1, dim(CH.m)[1]))

f <- apply(CH, 1, get.first)

sink('cjs-group-difference.jags')
cat("
model {
# priors
for (i in 1:nind) {
  for (t in f[i]:(n.occasions-1)) {
  logit(phi[i,t])<-phi.g+phi.male.difference*group[i]
  logit(p[i,t])<-p.g+p.male.difference*group[i]
  } #t
} #i

phi.g~dunif(-5,5)
p.g~dunif(-5,5)
phi.male.difference ~ dunif(-1,1)
p.male.difference ~ dunif(-1,1)

# Likelihood
for (i in 1:nind) {
  #define the latent state 
  z[i,f[i]]<-1
  for (t in (f[i]+1):n.occasions) {
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t]<-phi[i, t-1]*z[i, t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <-p[i, t-1]*z[i,t]
    } #t 
} #i
}
", fill=TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),  group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(1, -3, 3), p.g =  runif(1, -3, 3), phi.male.difference=runif(1, -1,1), p.male.difference=runif(1, -1,1))}  

# Parameters monitored
parameters <- c("phi.g", "p.g", "phi.male.difference", "p.male.difference")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 2 min)
cjs.group <- jags(jags.data, inits, parameters, "cjs-group-difference.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

cjs.group

####################################################
# age effects ibn survival
n.occasions<-10
marked.j<-rep(200, n.occasions-1)
marked.a<-rep(30, n.occasions-1)
phi.j<-0.3
phi.ad<-0.7
p<-rep(0.5, n.occasions-1)

PHI.J<-matrix(0, ncol=n.occasions-1, nrow=sum(marked.j))
for (i in 1:length(marked.j)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
  PHI.J[is.na(PHI.J[,i]),i]<-phi.ad
   }
   
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.ad, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)

# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.J)[1]){
   for (t in f.j[i]:(dim(CH.J)[2]-1)){
      x.j[i,t] <- 2
      x.j[i,f.j[i]] <- 1   
      } #t
   } #i
for (i in 1:dim(CH.A)[1]){
   for (t in f.a[i]:(dim(CH.A)[2]-1)){
      x.a[i,t] <- 2
      } #t
   } #i

CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

# Specify model in BUGS language
sink("cjs-age.jags")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
      } #t
   } #i
for (u in 1:2){
   beta[u] ~ dunif(0, 1)              # Priors for age-specific survival
   }
mean.p ~ dunif(0, 1)                  # Prior for mean recapture
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.age <- jags(jags.data, inits, parameters, "cjs-age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())





P.J<-matrix(0, ncol=n.occasions-1, nrow=sum(marked.j))












