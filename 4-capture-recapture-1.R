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
  for (t in f[i]:n.occasions-1) {
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

