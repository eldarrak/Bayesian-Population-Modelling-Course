# 1. Linear models and Poisson GLM in R and JAGS
## required data falcons.txt

y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(X, y, col = c(rep("red", 3), rep("blue", 3), rep("green", 3)), xlim = c(-1, 25), ylim = c(0, 140), pch=19, cex=3)

LM1<-lm(y~A-1+X)
summary(LM1)

model.matrix(~A-1+X)
model.matrix(~A+X)

LM2<-lm(y~A+X)
summary(LM2)

SUM<-function(a=1,b=3) {
 c=a+b
 return(c)
}

SUM(5,2)

data.fn<-function(n=40, alpha=3.5576, beta1=-0.0912, beta2=0.0091, beta3= -0.00014 ) {
# n number of Years
# alpha, beta, beta2, beta3 
#  Generate sequence
years<-1:n

log.expected.count<-alpha + 
                  beta1*years +
                  beta2*years^2 + 
                  beta3*years^3
expected.count<-exp(log.expected.count)

C<-rpois(n=n, lambda=expected.count)
return(list(n=n, alpha=alpha, beta1=beta1, beta2=beta2, beta3=beta3, years=years, expected.count=expected.count, C=C))
}

                  
data<-data.fn()

data

plot(data$C~data$years)

GLM1<-glm(C~ years + 
            I(years^2)+
            I(years^3) , data= data, family=poisson)
summary(GLM1)            
            
# Specify model in JAGS language
getwd()
# Specify model in BUGS language
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * years[i] + beta2 * pow(years[i],2) + beta3 * pow(years[i],3)                      # 3. Linear predictor
   } #i
}
",fill = TRUE)
sink()


# Bundle data
mean.year <- mean(data$years)             # Mean of year covariate
sd.year <- sd(data$years)                 # SD of year covariate
jags.data <- list(C = data$C, n = length(data$C), years = (data$years - mean.year) / sd.year)

   
# Initial values   
inits<-function() {
      list(
      alpha = runif(1, -2, 2), 
      beta1 = runif(1,-1,1),
      beta2 = runif(1, -1,1),
      beta3 = runif(1, -1,1)
)}
      
# Parameters to monitor      
params<-c('alpha', 'beta1', 'beta2', 'beta3', 'lambda')
      
#MCMC settings            
 
ni<-2000
nt=10
nb=80
nc=3

library(R2jags)
out<-jags(data=jags.data,
          model.file = "GLM_Poisson.jags", 
          inits=inits, 
          parameters.to.save=params,
          n.chains=nc, n.thin=nt,
          n.burnin=nb, n.iter=ni)

 out

traceplot(out) 


# Call JAGS from R (BRT < 1 min)
tmp <- jags(data = jags.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ years + I(years^2) + I(years^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")
JAGS.predictions <- out$BUGSoutput$mean$lambda
lines(1:40, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, JAGS.predictions)

#####################################################
# 
peregrine <- read.table("falcons.txt", header = TRUE)

attach(peregrine)

plot(Year, Pairs, type = "b", lwd = 2, main = "", las = 1, ylab = "Pair count", xlab = "Year", ylim = c(0, 200), pch = 16)

# Bundle data
mean.year <- mean(1:length(Year))        # Mean of year covariate
sd.year <- sd(1:length(Year))            # SD of year covariate

jags.data<-list(C=Pairs, n=length(Pairs), years=(1:length(Year)-mean.year)/sd.year)

 
# Initial values   
inits<-function() {
      list(
      alpha = runif(1, -2, 2), 
      beta1 = runif(1,-1,1),
      beta2 = runif(1, -1,1),
      beta3 = runif(1, -1,1)
)}
      
# Parameters to monitor      
params<-c('alpha', 'beta1', 'beta2', 'beta3', 'lambda')
      
#MCMC settings            
 
ni<-2000
nt=10
nb=80
nc=3

library(R2jags)
out_falcons<-jags(data=jags.data,
          model.file = "GLM_Poisson.jags", 
          inits=inits, 
          parameters.to.save=params,
          n.chains=nc, n.thin=nt,
          n.burnin=nb, n.iter=ni)

#####################################################
##################################################
#  Binomial GLM
# 
data.fn<-function(nyears=40, alpha=0, beta1=-0.1, beta2=-0.9) {

years=1:nyears

Yr<-(years-mean(years))/sd(years)

N<-round(runif(nyears, min=20, max=100))

exp.p<-plogis(alpha + 
              beta1*Yr +
              beta2*Yr^2)
C<-rbinom(n=nyears, size=N, prob=exp.p)

return(list(nyears = nyears, alpha = alpha, beta1 = beta1, beta2 = beta2, year = years, YR = Yr, exp.p = exp.p, C = C, N = N))
}

data<-data.fn()














}





