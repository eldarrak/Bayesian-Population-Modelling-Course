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
             sd.site=sd.site,
             sd.year=sd.year,
             expected.count=expected.count,
             C=C))
}

data<-data.fn()
data
