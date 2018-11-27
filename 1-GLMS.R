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
            
            
            
            
            
            
            
            




