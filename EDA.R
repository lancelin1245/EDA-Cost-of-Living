#+eval=FALSE
library(mvtnorm)
library(coda)
library(MASS)
library(RSQLite)

### Setting up inverse wishart and multivariate normal functions
## mvnormal simulation
rmvnorm<-function(n,mu,Sigma)
{ 
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

## Wishart simulation
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

### Setting up data
db <- dbConnect(RSQLite::SQLite(), "cost_of_living.db")
Cost_of_Living <- dbGetQuery(db, "SELECT * FROM cost_of_living")
Cost_of_Living <- Cost_of_Living[,-1]
dbDisconnect(db)

# attach(Cost_of_Living)

# COLI = cost_of_living_index # Y
# RI = rent_index # X1
# GI = groceries_index # X2
# RIP = restaurant_price_index # X3        
# LPPI = local_purchasing_power_index # X4
# Continent = continent

names = c("Intercept", "Rent Index", "Groceries Index", "Restaurant Price Index", "Local Purchasing Power Index")
continents_names = c("Africa", "Asia", "Europe", "Oceania", "North America", "South America")
# 510 cities in total
# 1) Africa 27 cities
# 2) Asia 107 cities
# 3) Europe 218 cities
# 4) Oceania 13 cities
# America 145 cities
# 5) North America + Caribbean + Central America = 103 + 7 + 11 = 121 cities
# 6) South America 24 cities

### Setting up data
regions<-sort(unique(Cost_of_Living$continent)) 
m<-length(regions)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) {
  Y[[j]]<-Cost_of_Living[Cost_of_Living$continent==regions[j], 2] 
  N[j]<- sum(Cost_of_Living$continent==regions[j])
  x1j<-Cost_of_Living[Cost_of_Living$continent==regions[j], 3] 
  x2j<-Cost_of_Living[Cost_of_Living$continent==regions[j], 4]
  x3j<-Cost_of_Living[Cost_of_Living$continent==regions[j], 5] 
  x4j<-Cost_of_Living[Cost_of_Living$continent==regions[j], 6] 
  X[[j]]<-cbind( rep(1,N[j]), x1j, x2j, x3j, x4j  )
}

#### OLS fits
S2.LS<-BETA.LS<-NULL
# p<-dim(X[[1]])[2]
# BETA.LS <- lapply(1:m, function(x) matrix(NA, ncol=p))
# lapply(BETA.LS, setNames, nm = c("Intercept", "Rent Index", "Groceries Index", "Restaurant Price Index", "Local Purchasing Power Index"))
for(j in 1:m) {
  rent_index = X[[j]][,2]
  groceries_index = X[[j]][,3]
  restaurant_price_index = X[[j]][,4]
  local_purchasing_power_index =  X[[j]][,5]
  fit<-lm(Y[[j]]~rent_index + groceries_index + restaurant_price_index + local_purchasing_power_index)
  BETA.LS<- rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 

#=======================================================================================================
#=======================================================================================================
### Exploratory Data Analysis
# Summary
summary(Cost_of_Living[,-1])

# Pairs plot
library(GGally)
ggpairs(Cost_of_Living[,-1])

# OLS estimates across continents
par(mfrow=c(1,4))
BETA.MLS<-apply(BETA.LS,2,mean)

# Cost of Living index vs Rent index
plot( range(Cost_of_Living[,3]),range(-2, 12),type="n",xlab="Rent Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Rent Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

abline(BETA.MLS[1],BETA.MLS[2],lwd=2)

# Cost of Living index vs Groceries index
plot( range(Cost_of_Living[,4]),range(Cost_of_Living[,2]),type="n",xlab="Groceries Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Groceries Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,3],col="gray")  }

abline(BETA.MLS[1],BETA.MLS[3],lwd=2)

# Cost of Living index vs Restaurant Price index
plot( range(Cost_of_Living[,5]),range(0, 8),type="n",xlab="Restaurant Price Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Restaurant Price Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,4],col="gray")  }

abline(BETA.MLS[1],BETA.MLS[4],lwd=2)

# Cost of Living index vs Local Purchasing Power index
plot( range(Cost_of_Living[,6]),range(BETA.LS[,1]),type="n",xlab="Local Purchasing Power Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Local Purchasing Power Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,5],col="gray")  }

abline(BETA.MLS[1],BETA.MLS[5],lwd=2)

par(mfrow=c(1,5))

# Plot of OLS Estimate variation from the mean across different sample sizes
plot(N,BETA.LS[,1],xlab="sample size",ylab="intercept")
abline(h= BETA.MLS[1],col="black",lwd=2)
plot(N,BETA.LS[,2],xlab="sample size",ylab="OLS estimates for Rent Index")
abline(h= BETA.MLS[2],col="black",lwd=2)
plot(N,BETA.LS[,3],xlab="sample size",ylab="OLS estimates for Groceries Index")
abline(h= BETA.MLS[3],col="black",lwd=2)
plot(N,BETA.LS[,4],xlab="sample size",ylab="OLS estimates for Restaurant Price Index")
abline(h= BETA.MLS[4],col="black",lwd=2)
plot(N,BETA.LS[,5],xlab="sample size",ylab="OLS estimates for Local Purchasing Power Index")
abline(h= BETA.MLS[5],col="black",lwd=2)

par(mfrow=c(1,1))


#=======================================================================================================
#=======================================================================================================
### Start of MCMC
## Setup
S = 510000
p<-dim(X[[1]])[2]
theta<-mu0<-apply(BETA.LS,2,mean)
nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) ; BETA<-as.matrix(BETA.LS)
THETA.b<-S2.b<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps<-matrix(0,p,p)
SIGMA.PS<-NULL
thin = 10
BETA.ps<- z.ps<-lapply(1:m, function(x) matrix(NA, nrow=S/thin, ncol=p))
BETA.pp<-NULL
z <- NULL
for (j in 1:m){
  z <- rbind(z, rep(1, p))
}

# A function to compute the marginal probability
lpy.X<-function(y,X,BETA,silent=TRUE)
{
  n<-dim(X)[1] ; p<-dim(X)[2] 
  (-1/(2*s2))*(t(y)%*%y-2*BETA%*%t(X)%*%y+BETA%*%t(X)%*%X%*%t(BETA))
} # here transpose of Beta is just normal Beta as Beta is a row vector

start.time <- Sys.time()
## MCMC
for(s in 1:S) {
  # sampling z
  for (j in 1:m){
    lpy.c<-lpy.X(Y[[j]],X[[j]][,z[j, ]==1,drop=FALSE],BETA[j,z[j, ]==1,drop=FALSE])
    
    # while (z[j, ] == c(0,0,0,0,0)) { # To get rid of z = (0,0,0,0,0)
      for(i in sample(1:p))
      {
        zp <- z
        zp[j, i] <- 1 - zp[j, i]
        lpy.p<-lpy.X(Y[[j]],X[[j]][,zp[j, ]==1,drop=FALSE],BETA[j,zp[j, ]==1,drop=FALSE])
        r<- (lpy.p - lpy.c)*(-1)^(zp[j,i]==0)
        z[j,i]<-rbinom(1,1,1/(1+exp(-r)))
        if(z[j,i]==zp[j,i]) {lpy.c<-lpy.p}
      }
    # }
  }
  
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma[z[j,]==1, z[j,]==1] + t(X[[j]][, z[j,]==1])%*%X[[j]][,z[j,]==1]/s2 )
    Ej<-Vj%*%( iSigma[z[j,]==1, z[j,]==1]%*%theta[z[j,]==1] + t(X[[j]][, z[j,]==1])%*%Y[[j]]/s2 )
    beta = numeric(p)
    coeffs = rmvnorm(1,Ej,Vj) 
    index = 1
    for (k in 1:length(z[j, ])) {
      if (z[j, k] == 1) {
        beta[k] = coeffs[, index]
        index = index + 1
      }  else {
        beta[k] = 0
      }
    }
    BETA[j,]<-beta
  } 
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]][, z[j,]==1]%*%t(BETA[j, z[j, ]==1,drop=FALSE]))^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  
  ##update theta
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
  theta<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) ) )
  ##
  
  ##store results
  if(s %% thin == 0){
  # cat(s,s2,"\n")
  S2.b<-c(S2.b,s2)
  THETA.b<-rbind(THETA.b,t(theta))
  Sigma.ps<-Sigma.ps+solve(iSigma)
  SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
  BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
  
  # Add posterior draws of Beta for each continent
  for (j in 1:m){
    BETA.ps[[j]][s/thin,] <- BETA[j,]
    z.ps[[j]][s/thin,] <- z[j,]
    }
  }
}
end.time <- Sys.time()
runtime <- end.time - start.time
runtime

# ### Applying burn in period
# # For theta
# THETA.b = THETA.b[1001:nrow(THETA.b),]
# 
# # For Sigma
# SIGMA.PS = SIGMA.PS[1001:nrow(SIGMA.PS),]
# 
# # For Beta
# for (j in 1:m){
#     BETA.ps[[j]] <- BETA.ps[[j]][1001:51000,]
# }
# 
# # For sigma^2
# S2.b = S2.b[1001:length(S2.b)]
# 
# # For z
# for (j in 1:m){
#   z.ps[[j]] <- z.ps[[j]][1001:51000,]
# }

#=======================================================================================================
#=======================================================================================================

### MCMC diagnostics
### Effective Sizes
library(coda)
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])

par(mfrow=c(1,3))
barplot(effectiveSize(THETA.b), main="Effective Sample Sizes of 5 Elements in Theta", ylab = "Effective Size", xlab = "Theta_i")
abline(h=1000)

barplot(apply(SIGMA.PS,2,effectiveSize), main="Effective Sample Sizes of 5*5 Elements in Sigma", ylab = "Effective Size", xlab = "Sigma") 
abline(h=1000)

barplot(effectiveSize(S2.b), main="Effective Sample Sizes of Sigma^2", ylab = "Effective Size", xlab = "Sigma^2") 
abline(h=1000)

par(mfrow=c(2,3))
barplot(apply(BETA.ps[[1]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in Africa", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

barplot(apply(BETA.ps[[2]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in Asia", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

barplot(apply(BETA.ps[[3]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in Europe", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

barplot(apply(BETA.ps[[4]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in Oceania", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

barplot(apply(BETA.ps[[5]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in North America", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

barplot(apply(BETA.ps[[6]],2,effectiveSize), main="Effective Sample Sizes of Coefficients in South America", ylab = "Effective Size", xlab = "Betas") 
abline(h=1000)

### ACF
# Theta
par(mfrow=c(2,3))
acf(THETA.b[,1])
acf(THETA.b[,2])
acf(THETA.b[,3])
acf(THETA.b[,4])
acf(THETA.b[,5])

# Sigma
par(mfrow=c(5,5))
for (i in 1:25) {
  acf(SIGMA.PS[,i])
}

# Sigma^2
par(mfrow=c(1,1))
acf(S2.b)

## Beta0 intercept
par(mfrow=c(2,3))
for (i in 1:m) {
  acf(BETA.ps[[i]][,1], main = paste("ACF plot for Intercept in", continents_names[i]))
}

## Beta1 Rent Index
par(mfrow=c(2,3))
for (i in 1:m) {
  acf(BETA.ps[[i]][,2], main = paste("ACF plot for Rent Index in", continents_names[i]))
}

## Beta2 Groceries Index
par(mfrow=c(2,3))
for (i in 1:m) {
  acf(BETA.ps[[i]][,3], main = paste("ACF plot for Groceries Index", continents_names[i]))
}

## Beta3 Restaurant Price Index
par(mfrow=c(2,3))
for (i in 1:m) {
  acf(BETA.ps[[i]][,4], main = paste("ACF plot for Restaurant Price Index", continents_names[i]))
}

## Beta4 Local Purchasing Power Index
par(mfrow=c(2,3))
for (i in 1:m) {
  acf(BETA.ps[[i]][,5], main = paste("ACF plot for Local Purchasing Power Index", continents_names[i]))
}

### Trace Plots
# Theta
par(mfrow=c(2,3))
for (i in 1:5) {
  plot(THETA.b[,i], main = "Traceplot of THETA")
}

# Sigma
par(mfrow=c(5,5))
for (i in 1:25) {
  plot(SIGMA.PS[,i], main = "Traceplot of SIGMA")
}

# Sigma^2
par(mfrow=c(1,1))
plot(S2.b, main = "Traceplot of SIGMA^2")

## Beta0 intercept
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,1], main = paste("Traceplot for Intercept in", continents_names[i]))
}

## Beta1 Rent Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,2], main = paste("Traceplot for Rent Index", continents_names[i]))
}

## Beta2 Groceries Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,3], main = paste("Traceplot for Groceries Index", continents_names[i]))
}

## Beta3 Restaurant Price Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,4], main = paste("Traceplot for Restaurant Price Index", continents_names[i]))
}

## Beta4 Local Purchasing Power Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,5], main = paste("Traceplot for Local Purchasing Power Index", continents_names[i]))
}


### Marginal posterior densities
par(mfrow = c(1,4))
# Rent Index
plot(density(THETA.b[,2],adj=2),xlim=range(BETA.pp[,2]), 
     main="",xlab="Rent Index Slope Parameter",ylab="Posterior Density",lwd=2)
lines(density(BETA.pp[,2],adj=2),col="gray",lwd=2)
legend(-0.2, 8 ,legend=c( expression(theta[1]),expression(tilde(beta)[1])), 
        lwd=c(2,2),col=c("black","gray"),bty="n") 

# Groceries Index
plot(density(THETA.b[,3],adj=2),xlim=range(BETA.pp[,3]), 
     main="",xlab="Groceries Index Slope Parameter",ylab="Posterior Density",lwd=2)
lines(density(BETA.pp[,3],adj=2),col="gray",lwd=2)
legend(0, 4 ,legend=c( expression(theta[2]),expression(tilde(beta)[2])), 
       lwd=c(2,2),col=c("black","gray"),bty="n") 

# Restaurant Price Index
plot(density(THETA.b[,4],adj=2),xlim=range(BETA.pp[,4]), 
     main="",xlab="Restaurant Price Index Slope Parameter",ylab="Posterior Density",lwd=2)
lines(density(BETA.pp[,4],adj=2),col="gray",lwd=2)
legend(0.6, 6 ,legend=c( expression(theta[3]),expression(tilde(beta)[3])), 
       lwd=c(2,2),col=c("black","gray"),bty="n") 

# Local Purchasing Power Index
plot(density(THETA.b[,5],adj=2),xlim=range(BETA.pp[,5]), 
     main="",xlab="Local Purchasing Power Index Slope Parameter",ylab="Posterior Density",lwd=2)
lines(density(BETA.pp[,5],adj=2),col="gray",lwd=2)
legend(-0.08, 20 ,legend=c( expression(theta[4]),expression(tilde(beta)[4])), 
       lwd=c(2,2),col=c("black","gray"),bty="n") 


### OLS estimates across continents
par(mfrow=c(1,4))

# Cost of Living index vs Rent index
plot( range(Cost_of_Living[,3]),range(0,20),type="n",xlab="Rent Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Rent Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(mean(BETA.ps[[j]][,1]),mean(BETA.ps[[j]][,2]),col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,2]),lwd=2 )

# Cost of Living index vs Groceries index
plot( range(Cost_of_Living[,4]),range(Cost_of_Living[,2]),type="n",xlab="Groceries Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Groceries Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(mean(BETA.ps[[j]][,1]),mean(BETA.ps[[j]][,3]),col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,3]),lwd=2 )

# Cost of Living index vs Restaurant Price index
plot( range(Cost_of_Living[,5]),range(6, 20),type="n",xlab="Restaurant Price Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Restaurant Price Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(mean(BETA.ps[[j]][,1]),mean(BETA.ps[[j]][,4]),col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,4]),lwd=2 )

# Cost of Living index vs Local Purchasing Power index
plot( range(Cost_of_Living[,6]),range(0,20),type="n",xlab="Local Purchasing Power Index", 
      ylab="Cost of Living Index", main = "OLS estimates for Local Purchasing Power Index across Continents", cex.main = 0.8)
for(j in 1:m) {    abline(mean(BETA.ps[[j]][,1]),mean(BETA.ps[[j]][,5]),col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,5]),lwd=2 )

### Posterior predictive checks
BETA.pp.avg.rent < apply(BETA.pp,1,mean) 

#=======================================================================================================
#=======================================================================================================
### Convergence Diagnostics
### Setting up for chain 2
# Chain 2
### Start of MCMC
## Setup
# Change theta, mu0, beta, and z

S = 510000
p<-dim(X[[1]])[2]

nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) 
THETA.b2<-S2.b2<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps2<-matrix(0,p,p)
SIGMA.PS2<-NULL
thin = 10
BETA.ps2<- z.ps2<-lapply(1:m, function(x) matrix(NA, nrow=S/thin, ncol=p))
BETA.pp2<-NULL

# Change these 3 starting values
# theta2<-mu02<-rnorm(mu0, 0.3)
# for (j in 1:m) {
#   BETA2<-rbind(BETA2, rep(0,p))
# }
theta2<-mu02<-apply(BETA.LS,2,mean)-0.05
BETA2 = as.matrix(BETA.LS)-0.05

z2 <- NULL
for (j in 1:m){
  z2 <- rbind(z2, rbinom(p, 1, 0.3))
}

# A function to compute the marginal probability
lpy.X<-function(y,X,BETA,silent=TRUE)
{
  n<-dim(X)[1] ; p<-dim(X)[2] 
  (-1/(2*s2))*(t(y)%*%y-2*BETA%*%t(X)%*%y+BETA%*%t(X)%*%X%*%t(BETA))
} # here transpose of Beta is just normal Beta as Beta is a row vector

start.time <- Sys.time()
## MCMC
for(s in 1:S) {
  # sampling z
  for (j in 1:m){
    lpy.c<-lpy.X(Y[[j]],X[[j]][,z2[j, ]==1,drop=FALSE],BETA2[j,z2[j, ]==1,drop=FALSE])
    
    # while (z2[j, ] == c(0,0,0,0,0)) { # To get rid of z2 = (0,0,0,0,0)
    for(i in sample(1:p))
    {
      zp <- z2
      zp[j, i] <- 1 - zp[j, i]
      lpy.p<-lpy.X(Y[[j]],X[[j]][,zp[j, ]==1,drop=FALSE],BETA2[j,zp[j, ]==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j,i]==0)
      z2[j,i]<-rbinom(1,1,1/(1+exp(-r)))
      if(z2[j,i]==zp[j,i]) {lpy.c<-lpy.p}
    }
    # }
  }
  
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma[z2[j,]==1, z2[j,]==1] + t(X[[j]][, z2[j,]==1])%*%X[[j]][,z2[j,]==1]/s2 )
    Ej<-Vj%*%( iSigma[z2[j,]==1, z2[j,]==1]%*%theta2[z2[j,]==1] + t(X[[j]][, z2[j,]==1])%*%Y[[j]]/s2 )
    beta = numeric(p)
    coeffs = rmvnorm(1,Ej,Vj) 
    index = 1
    for (k in 1:length(z2[j, ])) {
      if (z2[j, k] == 1) {
        beta[k] = coeffs[, index]
        index = index + 1
      }  else {
        beta[k] = 0
      }
    }
    BETA2[j,]<-beta
  } 
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]][, z2[j,]==1]%*%t(BETA2[j, z2[j, ]==1,drop=FALSE]))^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  
  ##update theta2
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu02 + iSigma%*%apply(BETA2,2,sum))
  theta2<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta2<-matrix(theta2,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA2-mtheta2)%*%(BETA2-mtheta2) ) )
  ##
  
  ##store results
  if(s %% thin == 0){
    # cat(s,s2,"\n")
    S2.b2<-c(S2.b2,s2)
    THETA.b2<-rbind(THETA.b2,t(theta2))
    Sigma.ps2<-Sigma.ps2+solve(iSigma)
    SIGMA.PS2<-rbind(SIGMA.PS2,c(solve(iSigma)))
    BETA.pp2<-rbind(BETA.pp2,rmvnorm(1,theta2,solve(iSigma)) )
    
    # Add posterior draws of Beta for each continent
    for (j in 1:m){
      BETA.ps2[[j]][s/thin,] <- BETA2[j,]
      z.ps2[[j]][s/thin,] <- z2[j,]
    }
  }
}
end.time <- Sys.time()
runtime <- end.time - start.time
runtime

### Applying burn in period
# For theta
THETA.b2 = THETA.b2[1001:nrow(THETA.b2),]

# For Sigma
SIGMA.PS2 = SIGMA.PS2[1001:nrow(SIGMA.PS2),]

# For Beta
for (j in 1:m){
  BETA.ps2[[j]] <- BETA.ps2[[j]][1001:51000,]
}

# For sigma^2
S2.b2 = S2.b2[1001:length(S2.b2)]

# For z2
for (j in 1:m){
  z.ps2[[j]] <- z.ps2[[j]][1001:51000,]
}
### Setting up for chain 3
# Chain 3
### Start of MCMC
## Setup
# Change theta, mu0, beta, and z

S = 510000
p<-dim(X[[1]])[2]

nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) 
THETA.b3<-S2.b3<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps3<-matrix(0,p,p)
SIGMA.PS3<-NULL
thin = 10
BETA.ps3<- z.ps3<-lapply(1:m, function(x) matrix(NA, nrow=S/thin, ncol=p))
BETA.pp3<-NULL

# Change these 3 starting values
# theta3<-mu03<-rnorm(mu0, 0.6)
# BETA3 = NULL
# for (j in 1:m) {
# BETA3<-rbind(BETA3, rep(1,p)) 
# }
theta3<-mu03<-apply(BETA.LS,2,mean)+0.05
BETA3 = as.matrix(BETA.LS)+0.05

z3 <- NULL
for (j in 1:m){
  z3 <- rbind(z3, rbinom(p, 1, 0.6))
}

# A function to compute the marginal probability
lpy.X<-function(y,X,BETA,silent=TRUE)
{
  n<-dim(X)[1] ; p<-dim(X)[2] 
  (-1/(2*s2))*(t(y)%*%y-2*BETA%*%t(X)%*%y+BETA%*%t(X)%*%X%*%t(BETA))
} # here transpose of Beta is just normal Beta as Beta is a row vector

start.time <- Sys.time()
## MCMC
for(s in 1:S) {
  # sampling z
  for (j in 1:m){
    lpy.c<-lpy.X(Y[[j]],X[[j]][,z3[j, ]==1,drop=FALSE],BETA3[j,z3[j, ]==1,drop=FALSE])
    
    # while (z3[j, ] == c(0,0,0,0,0)) { # To get rid of z3 = (0,0,0,0,0)
    for(i in sample(1:p))
    {
      zp <- z3
      zp[j, i] <- 1 - zp[j, i]
      lpy.p<-lpy.X(Y[[j]],X[[j]][,zp[j, ]==1,drop=FALSE],BETA3[j,zp[j, ]==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j,i]==0)
      z3[j,i]<-rbinom(1,1,1/(1+exp(-r)))
      if(z3[j,i]==zp[j,i]) {lpy.c<-lpy.p}
    }
    # }
  }
  
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma[z3[j,]==1, z3[j,]==1] + t(X[[j]][, z3[j,]==1])%*%X[[j]][,z3[j,]==1]/s2 )
    Ej<-Vj%*%( iSigma[z3[j,]==1, z3[j,]==1]%*%theta3[z3[j,]==1] + t(X[[j]][, z3[j,]==1])%*%Y[[j]]/s2 )
    beta = numeric(p)
    coeffs = rmvnorm(1,Ej,Vj) 
    index = 1
    for (k in 1:length(z3[j, ])) {
      if (z3[j, k] == 1) {
        beta[k] = coeffs[, index]
        index = index + 1
      }  else {
        beta[k] = 0
      }
    }
    BETA3[j,]<-beta
  } 
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]][, z3[j,]==1]%*%t(BETA3[j, z3[j, ]==1,drop=FALSE]))^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  
  ##update theta3
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu03 + iSigma%*%apply(BETA3,2,sum))
  theta3<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta3<-matrix(theta3,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA3-mtheta3)%*%(BETA3-mtheta3) ) )
  ##
  
  ##store results
  if(s %% thin == 0){
    # cat(s,s2,"\n")
    S2.b3<-c(S2.b3,s2)
    THETA.b3<-rbind(THETA.b3,t(theta3))
    Sigma.ps3<-Sigma.ps3+solve(iSigma)
    SIGMA.PS3<-rbind(SIGMA.PS3,c(solve(iSigma)))
    BETA.p3<-rbind(BETA.pp3,rmvnorm(1,theta3,solve(iSigma)) )
    
    # Add posterior draws of Beta for each continent
    for (j in 1:m){
      BETA.ps3[[j]][s/thin,] <- BETA3[j,]
      z.ps3[[j]][s/thin,] <- z3[j,]
    }
  }
}
end.time <- Sys.time()
runtime <- end.time - start.time
runtime

### Applying burn in period
# For theta
THETA.b3 = THETA.b3[1001:nrow(THETA.b3),]

# For Sigma
SIGMA.PS3 = SIGMA.PS3[1001:nrow(SIGMA.PS3),]

# For Beta
for (j in 1:m){
  BETA.ps3[[j]] <- BETA.ps3[[j]][1001:51000,]
}

# For sigma^2
S2.b3 = S2.b3[1001:length(S2.b3)]

# For z3
for (j in 1:m){
  z.ps3[[j]] <- z.ps3[[j]][1001:51000,]
}
### Setting up for chain 4
# Chain 4
### Start of MCMC
## Setup
# Change theta, mu0, beta, and z

S = 510000
p<-dim(X[[1]])[2]

nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) 
THETA.b4<-S2.b4<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps4<-matrix(0,p,p)
SIGMA.PS4<-NULL
thin = 10
BETA.ps4<- z.ps4<-lapply(1:m, function(x) matrix(NA, nrow=S/thin, ncol=p))
BETA.pp4<-NULL

# Change these 3 starting values
# theta4<-mu04<-rnorm(mu0, 0.9)
# BETA4 = NULL
# for (j in 1:m) {
# BETA4<-rbind(BETA4, rep(2,p)) 
# }

theta4<-mu04<-apply(BETA.LS,2,mean)+0.1
BETA4 = as.matrix(BETA.LS)+0.1

z4 <- NULL
for (j in 1:m){
  z4 <- rbind(z4, rbinom(p, 1, 0.9))
}

# A function to compute the marginal probability
lpy.X<-function(y,X,BETA,silent=TRUE)
{
  n<-dim(X)[1] ; p<-dim(X)[2] 
  (-1/(2*s2))*(t(y)%*%y-2*BETA%*%t(X)%*%y+BETA%*%t(X)%*%X%*%t(BETA))
} # here transpose of Beta is just normal Beta as Beta is a row vector

start.time <- Sys.time()
## MCMC
for(s in 1:S) {
  # sampling z
  for (j in 1:m){
    lpy.c<-lpy.X(Y[[j]],X[[j]][,z4[j, ]==1,drop=FALSE],BETA4[j,z4[j, ]==1,drop=FALSE])
    
    # while (z4[j, ] == c(0,0,0,0,0)) { # To get rid of z4 = (0,0,0,0,0)
    for(i in sample(1:p))
    {
      zp <- z4
      zp[j, i] <- 1 - zp[j, i]
      lpy.p<-lpy.X(Y[[j]],X[[j]][,zp[j, ]==1,drop=FALSE],BETA4[j,zp[j, ]==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j,i]==0)
      z4[j,i]<-rbinom(1,1,1/(1+exp(-r)))
      if(z4[j,i]==zp[j,i]) {lpy.c<-lpy.p}
    }
    # }
  }
  
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma[z4[j,]==1, z4[j,]==1] + t(X[[j]][, z4[j,]==1])%*%X[[j]][,z4[j,]==1]/s2 )
    Ej<-Vj%*%( iSigma[z4[j,]==1, z4[j,]==1]%*%theta4[z4[j,]==1] + t(X[[j]][, z4[j,]==1])%*%Y[[j]]/s2 )
    beta = numeric(p)
    coeffs = rmvnorm(1,Ej,Vj) 
    index = 1
    for (k in 1:length(z4[j, ])) {
      if (z4[j, k] == 1) {
        beta[k] = coeffs[, index]
        index = index + 1
      }  else {
        beta[k] = 0
      }
    }
    BETA4[j,]<-beta
  } 
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]][, z4[j,]==1]%*%t(BETA4[j, z4[j, ]==1,drop=FALSE]))^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  
  ##update theta4
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu04 + iSigma%*%apply(BETA4,2,sum))
  theta4<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta4<-matrix(theta4,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA4-mtheta4)%*%(BETA4-mtheta4) ) )
  ##
  
  ##store results
  if(s %% thin == 0){
    # cat(s,s2,"\n")
    S2.b4<-c(S2.b4,s2)
    THETA.b4<-rbind(THETA.b4,t(theta4))
    Sigma.ps4<-Sigma.ps4+solve(iSigma)
    SIGMA.PS4<-rbind(SIGMA.PS4,c(solve(iSigma)))
    BETA.pp3<-rbind(BETA.pp4,rmvnorm(1,theta4,solve(iSigma)) )
    
    # Add posterior draws of Beta for each continent
    for (j in 1:m){
      BETA.ps4[[j]][s/thin,] <- BETA4[j,]
      z.ps4[[j]][s/thin,] <- z4[j,]
    }
  }
}
end.time <- Sys.time()
runtime <- end.time - start.time
runtime

### Applying burn in period
# For theta
THETA.b4 = THETA.b4[1001:nrow(THETA.b4),]

# For Sigma
SIGMA.PS4 = SIGMA.PS4[1001:nrow(SIGMA.PS4),]

# For Beta
for (j in 1:m){
  BETA.ps4[[j]] <- BETA.ps4[[j]][1001:51000,]
}

# For sigma^2
S2.b4 = S2.b4[1001:length(S2.b4)]

# For z4
for (j in 1:m){
  z.ps4[[j]] <- z.ps4[[j]][1001:51000,]
}
### Trace Plots
# Theta
par(mfrow=c(2,3))
for (i in 1:5) {
  plot(THETA.b[,i], main = "Traceplot of THETA",type='l', col="red")
  lines(THETA.b2[,i], col="blue")
  lines(THETA.b3[,i], col="green")
  lines(THETA.b4[,i], col="purple")
}

# Sigma
par(mfrow=c(5,5))
for (i in 1:25) {
  plot(SIGMA.PS[,i], main = "Traceplot of SIGMA",type='l', col="red")
  lines(SIGMA.PS2[,i], col="blue")
  lines(SIGMA.PS3[,i], col="green")
  lines(SIGMA.PS4[,i], col="purple")
}

# Sigma^2
par(mfrow=c(1,1))
plot(S2.b, main = "Traceplot of SIGMA^2",type='l', col="red")
lines(S2.b2, col="blue")
lines(S2.b3, col="green")
lines(S2.b4, col="purple")


## Beta0 intercept
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,1], main = paste("Traceplot for Intercept in", continents_names[i]),type='l', col="red")
  # lines(BETA.ps2[[i]][,1], col="blue")
  # lines(BETA.ps3[[i]][,1], col="green")
  lines(BETA.ps4[[i]][,1], col="purple")
}

## Beta1 Rent Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,2], main = paste("Traceplot for Rent Index", continents_names[i]),type='l', col="red")
  lines(BETA.ps2[[i]][,2], col="blue")
  lines(BETA.ps3[[i]][,2], col="green")
  lines(BETA.ps4[[i]][,2], col="purple")
}

## Beta2 Groceries Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,3], main = paste("Traceplot for Groceries Index", continents_names[i]),type='l', col="red")
  lines(BETA.ps2[[i]][,3], col="blue")
  lines(BETA.ps3[[i]][,3], col="green")
  lines(BETA.ps4[[i]][,3], col="purple")
}

## Beta3 Restaurant Price Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,4], main = paste("Traceplot for Restaurant Price Index", continents_names[i]),type='l', col="red")
  lines(BETA.ps2[[i]][,4], col="blue")
  lines(BETA.ps3[[i]][,4], col="green")
  lines(BETA.ps4[[i]][,4], col="purple")
}

## Beta4 Local Purchasing Power Index
par(mfrow=c(2,3))
for (i in 1:m) {
  plot(BETA.ps[[i]][,5], main = paste("Traceplot for Local Purchasing Power Index", continents_names[i]),type='l', col="red")
  lines(BETA.ps2[[i]][,5], col="blue")
  lines(BETA.ps3[[i]][,5], col="green")
  lines(BETA.ps4[[i]][,5], col="purple")
}
#=======================================================================================================
#=======================================================================================================
### Confidence intervals
apply(THETA.b, 2, quantile, c(0.025,0.975))

for (i in 1:m) {
  print(paste(paste0("Confidence Interval for coefficients in ", continents_names[i])))
  print(apply(BETA.ps[[i]], 2, quantile, c(0.025,0.975)))
}
#=======================================================================================================
#=======================================================================================================


### Significant covariates
# Across all continents
colnames(z.ps[[1]]) <- names 
colnames(z.ps[[2]]) <- names 
colnames(z.ps[[3]]) <- names 
colnames(z.ps[[4]]) <- names 
colnames(z.ps[[5]]) <- names 
colnames(z.ps[[6]]) <- names 

z.allcontinents <- c()
for (i in 1:p){
  all.z <- c()
  for (j in 1:m){
    all.z <- c(all.z, z.ps[[j]][,i])
  }
  z.allcontinents <- c(z.allcontinents, mean(all.z))
}
z.allcontinents > 0.6

# Across each continent
z.avg <- list()
significant_covariates_continent <- list()
for (j in 1:m){
  z.avg[[j]] <- apply(z.ps[[j]], 2, mean)
  significant_covariates_continent[[j]] <- z.avg[[j]] > 0.6
}
for (j in 1:m){
  print(z.avg[[j]]  > 0.6)
}

### Plot which covariates are significant across continents
par(mfrow=c(1,1))               
plot(z.allcontinents,xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in all Continents")
abline(h=0.6, col='grey')

par(mfrow=c(2,3))  
plot(z.avg[[1]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in Africa")
abline(h=0.6, col='grey')

plot(z.avg[[2]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in Asia")
abline(h=0.6, col='grey')

plot(z.avg[[3]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in Europe")
abline(h=0.6, col='grey')

plot(z.avg[[4]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in Oceania")
abline(h=0.6, col='grey')

plot(z.avg[[5]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in North America")
abline(h=0.6, col='grey')

plot(z.avg[[6]],xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2,
  ylim=c(0,1), main="Overall Significance of Covariates in South America")
abline(h=0.6, col='grey')

### Means of covariates across continents
colnames(BETA.ps[[1]]) <- names 
colnames(BETA.ps[[2]]) <- names 
colnames(BETA.ps[[3]]) <- names 
colnames(BETA.ps[[4]]) <- names 
colnames(BETA.ps[[5]]) <- names 
colnames(BETA.ps[[6]]) <- names 

sort(abs(colMeans(BETA.ps[[1]])))
sort(abs(colMeans(BETA.ps[[2]])))
sort(abs(colMeans(BETA.ps[[3]])))
sort(abs(colMeans(BETA.ps[[4]])))
sort(abs(colMeans(BETA.ps[[5]])))
sort(abs(colMeans(BETA.ps[[6]])))

