############################################################################
# tables
############################################################################

source('xtab.r') # Yang Feng's package

# sensetivity #####################################################################

 

###############table 1########################### 
load(file='result/sensitivity.Rdata')  

mean11  <-  apply(sensitivity11$sigmahat,2,mean)  
se11  <-  apply(sensitivity11$sigmahat,2,sd)

mean12 <-  apply(sensitivity12$sigmahat,2,mean)  
se12 <-  apply(sensitivity12$sigmahat,2,sd)

mean13 <-  apply(sensitivity13$sigmahat,2,mean)  
se13 <-  apply(sensitivity13$sigmahat,2,sd)

mean21 <-  apply(sensitivity21$sigmahat,2,mean)  
se21 <-  apply(sensitivity21$sigmahat,2,sd)

mean22 <-  apply(sensitivity22$sigmahat,2,mean)  
se22 <-  apply(sensitivity22$sigmahat,2,sd)

mean23 <-  apply(sensitivity23$sigmahat,2,mean)  
se23 <-  apply(sensitivity23$sigmahat,2,sd)

mean31 <-  apply(sensitivity31$sigmahat,2,mean)  
se31 <-  apply(sensitivity31$sigmahat,2,sd)

mean32 <-  apply(sensitivity32$sigmahat,2,mean)  
se32 <-  apply(sensitivity32$sigmahat,2,sd)

mean33 <-  apply(sensitivity33$sigmahat,2,mean)  
se33 <-  apply(sensitivity33$sigmahat,2,sd)


xmean <-  rbind(mean11,mean12,mean13,mean21,mean22,mean23,mean31,mean32,mean33) 
xse <-  rbind(se11,se12,se13,se21,se22,se23,se31,se32,se33) 

colnames(xmean) <- c("K=5","K=10", "K=15", "K=20", "tuned",  "Oracle")
rownames(xmean) <- c("S1-G", "S1-T", "S1-E", "S2-G", "S2-T", "S2-E", "S3-G", "S3-T", "S3-E")
xtab(xmean,xse,digits=3)

# relative efficency
#mse <-  (xmean-1)^2+xse^2
#mse[,1:5]/mse[,6]

# choice of K

sum(sensitivity31$selectedK==10)/500       #[1] 0.968
sum(sensitivity32$selectedK==10)/500       #[1] 0.96
sum(sensitivity33$selectedK==10)/500       #[1] 0.952


# comparison #####################################################################
load(file='result/sigmahat.Rdata') 


###############table 2########################### 

mean11  <-  apply(sigmahat11,2,mean)  
se11 <-  apply(sigmahat11,2,sd)

mean12 <-  apply(sigmahat12,2,mean)  
se12 <-  apply(sigmahat12,2,sd)

mean13 <-  apply(sigmahat13,2,mean)  
se13 <-  apply(sigmahat13,2,sd)

mean21 <-  apply(sigmahat21,2,mean)  
se21 <-  apply(sigmahat21,2,sd)

mean22 <-  apply(sigmahat22,2,mean)  
se22 <-  apply(sigmahat22,2,sd)

mean23 <-  apply(sigmahat23,2,mean)  
se23 <-  apply(sigmahat23,2,sd)

mean31 <-  apply(sigmahat31,2,mean)  
se31 <-  apply(sigmahat31,2,sd)

mean32 <-  apply(sigmahat32,2,mean)  
se32 <-  apply(sigmahat32,2,sd)

mean33 <-  apply(sigmahat33,2,mean)  
se33 <-  apply(sigmahat33,2,sd)

xmean <-  rbind(mean11,mean12,mean13,mean21,mean22,mean23,mean31,mean32,mean33) 
xse <-  rbind(se11,se12,se13,se21,se22,se23,se31,se32,se33) 
colnames(xmean) <- c("EVE","EVE(K=10)", "MS(K=10)", "MAD", "DK", "Rice", "Oracle")
rownames(xmean) <- c("S1-G", "S1-T", "S1-E", "S2-G", "S2-T", "S2-E", "S3-G", "S3-T", "S3-E")
xtab(xmean,xse,digits=3)

###############table 3########################### 

# relative efficency
mse <-  (xmean-1)^2+xse^2
relEf <-  mse[,1:6]/mse[,7]
xtab(relEf,digits=2)

# relative efficiency
# distance to the oracle estimator


# noise from real data #####################################################################
load(file='result/realdata.Rdata') 

###############table 4########################### 

mean1  <-  apply(sigmahat1,2,mean)  
se1  <-  apply(sigmahat1,2,sd)

mean2  <-  apply(sigmahat2,2,mean)  
se2  <-  apply(sigmahat2,2,sd)

mean3  <-  apply(sigmahat3,2,mean)  
se3 <-  apply(sigmahat3,2,sd)


xmean  <-  rbind(mean1,mean2,mean3) 
xse  <-  rbind(se1,se2,se3) 
colnames(xmean) <- c("EVE","EVE(K=10)", "MS(K=10)", "MAD", "DK", "Rice", "Oracle")
rownames(xmean) <- c("S1", "S2", "S3")
xtab(xmean,xse,digits=3)

###############table 5########################### 

# relative efficiency
mse <-  (xmean-1)^2+xse^2
relEf <-  mse[,1:6]/mse[,7]
xtab(relEf,digits=2)


# Real data #####################################################################
load(file='result/realdata2.Rdata') 

###############table 6########################### 
xtab(var.tab, digit=2)


############################################################################
# figure 1
############################################################################
utm1 <- function(k){
  aa <- matrix(0,k,k)
  for (i in 1:k){
    aa[i,i:k] <- 1
  }
  return(aa)
}

K <- 15
X <- cbind(rep(1,K),1:K)
U <- utm1(K)
lambda <- 0:1600/400
rs <- array(0,c(3,length(lambda)))
for (i in 1:length(lambda)){
  Sigma <- diag(K)+lambda[i]*t(U)%*%U
  Sigma.inv <- solve(Sigma)
  temp <- t(X)%*%Sigma.inv%*%X
  rs[1,i] <- (solve(t(X)%*%X)%*%t(X)%*%Sigma%*%X%*%solve(t(X)%*%X))[1,1] ## OLS on Theta_2L
  rs[2,i] <- solve(temp)[1,1] ## GLS on Theta_2L
  Sigma <- diag(K)+2*lambda[i]*t(U)%*%U
  Sigma.inv <- solve(Sigma)
  temp <- t(X)%*%Sigma.inv%*%X
  rs[3,i] <- solve(temp)[1,1] ## Upper bound on Theta_L using g()
}
yr <- range(c(2,rs+2))
par(mar=c(4,4,.5,.5))
plot(lambda/2,rs[1,]+2,type="l",ylim=yr,lty=3,xlab=expression(w),ylab=expression(r)) ## OLS-2L
legend("topleft", legend=c("OLS-2L","OLS-L", "GLS-2L","GLS-L"),
       lty=c(3,1,4,5), cex=0.8)
lines(lambda/2,rs[2,]+2,type="l",lty=4) ## GLS-2L
lines(lambda/2,rs[3,]+2,type="l",lty=5) ## GLS-L
b1 <- (K+1)*(K+2)*(2*K+1)/15/K/(K-1)
b2 <- K/(K-1)
b3 <- (K+1)*(K+2)^2/6/K/(K-1)
#abline(rs[1,1],b1)
#abline(rs[1,1],b2)
#abline(rs[1,1]+2,2*b3,lty=1)
lines(range(lambda/2),c(rs[1,1]+2,rs[1,1]+2+2*b3*max(lambda/2)),lty=1) ## OLS-L
abline(h=2) ## Oracle
