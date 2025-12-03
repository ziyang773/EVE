library(EVE)

################################################################################
# Other variance (standard deviation) methods used for comparison
################################################################################
var_est_methods <- function(x){
  n <- length(x)
  sigmahat1 <- 1.4826*median(abs(x-median(x)))             # MAD
  sigmahat2 <- 1.48*sqrt(0.5)*median(abs(x[-1]-x[-n]))     # DK
  sigmahat3 <- sqrt(mean((x[-1]-x[-n])^2)/2)               #Rice
  sigmahat <- c(sigmahat1,sigmahat2,sigmahat3)
  names(sigmahat) = c("MAD", "DK", "Rice")
  return(list(sigmahat = sigmahat, K = names(sigmahat)))
}

################################ helpers #######################################
################################################################################
#     sequence generator under 3 scenarios with 3 noise types.
################################################################################

seqgenerator = function(scenario = 1,noise = 1,n = 1000){
  if (noise == 1) ve <- rnorm(n)                            #kappa4=3
  if (noise == 2) ve <- rt(n,df=6)/sqrt(1.5)                #kappa4=6
  if (noise == 3) ve <- rexp(n)-1                           #Asymmetry
  if (scenario == 1) mu <- rep(0,n)
  if (scenario == 2) {
     mu <- rep(0,n)
     mu[c(101:110,201:210,301:310,401:410,501:510,601:610)] <- 1
     mu[801:820] <- -3
  }
  if (scenario == 3) {
     mu <- rep(c(rep(1,10),rep(-1,10)),n/20)
  }
  return(list(mu = mu, ve = ve, x = mu+ve))
}

################################################################################
# All methods used for numerical analysis
################################################################################

varMCM = function(x){
  ob1 <- as.data.frame(eve(x))
  ob2 <- as.data.frame(eve(x, K = 10))
  ob3 <- as.data.frame(MSve(x, K = 10))
  ob4 <- as.data.frame(var_est_methods(x))
  return(rbind(ob1,ob2,ob3,ob4))
}

############################### Simulation ######################################
################################################################################
# Sensitivity Analysis: compare the bias/variance of estimators using different K's, Table 1
################################################################################

simulator1 <- function(scenario = 1, noise = 1, K = c(5,10,15,20), simT = 500, n = 1000, seed = 1){
   set.seed(seed)
   sigmahat <- matrix(0,simT,length(K)+2)
   selectedK <- numeric(simT)
   for (i in 1:simT){
     temp <- seqgenerator(scenario,noise,n)
     x <- temp$x
     ve <- temp$ve
     sigmahat[i,1:length(K)] <- eve(x,K)$sigmahat
     selKobj <- eve(x, K = NULL, Kmax = 20)
     sigmahat[i,length(K)+1] <-  selKobj$sigmahat
     sigmahat[i,length(K)+2] <- sd(ve)
     selectedK[i] <- selKobj$K
   }
   colnames(sigmahat) <- c(K,'tuned','Oracle')
   bias <- apply(sigmahat,2,mean)-1
   se <- apply(sigmahat,2,sd)
   mse <- bias^2+se^2
   # print(rbind(bias,se,mse))
   return(list(sigmahat=sigmahat,selectedK=selectedK))
}

## labels for noise types (adjust if needed)
noise_lab <- c("G", "T", "E")
## all (scenario, noise) combinations
pars <- expand.grid(noise = 1:3, scenario = 1:3)

## run all simulations
res_list <- lapply(seq_len(nrow(pars)), function(i) {
  simulator1(scenario = pars$scenario[i], noise = pars$noise[i])
})

# get mean and se
extract_mean_se <- function(sim_obj) {
  mean_vec <- apply(sim_obj$sigmahat, 2, mean)
  se_vec   <- apply(sim_obj$sigmahat, 2, sd)
  list(mean = mean_vec, se = se_vec)
}
stats_list <- lapply(res_list, extract_mean_se)

## row names like S1-G, S1-T, S1-E, ...
names(stats_list) <- paste0("S", pars$scenario, "-", noise_lab[pars$noise])

## turn each result into "mean(se)" strings and stack them
table1 <- do.call(rbind, lapply(stats_list, function(mat) {
  mean_est <- mat$mean
  se_est   <- mat$se
  sprintf("%.3f(%.3f)", mean_est, se_est)
}))

colnames(table1) <- names(stats_list[[1]]$mean)
table1

################################################################################
# Comparison of EVE against other variance estimators, Table 2-3
################################################################################

simulator2  <-  function(scenario=1,noise=1, simT=500, n=1000, seed=1){
   set.seed(seed)
   sigmahat=matrix(0,simT,7)
   for (i in 1:simT){
     temp <-  seqgenerator(scenario,noise,n)
     x <-  temp$x
     ve <-  temp$ve
     obj <- varMCM(x)
     sigmahat[i,1:6] <- obj$sigmahat
     sigmahat[i,7] <- sd(ve)
   }
   bias <-  apply(sigmahat,2,mean)-1
   se <-  apply(sigmahat,2,sd)
   mse <-  bias^2+se^2
   # print(rbind(bias,se,mse))
   return(sigmahat)
}

res_list <- lapply(seq_len(nrow(pars)), function(i) {
  simulator2(scenario = pars$scenario[i], noise = pars$noise[i])
})

extract_mean_se2 <- function(mat) {
  mean_vec <- colMeans(mat)
  se_vec   <- apply(mat, 2, sd)
  list(mean = mean_vec, se = se_vec)
}

stats_list <- lapply(res_list, extract_mean_se2)

## row names like S1-G, S1-T, S1-E, ...
names(stats_list) <- paste0("S", pars$scenario, "-", noise_lab[pars$noise])

## turn each result into "mean(se)" strings and stack them
table2 <- do.call(rbind, lapply(stats_list, function(mat) {
  mean_est <- mat$mean
  se_est   <- mat$se
  sprintf("%.3f(%.3f)", mean_est, se_est)
}))

colnames(table2) <- c('EVE', 'EVE(K=10)', 'MS(K=10)', 'MAD', 'DK', 'Rice', 'Oracle')
table2

## means and sds from the same sims used for table2
xmean <- do.call(rbind, lapply(stats_list, `[[`, "mean"))
xse   <- do.call(rbind, lapply(stats_list, `[[`, "se"))

## MSE and relative efficiency vs Oracle
mse   <- (xmean - 1)^2 + xse^2
relEf <- sweep(mse[, 1:6, drop = FALSE], 1, mse[, 7], "/")

# table 3
table3 <- round(relEf, 2)
colnames(table3) <- colnames(table2)[1:6]
table3

################################################################################
#     Noise from real data, Table 4-5
################################################################################
load("trios.Rdata")

# # we focus on chr11 of the father
Chr       <-  father[,2]
father11 <-  father[which(Chr==11),5]

newfather11 <- father11[-which(is.na(father11))]
newfather11 <- scale(newfather11)

simulator3 = function(scenario = 1, noise = NULL, simT = 500, n = 1000, seed = 1){
   if (scenario==1) mu=rep(0,n)
   if (scenario==2) {
      mu <- rep(0,n)
      mu[c(101:110,201:210,301:310,401:410,501:510,601:610)]=1
      mu[801:820]=-3
   }
   if (scenario==3) {
      mu <- rep(c(rep(1,10),rep(-1,10)),n/20)
   }
   set.seed(seed)
   sigmahat <- matrix(0,simT,7)
   for (i in 1:simT){
     ve  <-  noise[sample(length(noise),n)]
     x  <-  mu+ve
     sigmahat[i,1:6] <- varMCM(x)$sigmahat
     sigmahat[i,7] <- sd(ve)
   }
   bias  <-  apply(sigmahat,2,mean)-1
   se  <-  apply(sigmahat,2,sd)
   mse  <-  bias^2+se^2
   # print(rbind(bias,se,mse))
   return(sigmahat)
}

res_list <- lapply(1:3, function(s) {
  simulator3(scenario = s, noise = newfather11)
})

## means and sds across the 500 replicates
mean_mat <- do.call(rbind, lapply(res_list, colMeans))
se_mat   <- do.call(rbind, lapply(res_list, function(m) apply(m, 2, sd)))

## table 4
table4 <- t(sapply(1:3, function(i) {
  sprintf("%.3f(%.3f)", mean_mat[i, ], se_mat[i, ])
}))

rownames(table4) <- paste0("S", 1:3)
colnames(table4) <- c("EVE", "EVE(K=10)", "MS(K=10)", "MAD", "DK", "Rice", "Oracle")
table4

## mse and relative efficiency VS Oracle
mse <- (mean_mat - 1)^2 + se_mat^2
relEf <- mse[, 1:6, drop = FALSE] / mse[, 7]

table5 <- round(relEf, 2)
rownames(table5) <- paste0("S", 1:3)
colnames(table5) <- c("EVE", "EVE(K=10)", "MS(K=10)", "MAD", "DK", "Rice")
table5

################################################################################
#    Real data, table 6, Figure 2
################################################################################

aaa <- read.csv("lp1.csv",header=T)
aaa <- aaa[1:131,]
n <- dim(aaa)[1]
DUR <- aaa$DUR
NDUR <- aaa$NDUR
bbb <- read.csv("lp2.csv",header=T)
bbb <- bbb[161:291,]
NFC <- bbb$NFC
NFBUS <- bbb$NFBUS
BUS <- bbb$BUS

var.tab <- c(varMCM(DUR)$sigmahat[c(1,4,5,6)],sd(DUR))
var.tab <- rbind(var.tab, c(varMCM(NDUR)$sigmahat[c(1,4,5,6)], sd(NDUR)))
var.tab <- rbind(var.tab, c(varMCM(BUS)$sigmahat[c(1,4,5,6)], sd(BUS)))
var.tab <- rbind(var.tab, c(varMCM(NFBUS)$sigmahat[c(1,4,5,6)], sd(NFBUS)))
var.tab <- rbind(var.tab, c(varMCM(NFC)$sigmahat[c(1,4,5,6)], sd(NFC)))
rownames(var.tab) <- c("DUR","NDUR","BUS","NFBUS","NFC")
colnames(var.tab) <- c("EVE", "MAD", "DK", "Rice", "SD")
table6 <- round(var.tab, 2)
table6

# # Figure 2
var.tab <- c(as.numeric(varMCM(DUR)$K[1]),varMCM(DUR)$sigmahat[c(1,4,5,6)],sd(DUR))
var.tab <- rbind(var.tab, c(as.numeric(varMCM(NDUR)$K[1]),varMCM(NDUR)$sigmahat[c(1,4,5,6)], sd(NDUR)))
var.tab <- rbind(var.tab, c(as.numeric(varMCM(BUS)$K[1]),varMCM(BUS)$sigmahat[c(1,4,5,6)], sd(BUS)))
var.tab <- rbind(var.tab, c(as.numeric(varMCM(NFBUS)$K[1]),varMCM(NFBUS)$sigmahat[c(1,4,5,6)], sd(NFBUS)))
var.tab <- rbind(var.tab, c(as.numeric(varMCM(NFC)$K[1]),varMCM(NFC)$sigmahat[c(1,4,5,6)], sd(NFC)))
rownames(var.tab) <- c("DUR","NDUR","BUS","NFBUS","NFC")
colnames(var.tab) <- c("K","EVE", "MAD", "DK", "Rice", "SD")

par(mfrow=c(5,2))
par(mar=c(.5,.5,.5,.5))

plot(DUR,type="l",xaxt="n",yaxt="n")
legend("topleft",legend=c("DUR "), cex = 0.75)
abline(v=84,lty=2)
abline(v=88,lty=2)
abline(v=93,lty=2)
DUR1=DUR[1:84]
DUR1=DUR1-mean(DUR1)
DUR2=DUR[85:88]
DUR2=DUR2-mean(DUR2)
DUR3=DUR[89:93]
DUR3=DUR3-mean(DUR3)
DUR4=DUR[94:131]
DUR4=DUR4-mean(DUR4)
DUR.new=c(DUR1,DUR2,DUR3,DUR4)
## plot(DUR.new,type="l")
acf(DUR.new,xaxt="n",yaxt="n")
sd(DUR.new)

plot(NDUR,type="l",xaxt="n",yaxt="n")
legend("topleft",legend=c("NDUR "), cex = 0.75)
abline(v=72,lty=2)
NDUR1=NDUR[1:72]
NDUR1=NDUR1-mean(NDUR1)
NDUR2=NDUR[73:131]
NDUR2=NDUR2-mean(NDUR2)
NDUR.new=c(NDUR1,NDUR2)
## plot(NDUR.new,type="l")
acf(NDUR.new,xaxt="n",yaxt="n")
sd(NDUR.new)

plot(BUS,type="l",xaxt="n",yaxt="n")
legend("topleft",legend=c("BUS "), cex = 0.75)
acf(BUS,xaxt="n",yaxt="n")
plot(NFBUS,type="l",xaxt="n",yaxt="n")
legend("topleft",legend=c("NFBUS "), cex = 0.75)
acf(NFBUS,xaxt="n",yaxt="n")
plot(NFC,type="l",xaxt="n",yaxt="n")
legend("topleft",legend=c("NFC "), cex = 0.75)
acf(NFC,xaxt="n",yaxt="n")

dev.off()
plot.new()

################################################################################
#    Figure 1
################################################################################
utm1 <- function(k){
  aa <- matrix(0,k,k)
  for (i in 1:k){
    aa[i,i:k] <- 1
  }
  return(aa)
}

par(mfrow=c(1,2))

K <- 10
X <- cbind(rep(1,K),1:K)
U <- utm1(K)
lambda <- 0:640/400
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



plot(lambda/2,rs[1,]+2,type="l",ylim=yr,lty=3,xlab=expression(w),ylab=expression(r)) ## OLS-2L
legend("topleft", legend=c("OLS-2L","OLS-L", "GLS-2L","GLS-L"),
       lty=c(3,1,4,2), cex=0.8, x.intersp = 0.4)
lines(lambda/2,rs[2,]+2,type="l",lty=4) ## GLS-2L
lines(lambda/2,rs[3,]+2,type="l",lty=5) ## GLS-L
b1 <- (K+1)*(K+2)*(2*K+1)/15/K/(K-1)
b2 <- K/(K-1)
b3 <- (K+1)*(K+2)^2/6/K/(K-1)
lines(range(lambda/2),c(rs[1,1]+2,rs[1,1]+2+2*b3*max(lambda/2)),lty=1) ## OLS-L
abline(h=2) ## Oracle


K <- 15
X <- cbind(rep(1,K),1:K)
U <- utm1(K)
lambda <- 0:640/400
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



plot(lambda/2,rs[1,]+2,type="l",ylim=yr,lty=3,xlab=expression(w),ylab=expression(r)) ## OLS-2L
legend("topleft", legend=c("OLS-2L","OLS-L", "GLS-2L","GLS-L"),
       lty=c(3,1,4,2), cex=0.8, x.intersp = 0.4)
lines(lambda/2,rs[2,]+2,type="l",lty=4) ## GLS-2L
lines(lambda/2,rs[3,]+2,type="l",lty=5) ## GLS-L
b1 <- (K+1)*(K+2)*(2*K+1)/15/K/(K-1)
b2 <- K/(K-1)
b3 <- (K+1)*(K+2)^2/6/K/(K-1)
lines(range(lambda/2),c(rs[1,1]+2,rs[1,1]+2+2*b3*max(lambda/2)),lty=1) ## OLS-L
abline(h=2) ## Oracle
