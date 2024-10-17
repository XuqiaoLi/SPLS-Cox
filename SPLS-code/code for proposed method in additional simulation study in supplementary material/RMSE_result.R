rm(list=ls())
############### RMSE table###################
simutimes=500
n=200 # 500 or 1000
q=4
gindex=1 # 2 or 3
cvalue=0.1
load(paste('newPointwise','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))

h <- calh <- 0.01
p <- 1/h+1
t <- seq(0,1,calh)

#eigenfunction, FPC
eigenf1 <- function(x) sin(2*pi*x)*sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)*sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)*sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)*sqrt(2)

#beta(s)
betafun <- function(x){
  index=c(1,1/4,1/9,1/16)
  true_index=index/sqrt(sum(index^2))
  beta=true_index[1]*eigenf1(x)+true_index[2]*eigenf2(x)+
    true_index[3]*eigenf3(x)+true_index[4]*eigenf4(x)
}
truebeta<-betafun(t)


# linkfunction
g1 <- function(x) sin(2 * x) + 2 * cos(2 + x) - 2 * cos(2)
g2 <- function(x) x
g3 <- function(x) exp(x) - 1

# record valid result, exclude which does not converge or is singular
valid=rep(NA,length(res))
j=1
for(i in 1:length(res)){
  if(!is.character(res[[i]])){valid[j]=i;j=j+1}
}

valid=na.omit(valid)
if(length(valid)>simutimes){valid=valid[1:simutimes]}

# cindex loglik betahat gammahat result, here the gamma is the gamma in paper,but it is named by alpha in Rdata saved from simulation.R
cindex_fit=matrix(NA,nrow=simutimes,2)
beta_hat=matrix(0,nrow=simutimes,ncol=p)
gamma_hat=matrix(0,nrow=simutimes,ncol=q)

for (count in 1:length(valid)){
  res_count=res[[valid[count]]]
  cindex_fit[count,]=res_count$corindex
  beta_hat[count,]=t(res[[valid[count]]]$eigfun%*%res[[valid[count]]]$beta)
  gamma_hat[count,]=res[[valid[count]]]$alpha
}

best_res=which.max(cindex_fit[,1])
cindex_mean=round(mean(cindex_fit[,1]),3)
cindex_sd=round(sd(cindex_fit[,1]),3)

beta_true=betafun(t)
beta_true=as.matrix(t(beta_true))
beta_true_rep=matrix(rep(t(beta_true),simutimes),ncol=ncol(beta_true),byrow=TRUE)

denominator=matrix(0,simutimes)
numerator=matrix(0,simutimes)
RMSE_beta=matrix(0,simutimes)


for (i in 1:simutimes){
  denominator[i]=sum(beta_true_rep[i,]^2 * h)
  numerator[i]=sum((beta_hat[i,]-beta_true_rep[i,])^2*h)
  RMSE_beta[i]=numerator[i]/denominator[i]
}


RMSE_beta_mean=round(mean(RMSE_beta),3)
RMSE_beta_sd=round(sd(RMSE_beta),3)


gamma0=rep(0.2,q)
gamma0=as.matrix(t(gamma0))
gamma0_rep=matrix(rep(t(gamma0),simutimes),ncol=ncol(gamma0),byrow=TRUE)

RMSE_gamma=matrix(0,simutimes)
for (i in 1:simutimes){
  RMSE_gamma[i]=sum((gamma_hat[i,]-gamma0_rep[i,])^2)/sum(gamma0_rep[i,]^2)
}

RMSE_gamma_mean=round(mean(RMSE_gamma),3)
RMSE_gamma_sd=round(sd(RMSE_gamma),3)


### RMSE result and C-index
paste(RMSE_beta_mean,RMSE_beta_sd,RMSE_gamma_mean,RMSE_gamma_sd,cindex_mean,cindex_sd,sep=" & ")
