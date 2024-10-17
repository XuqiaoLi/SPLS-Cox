rm(list=ls())

cleandata=read.csv('ccdata.csv',header = F)
observeX=as.matrix(cleandata)
############### RMSE table ###################
simutimes=500
gindex=1 # 2 or 3
cvalue=0.1

n=372
p=200
q=6
h <- calh <-1 
t <- seq(1,ncol(cleandata),h)
mygamma<-rep(0.2,q)
LengthAlpha=q
load(paste('newShape','g',gindex,'cr',cvalue,'.Rdata',sep = ''))

# linkfunction
g1 <- function(x) sin(2 * x) + 2 * cos(2 + x) - 2 * cos(2)
g2 <- function(x) x
g3 <- function(x) exp(x) - 1

### CC data for FPCA
X<-matrix(NA,n,p)
for (i in 1:n){
  dataframe <- data.frame(t)
  dataframe$observeX<-observeX[i,]
  lpfit <- locpol(observeX~t,dataframe, xeval=t)
  X[i,]<-lpfit$lpFit[,2]
}

scale_shape<-scale(X,center=TRUE,scale=FALSE) #-μ(s)
mu_s=attributes(scale_shape)[["scaled:center"]] #μ(s)
shape_SVDresult=svd(scale_shape)
shape_Eigvec=shape_SVDresult$v 
shape_Singval=shape_SVDresult$d 
shape_est_Eigval=shape_Singval^2*h/(n-1) 
shape_est_Eigfun=shape_Eigvec/sqrt(h)
shape_est_Eigfun_number=length(shape_est_Eigval)

# Only take the first 9 eigen functions that contribute to total 73% variance, to generate data
shape_est_Eigfun_number=9
shape_est_Eigfun=shape_est_Eigfun[,1:shape_est_Eigfun_number]

#beta coefficient
beta_coef=rep(1,shape_est_Eigfun_number)/(seq(1,shape_est_Eigfun_number))^2
beta_coef=beta_coef/sqrt(sum(beta_coef**2)) 

truebeta <- c(beta_coef %*% t(shape_est_Eigfun))

remove(observeX,X,cleandata)

# record valid result, exclude which does not converge or is singular
valid=rep(NA,length(res))
j=1
for(i in 1:length(res)){
  if(!is.character(res[[i]])){valid[j]=i;j=j+1}
}

valid=na.omit(valid)
if(length(valid)>simutimes){valid=valid[1:simutimes]}



# cindex loglik betahat gammahat result, here the gamma is the gamma in paper,but it is named by alpha in Rdata saved from simulation_shape.R
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

beta_true=as.matrix(t(truebeta))
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
