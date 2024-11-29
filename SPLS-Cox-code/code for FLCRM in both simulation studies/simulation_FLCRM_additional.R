###############################
rm(list=ls())
options(warn=-1)
library(locpol)
library(MASS)
library(survival)
library(survAUC)
library(survcomp)
library(fda)
library(splines)
library(statmod)
library(snowfall)
###############################

thisseed=2020
cumulative_var=0.95 
errorvariance= 0.5

g1 <- function(x) sin(2*x)+2*cos(2+x)-2*cos(2)
g2 <- function(x) x
g3 <- function(x) exp(x)-1

###################
n=200 
g=g1
gindex=1 
crindex=1 
censoringvalues<-c(0.1, 0.3) 
cvalue=censoringvalues[crindex] 

taulist1<-c(9.01,2.5) 
taulist2<-c(14.2,4.1) 
taulist3<-c(11.5,3.4) 

taulist=taulist1 

h <- calh <- 0.01
p <- 1/h+1 
t <- seq(0,1,calh)

lambda=c(1,1/2,1/4,1/8)

q<-4
mygamma<-rep(0.2,q) 

eigenf1 <- function(x) sin(2*pi*x)*sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)*sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)*sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)*sqrt(2)

betafun <- function(x){ 
  index=c(1,1/4,1/9,1/16)
  true_index=index/sqrt(sum(index^2))
  beta=true_index[1]*eigenf1(x)+true_index[2]*eigenf2(x)+
    true_index[3]*eigenf3(x)+true_index[4]*eigenf4(x)
}
truebeta<-betafun(t)

LengthAlpha=q

############ functions defined ##############
generator <- function(g,crindex,taulist) 
{
  Sigma<-matrix(0,4,4)
  for (j in 1:4) Sigma[j,j]=sqrt(lambda[j])
  
  SigmaZ<-0.5^t(sapply(1:q, function(i, j) abs(i-j), 1:q))
  
  correlationmatrix<-matrix(0,nrow(Sigma),q)

  correlationmatrix[1,1:q]<-0.1
  correlationmatrix[1:q,1]<-0.1
  
  Bigsigma<-rbind(cbind(Sigma, correlationmatrix), 
                  cbind(t(correlationmatrix), SigmaZ))
  
  mu<-rep(0,nrow(Bigsigma))
  data<-mvrnorm(n,mu,Bigsigma)
  
  B <- cbind(eigenf1(t), eigenf2(t), eigenf3(t), eigenf4(t))
  score <- matrix(0, nrow=n, ncol=length(lambda))
  
  for(j in 1:4)
    score[,j] <- data[,j]
  rawX <- t(matrix(rep(t,n),nrow=p,ncol=n)) + score %*% t(B)
  
  ZS=data[,(nrow(Sigma)+1):(nrow(Sigma)+q)]

  error<-matrix(rnorm(n*p,0,sqrt(errorvariance)),n,p)
  observeX<-rawX+error
  
  X<-matrix(NA,n,p)
  for (i in 1:n){
    dataframe <- data.frame(t)
    dataframe$observeX<-observeX[i,]
    lpfit <- locpol(observeX~t,dataframe, xeval=t)
    X[i,]<-lpfit$lpFit[,2]
  }
  
  scaleX<-scale(X,center=TRUE,scale=FALSE)
  SVDresult=svd(scaleX)
  Eigvec=SVDresult$v 
  Singval=SVDresult$d 
  est_Eigval=Singval^2*h/(n-1)
  
  for(i in 1:10){
    culVar=sum(est_Eigval[1:i])/sum(est_Eigval)
    if(culVar>=cumulative_var){VARselectindex=i;break}
  }
  
  tau=taulist[crindex]
  
  parameter<-exp(g((rawX-t(matrix(rep(t,n),nrow=p,ncol=n)))%*%truebeta*calh)+
                   ZS%*%mygamma)
  
  failuretime<-rexp(n, rate = parameter)
  censoringtime<-runif(n,0,tau)
  
  event<-rep(0,n)
  event<-as.numeric(failuretime<censoringtime)
  
  time<-failuretime*event+censoringtime*(rep(1,n)-event)
  
  dataframe<-list(time=time,event=event,scaleX=scaleX,ZS=ZS,Eigvec=Eigvec,VARselectindex)
}

paraFLCRM<-function(count){
  
  set.seed(thisseed+count)
  
  dataframe <- generator(g,crindex,taulist) #time,event,scaleX,ZS,Eigvec
  time=dataframe[1][[1]]
  event=dataframe[2][[1]]
  scaleX=dataframe[3][[1]]
  ZS=dataframe[4][[1]]
  Eigvec=dataframe[5][[1]]
  
  VARselectindex=dataframe[6][[1]] 
  
  
  est_Eigfun=Eigvec[,1:VARselectindex]/sqrt(h)
  
  eigenscore<-scaleX%*%est_Eigfun*h
  
  designmatrix<-cbind(eigenscore,ZS)
  
  newdataframe<-data.frame(cbind(time, event, designmatrix))
  
  totalvariable<-"V3"
  for (i in 4:(ncol(designmatrix)
               +2)){
    totalvariable<-paste(totalvariable,"+V", i, sep="")
  }
  
  totalformular<-as.formula(paste("Surv(time,event)~", totalvariable, sep=""))
  
  fullcoxresult<-coxph(formula=totalformular, data=newdataframe)
  
  cindex<-concordance.index(predict(fullcoxresult),surv.time = time, 
                            surv.event = event ,method = "noether")
  corindex<-c(cindex$c.index,cindex$se)
  
  return(list=c(corindex=corindex))
}

#############################
replicate_time=500

n=200 # 500 or 1000
g=g1 # g2 or g3
gindex=1 # 2 or 3

taulist=switch(gindex,taulist1,taulist2,taulist3)
crindex=1
cvalue=censoringvalues[crindex] 

sfInit(parallel = TRUE,cpus = 14)
sfLibrary(snowfall)
sfLibrary(locpol)
sfLibrary(MASS)
sfLibrary(survival)
sfLibrary(survAUC)
sfLibrary(survcomp)
sfLibrary(fda)
sfLibrary(splines)
sfLibrary(statmod)
sfExportAll()
res=sfClusterApplyLB(1:replicate_time,paraFLCRM)
sfStop()
save(res, file=paste('FLCRM','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))

#############################################
cindex_res=rep(NA,replicate_time)
for(i in 1:replicate_time){cindex_res[i]=res[[i]][1]}
cindex_mean=round(mean(cindex_res),3)
cindex_sd=round(sd(cindex_res),3)
c(cindex_mean,cindex_sd)