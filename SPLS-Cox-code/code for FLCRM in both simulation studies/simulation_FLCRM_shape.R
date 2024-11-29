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
library(mvtnorm)
###############################
thisseed=2024
cleandata=read.csv('ccdata.csv',header = F)
observeX=as.matrix(cleandata)

n=372
p=200
q=6
h <- calh <-1  
t <- seq(1,ncol(cleandata),h)
mygamma<-rep(0.2,q)
cumulative_var=0.85
errorvariance <- 0.001

g1 <- function(x) sin(2*x)+2*cos(2+x)-2*cos(2)
g2 <- function(x) x
g3 <- function(x) exp(x)-1

###################
g=g1
gindex=1 
crindex=1 
censoringvalues<-c(0.1)
cvalue=censoringvalues[crindex] 

taulist1<-c(11.3) 
taulist2<-c(18.7) 
taulist3<-c(14.3) 

taulist=switch(gindex,taulist1,taulist2,taulist3) 

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

############ functions defined ##############
generator <- function(g,crindex,taulist) 
{ 
  #Use the covariance matrix of first six FPCs and six-dimension Zi (12*12) to generate this part of data.
  #Then the remain FPCs are independently generated from N(0,λj) using rnorm
  
  SigmaZ<-0.5^t(sapply(1:q, function(i, j) abs(i-j), 1:q))
  Sigma<-matrix(0,6,6)
  for (j in 1:6) Sigma[j,j]=shape_est_Eigval[j] #estimated eigenvalue
  correlationmatrix<-matrix(0,nrow(Sigma),q)
  correlationmatrix[1,1:q]<-0.1
  correlationmatrix[1:q,1]<-0.1
  
  data<-matrix(0,nrow = n,ncol = shape_est_Eigfun_number+q)
  
  Bigsigma<-rbind(cbind(Sigma, correlationmatrix),
                  cbind(t(correlationmatrix), SigmaZ))
  
  #first generate all the FPCs
  for(i in 1:shape_est_Eigfun_number){
    data[,i]=rnorm(n,mean = 0,sd=sqrt(shape_est_Eigval[i]))
  }
  
  #them generate the Zi by conditional distribution
  invSigma=solve(Sigma)
  ConditionalZ_Sigma= SigmaZ-correlationmatrix%*%invSigma%*%t(correlationmatrix)
  for(i in 1:n){
    data[i,(shape_est_Eigfun_number+1):(shape_est_Eigfun_number+q)]=
      rmvnorm(1,mean = correlationmatrix%*%invSigma%*%matrix(data[i,1:6],ncol = 1),sigma = ConditionalZ_Sigma)
  }
  
  #The generated FPC scores
  score <- matrix(0, nrow=n, ncol=shape_est_Eigfun_number)
  for(j in 1:shape_est_Eigfun_number)
    score[,j] <- data[,j]
  
  rawX <- t(matrix(rep(mu_s,n),nrow=p,ncol=n)) + score %*% t(shape_est_Eigfun)
  
  ZS=data[,(shape_est_Eigfun_number+1):(shape_est_Eigfun_number+q)]
  
  error<-matrix(rnorm(n*p,0,sqrt(errorvariance)),n,p)
  observeX<-rawX+error #W=X+epsilon
  
  tau=taulist[crindex]
  
  parameter<-exp(g((rawX-t(matrix(rep(mu_s,n),nrow=p,ncol=n)))%*%truebeta*calh)+
                   ZS%*%mygamma)
  
  failuretime<-rexp(n, rate = parameter)
  censoringtime<-runif(n,0,tau)
  
  event<-rep(0,n)
  event<-as.numeric(failuretime<censoringtime)
  sum(event)/n
  
  time<-failuretime*event+censoringtime*(rep(1,n)-event)
  
  #smooth Xi(s)
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
  
  for(i in 1:30){
    culVar=sum(est_Eigval[1:i])/sum(est_Eigval)
    if(culVar>=cumulative_var){VARselectindex=i;break}
  }
  
  dataframe<-list(time=time,event=event,scaleX=scaleX,ZS=ZS,
                  Eigvec=Eigvec,VARselectindex=VARselectindex)
}

paraFLCRM<-function(count){
  
  set.seed(thisseed+count)
  
  dataframe <- generator(g,crindex,taulist) #time,event,scaleX,ZS,Eigvec
  time=dataframe$time
  event=dataframe$event
  scaleX=dataframe$scaleX
  ZS=dataframe$ZS
  Eigvec=dataframe$Eigvec
  
  VARselectindex=dataframe$VARselectindex
  
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
sfLibrary(mvtnorm)
sfExportAll()
res=sfClusterApplyLB(1:replicate_time,paraFLCRM)
sfStop()
save(res, file=paste('shape_FLCRM','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))

#############################################
cindex_res=rep(NA,replicate_time)
for(i in 1:replicate_time){cindex_res[i]=res[[i]][1]}
cindex_mean=round(mean(cindex_res),3)
cindex_sd=round(sd(cindex_res),3)
c(cindex_mean,cindex_sd)