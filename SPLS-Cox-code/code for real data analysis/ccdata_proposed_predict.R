rm(list=ls())
options(warn=-1)
library(locpol)
library(MASS)
library(survival)
library(survAUC)
library(survcomp)
library(fda)
library(splines)
library(MASS)
library(statmod)
library(parallel)
library(snowfall)

load("clinical.dat") 
scoredata=scoredata[-331,] # clinical data
cleandata=read.csv('ccdata.csv',header = F) # CC data
CCinfo=read.csv('CCinfo.csv',header = T) # CC size data

observeX=as.matrix(cleandata)
ZS=scoredata[,4:15]
info=apply(CCinfo[,2:6]/CCinfo[,7],2,scale)

ZS=cbind(ZS,info)
time=scoredata[,1]
event=scoredata[,2]

n=372
n.test=122
n.train=n-n.test

p=200
q=17
t <- seq(0,2,length.out=200)
h=t[2]
###############################
thisseed=1000

crit <- 0.01
corr <- 1e-13  

nknots = 4  # number of knots
df = nknots+1  #cubic spline

###############data process####################
# smoothing
X<-matrix(NA,n,p)
for (i in 1:n){
  dataframe <- data.frame(t)
  dataframe$observeX<-observeX[i,]
  lpfit <- locpol(observeX~t,dataframe, xeval=t)
  X[i,]<-lpfit$lpFit[,2]
}

Eigvec=svd(scale(X,center=TRUE,scale=FALSE))$v # eigen vector
est_Eigval=svd(scale(X,center=TRUE,scale=FALSE))$d^2*h/(n-1) # eigen value
Rnindex=16 # cumulative variance 85%
############## Our method #######################
CC_estimate_predict<-function(init_seed){
  set.seed(thisseed+init_seed)
  index.train=sample(x = n,size = n.train,replace = F) #index in training set
  
  LengthAlpha=q
  LengthBeta=Rnindex
  
  est_Eigfun=Eigvec[,1:Rnindex]/sqrt(h)
  eigenscore<-scale(X,center=TRUE,scale=FALSE)%*%est_Eigfun*h  #FPC
  
  designmatrix<-cbind(eigenscore,ZS)
  newdataframe<-data.frame(cbind(time, event, designmatrix))
  newdataframe.train=newdataframe[index.train,]
  newdataframe.test=newdataframe[-index.train,]
  
  # following notation is named for utilizing the code of Sun et al.(2008) Polynomial spline estimation of partially linear single-index proportional hazards regression models
  data <- newdataframe.train
  x <- t(data[ ,3:(2+Rnindex)]) #FPC
  delta<-data[ ,2]
  v <- t(data[ , -1:-(2+Rnindex)]) #scalar Z
  Z<-data[,1]   #observed time
  
  tm<-sort(unique(Z[delta==1]))
  m<-length(tm) 
  n<-length(Z) 
  mat <- x
  
  count=0
  while(count<1){
    # initialized
    beta <- rnorm(Rnindex)
    alpha <- rnorm(q)
    beta <- beta/sqrt(sum(beta^2))
    gamma <- rnorm(df) 
    
    criterion <-5
    loop <- 0
    
    while(criterion > crit && loop < 50)
    {
      xbeta<-t(beta) %*% mat  # 1*n
      xbeta<-t(xbeta)  # n*1
      xbetas<-sort(xbeta)
      step_xbeta<- (range(xbetas)[2]-range(xbetas)[1]+corr)/(nknots-1)
      knotspos_xbeta<-rep(0,nknots)
      for(i in 1:nknots)  knotspos_xbeta[i]<- range(xbetas)[1]-corr/2+(i-1)*step_xbeta
      f <- function(x)  {  B <- bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0) }
      
      Btilda<-matrix(0,n,df)
      ab<- max(range(xbeta)[1],0)
      
      for(inte in 1:n)
      {
        aa<-xbeta[inte]
        for(intedf in 1:df)Btilda[inte,intedf]<-integrate(function(x)f(x)[,intedf],ab, aa)$value
      }
      
      B <- bsplineS(xbeta, breaks=knotspos_xbeta, norder=3, nderiv=0)
      Bder <- bsplineS(xbeta, breaks=knotspos_xbeta, norder=3, nderiv=1)
      
      alphaV <- t(v)%*% alpha 
      
      ######HessianBeta#########
      H1temp<-Bder %*% gamma
      H1 <-0
      H2 <-0
      for(i in 1:m)
      {
        indD <- as.numeric(Z==tm[i]) * delta # indicator for Di set
        Dsize <- sum(indD) 
        A<-diag(c(indD * H1temp)) # n*n diag matrix
        H1 <- H1+ mat %*% A %*% t(mat)
        
        indvec<-as.numeric(Z>=tm[i])  # 1*n indicator
        denom  <-sum(exp(Btilda %*% gamma+alphaV) * indvec)
        H2a <- 0
        H2b <- 0
        for(j in 1:n)    
        {
          ind<-1*(Z[j]>=tm[i])  # a number as indicator
          ex <- sum(exp(Btilda[j,] %*% gamma+alphaV[j]))       
          H2a <- H2a + ind * sum(Bder[i,] %*% gamma+ (B[j,] %*% gamma)**2) * ex * (mat[,j] %*% t(mat[,j])) 
          #### cuz its 1*1 matrix instead of a number, so use "sum" to make the matrix into a number, without changing anything ###
          H2b <- H2b + ind * ex * (B[j,] %*% gamma) * mat[,j]  
        }
        
        H21 <- H2a/denom
        H22 <- H2b %*% t(H2b)/denom**2
        H2 <- H2+ (H21-H22)*Dsize
      }
      H <- H1-H2
      ##########################
      HB <- H
      ######ScoreBeta#########
      s1p <- matvec(mat, (B %*% gamma))
      s1<- 0   ## first part
      ex<-exp(Btilda %*% gamma+alphaV)* (B %*% gamma) # n*1 
      numer <- matvec(mat,ex) 
      s2<-0
      for(i in 1:m)
      {
        indD <- as.numeric(Z==tm[i]) * delta # indicator for Di set
        Dsize <- sum(indD) 
        s1 <- s1 + matvec(s1p,indD) %*% rep(1,n)
        
        ind<-as.numeric(Z>=tm[i])  # 1*n indicator
        newnum <- matvec(numer, ind) %*% rep(1,n) 
        denom  <- sum(exp(Btilda %*% gamma+alphaV) * ind)
        s2<-s2 + Dsize*newnum/denom
      }
      s <- s1-s2
      ########################
      SB <- s
      
      inverseHB=try(solve(HB),TRUE)
      if(is.character(inverseHB)) {cat('Hessian Beta singular'); break}
      betanew <- beta - inverseHB %*%SB
      
      if (is.na(betanew[1])) {cat('is.na(betanew[1])'); break}
      if(betanew[1] <= 0) betanew <- -betanew
      betanew <- betanew/sqrt(sum(betanew**2))  # beta updated
      
      xbetanew<-t(betanew) %*% mat  # 1*n
      xbetanew<-t(xbetanew)  # n*1
      xbetanews<-sort(xbetanew)
      step_xbetanew<- (range(xbetanews)[2]-range(xbetanews)[1]+corr)/(nknots-1)
      knotspos_xbetanew<-rep(0,nknots)
      for(i in 1:nknots)  knotspos_xbetanew[i]<- range(xbetanews)[1]-corr/2+(i-1)*step_xbetanew
      fnew <- function(x) {  B <- bsplineS(x, breaks=knotspos_xbetanew, norder=3, nderiv=0) } 
      Btildanew<-matrix(0,n,df)
      ac<- max(range(xbetanew)[1],0)
      
      for(inte in 1:n)
      {
        aa<-xbetanew[inte]
        for(intedf in 1:df) Btildanew[inte,intedf]<-integrate(function(x)fnew(x)[,intedf],ac, aa)$value
      }
      
      ss <- try(coxph(Surv(Z,delta) ~ Btildanew + t(v),  
                      method="breslow", eps=1e-7, iter.max=20), TRUE)
      
      if(is.character(ss)) {cat('coxph something wrong'); break}
      coeff<-coef(ss)
      aa<-is.na(coeff)  
      if(any(aa)) coeff[aa] <-0
      
      gammanew <- coeff[1:df]
      alphanew <- coeff[(df+1): (LengthAlpha+df)]
      
      lastnorm2=sum((beta)**2)+sum((gamma)**2)+sum((alpha)**2)
      Thiserror2=sum((betanew-beta)**2)+sum((gammanew-gamma)**2)
      +sum((alphanew-alpha)**2)
      criterion <- sqrt(Thiserror2/lastnorm2)
      
      
      cat(count, "th simu: ","convergance criterion is: ", criterion,  
          "difference is:",criterion-crit," number of loops: ", loop,"\n")
      
      beta <- betanew
      gamma <- gammanew
      alpha<-alphanew
      loop <- loop+1
      if (is.na(criterion)) break     
    }  
    
    if(criterion<=crit) { 
      count=count+1
      #cindex of fitting
      cindex<-concordance.index(predict(ss),surv.time = Z, 
                                surv.event = delta ,method = "noether")
      corindex.train<-c(cindex$c.index,cindex$se)
      
      
      #####################################################################################
      
      #use the same code to rename test data
      data <- newdataframe.test
      x <- t(data[ ,3:(2+Rnindex)]) #FPC
      delta<-data[ ,2]
      v <- t(data[ , -1:-(2+Rnindex)]) #scalar Z
      Z<-data[,1]   #observed time
      
      tm<-sort(unique(Z[delta==1]))
      m<-length(tm) 
      n<-length(Z) 
      mat <- x
      
      # The previous code regarding score and hessian beta is not used. We only need the result from cox model and test data.
      # Use the test data to build up Btilda 
      xbetanew<-t(betanew) %*% mat  # 1*n
      xbetanew<-t(xbetanew)  # n*1
      xbetanews<-sort(xbetanew)
      step_xbetanew<- (range(xbetanews)[2]-range(xbetanews)[1]+corr)/(nknots-1)
      knotspos_xbetanew<-rep(0,nknots)
      for(i in 1:nknots)  knotspos_xbetanew[i]<- range(xbetanews)[1]-corr/2+(i-1)*step_xbetanew
      fnew <- function(x) {  B <- bsplineS(x, breaks=knotspos_xbetanew, norder=3, nderiv=0) } 
      Btildanew<-matrix(0,n,df)
      ac<- max(range(xbetanew)[1],0)
      
      for(inte in 1:n)
      {
        aa<-xbetanew[inte]
        for(intedf in 1:df) Btildanew[inte,intedf]<-integrate(function(x)fnew(x)[,intedf],ac, aa)$value
      }
      
      pre.new=predict(ss,newdata = data.frame(Btildanew,t(v)),type = 'lp') # use test data for linear prediction
      cindex.pre<-concordance.index(pre.new,surv.time = newdataframe.test$time, 
                                    surv.event = newdataframe.test$event ,method = "noether")
      corindex.test=c(cindex.pre$c.index,cindex.pre$se)
      
      ###########################################################
      return(list(corindex.train=corindex.train,corindex.test=corindex.test))
    }
  }
}



replicate_time=1000
sfInit(parallel = TRUE,cpus = 78)
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
res=sfClusterApplyLB(1:replicate_time,CC_estimate)
sfStop()
save(res,file=paste('realdata_',Rnindex,'CCfpca_Prediction2026_test',n.test,'.Rdata',sep = ''))

