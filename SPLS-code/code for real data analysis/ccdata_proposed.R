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

############## Our method #######################
CC_estimate<-function(init_seed){
  LengthAlpha=q
  LengthBeta=Rnindex
  
  est_Eigfun=Eigvec[,1:Rnindex]/sqrt(h)
  eigenscore<-scale(X,center=TRUE,scale=FALSE)%*%est_Eigfun*h  #FPC
  
  designmatrix<-cbind(eigenscore,ZS)
  newdataframe<-data.frame(cbind(time, event, designmatrix))
  
  # following notation is named for utilizing the code of Sun et al.(2008) Polynomial spline estimation of partially linear single-index proportional hazards regression models
  data <- newdataframe 
  x <- t(data[ ,3:(2+Rnindex)]) #FPC
  delta<-data[ ,2]
  v <- t(data[ , -1:-(2+Rnindex)]) #scalar Z
  Z<-data[,1]   #observed time
  
  tm<-sort(unique(Z[delta==1]))
  m<-length(tm) 
  n<-length(Z) 
  mat <- x
  
  set.seed(thisseed+init_seed)
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
      B=Bnew <- bsplineS(xbetanew, breaks=knotspos_xbetanew, norder=3, nderiv=0)
      Bder=Bdernew <- bsplineS(xbetanew, breaks=knotspos_xbetanew, norder=3, nderiv=1)
      
      Btilda=Btildanew
      alphaV=t(v)%*% alpha 
      
      ######HessianBeta#########
      H1temp<-Bder %*% gamma # n*1
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
      HB <- H
      ##########################################
      invHB = solve(HB)
      
      ########### HHH #########
      temp1<- (betanew/betanew[1])[2:LengthBeta]
      ll<-LengthBeta-1+df+LengthAlpha
      temp2<- rep(0,df+LengthAlpha)
      temp3<- diag(rep(1,(ll)))
      Gt<-matrix(c(temp1,temp2,temp3), ll,(ll+1)) 
      G<-t(Gt)
      
      HHH <- function(beta, gamma, alpha, B, Bder, Btilda)
      {
        alphaV <- t(v)%*% alpha 
        sigma<-beta[2:LengthBeta]
        norm <- 1 - sum(sigma**2)
        epsilon <- -mat[1,] %*% t(sigma)/sqrt(norm) + t(mat[2:LengthBeta, ])
        matA <- sigma%*%t(sigma)/norm+diag(rep(1, LengthBeta-1))
        H1temp <- B %*% gamma
        H2temp<-Bder %*% gamma 
        H1 <-0
        H2 <-0
        H3 <-0
        H4 <-0
        gH1 <-0
        gH2 <-0
        aH1 <-0
        aH2 <-0
        gaH1<-0
        gaH2<-0
        
        Heg1 <-0
        Heg2<-0
        Heg3<-0
        
        Hea1<-0
        Hea2<-0
        
        for(i in 1:m)
        {
          indD <- as.numeric(Z==tm[i]) * delta # indicator for Di set
          Dsize <- sum(indD) 
          
          H1<- H1+sum(indD * (-mat[1,]/sqrt(norm)) * H1temp) * matA
          
          A2<-diag(c(indD * H2temp)) # n*n diag matrix
          H2 <- H2+ t(epsilon) %*% A2 %*% epsilon
          
          Aeg1<-diag(indD) # n*n diag matrix
          Heg1 <- Heg1+ t(epsilon) %*% Aeg1 %*% B
          
          indvec<-as.numeric(Z>=tm[i])  # 1*n indicator
          denom  <-sum(exp(Btilda %*% gamma+alphaV) * indvec)
          
          H3a <- 0
          H3b <- 0
          H4a<-0  
          
          gH1a<-0
          gH2a <- 0
          aH1a<-0
          aH2a <- 0
          gaH1a <-0
          
          egH2<-0
          egH3a<-0
          egH3b<-0
          
          eaH1<-0
          eaH2a<-0
          eaH2b<-0
          
          for(j in 1:n)    
          {
            ind<-1*(Z[j]>=tm[i])  # a number as indicator for Ri set
            ex <- sum(exp(Btilda[j,] %*% gamma+alphaV[j]))       
            
            H3a <- H3a + ind * ((-mat[1,j]/sqrt(norm)) * H1temp[j,]) * ex * matA
            H3b <- H3b + ind * c(Bder[i,] %*% gamma+ (B[j,] %*% gamma)**2) * ex * (epsilon[j,] %*% t(epsilon[j,])) 
            H4a <- H4a + ind * ex * (B[j,] %*% gamma) * epsilon[j,]
            
            gH1a <- gH1a + ind * ex * (Btilda[j,] %*% t(Btilda[j,])) 
            gH2a <- gH2a + ind * ex * Btilda[j,]
            aH1a <- aH1a + ind * ex * (v[,j] %*% t(v[,j])) 
            aH2a <- aH2a + ind * ex * v[,j]
            gaH1a <- gaH1a+ ind * ex * (Btilda[j,] %*% t(v[,j]))
            
            egH2 <- egH2 + ind * ex * ( epsilon[j,] %*% t(B[j,]) + c(B[j,]%*%gamma) * (epsilon[j,] %*% t(Btilda[j,])))         
            egH3a <- egH3a + ind * ex * c(B[j,] %*% gamma) * epsilon[j,]
            egH3b <- egH3b + ind * ex * Btilda[j,]
            
            eaH1 <- eaH1 + ind * ex * c(B[j,]%*%gamma) * (epsilon[j,] %*% t(v[,j]))
            eaH2a <- eaH2a + ind * ex * c(B[j,]%*%gamma) * epsilon[j,]
            eaH2b <- eaH2b + ind * ex * v[,j]
          }
          
          H3 <- H3 + (H3a+H3b)*Dsize/denom
          H4 <- H4 + Dsize* H4a %*% t(H4a)/denom**2
          
          gH1 <- gH1+Dsize*gH1a/denom
          gH2 <- gH2+ Dsize*gH2a %*% t(gH2a)/(denom**2)
          aH1 <- aH1+Dsize*aH1a/denom
          aH2 <- aH2+ Dsize*aH2a %*% t(aH2a)/(denom**2)
          gaH1 <- gaH1 + Dsize*gaH1a/denom
          gaH2 <- gaH2 + Dsize*gH2a %*% t(aH2a)/(denom**2)
          
          Heg2 <- Heg2+ Dsize*egH2/denom
          Heg3 <- Heg3 + Dsize* egH3a %*% t(egH3b)/(denom**2)
          
          Hea1 <- Hea1 + Dsize * eaH1/denom
          Hea2 <- Hea2 + Dsize * eaH2a %*% t(eaH2b)/(denom**2)
        }
        Hee <- H1+H2-H3+H4
        Hgg <- -gH1+gH2
        Haa <- -aH1+aH2
        Hga <- -gaH1+gaH2
        Heg <- Heg1 - Heg2 +Heg3
        Hea <- -Hea1 + Hea2
        
        m1 <- matrix(c(Hee,Heg,Hea),(LengthBeta-1),ll)
        m2 <- matrix(c(t(Heg),Hgg,Hga),df,ll) 
        m3 <- matrix(c(t(Hea),t(Hga),Haa),LengthAlpha,ll) 
        
        H <- matrix(c(t(m1),t(m2),t(m3)),ll,ll)
        H <- t(H)
        H
      }
      Heth <- HHH(beta, gamma, alpha, Bnew, Bdernew, Btildanew)
      
      midmat <- try(-solve(Heth),TRUE)
      if(is.character(midmat)) return('is.character(midmat)') 
      varmat<-diag(G %*% midmat %*% Gt)
      if(any(is.na(sqrt(varmat)))) return('is.na(sqrt(varmat)')  #get rid of NA in SEs;
      varmat<-(G %*% midmat %*% Gt)
      varbeta<-varmat[1:LengthBeta,1:LengthBeta]
      
      #############ScoreBeta i###############
      scorebeta_i=matrix(NA,Rnindex,m)
      
      s1p <- matvec(mat, (B %*% gamma))
      s1<- 0   ## first part
      ex<-exp(Btilda %*% gamma+alphaV)* (B %*% gamma) # n*1 
      numer <- matvec(mat,ex) 
      s2<-0
      for(i in 1:m){
        indD <- as.numeric(Z==tm[i]) * delta # indicator for Di set
        Dsize <- sum(indD) 
        s1_i= matvec(s1p,indD) %*% rep(1,n)
        
        ind<-as.numeric(Z>=tm[i])  # 1*n indicator
        newnum <- matvec(numer, ind) %*% rep(1,n) # to sum, get a 8*1 vector
        denom  <- sum(exp(Btilda %*% gamma+alphaV) * ind)
        s2_i= Dsize*newnum/denom
        
        scorebeta_i[,i]=s1_i-s2_i
      }
      
      ############################################################
      linkxbeta <- t(eigenscore %*% as.matrix(beta))
      linkres <- t(Btildanew %*% as.matrix(gamma))
      
      cindex<-concordance.index(predict(ss),surv.time = Z, 
                                surv.event = delta ,method = "noether")
      corindex<-c(cindex$c.index,cindex$se)
      cr=1-mean(delta)
      
      ###########################################################
      return(list(beta=beta,gamma=gamma,alpha=alpha,inv_hessianbt=invHB,varbeta=varbeta,
                  scorebeta_i=scorebeta_i,eigfun=est_Eigfun,summary=summary(ss),
                  num_distin=m,rn=Rnindex,cr=cr,loglik=summary(ss)$loglik[2],
                  linkxbeta=linkxbeta,linkres=linkres,corindex=corindex))
    }
  }
}


Rnindex=16 # cumulative variance 85%
replicate_time=6000
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
save(res,file=paste('realdata_',Rnindex,'CCfpca','.Rdata',sep = ''))

