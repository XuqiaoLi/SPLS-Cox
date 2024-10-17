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

crit = 0.01 
corr = 1e-13

nknots = 4  # number of knots
df = nknots+1  #cubic spline

#link function
g1 <- function(x) sin(2*x)+2*cos(2+x)-2*cos(2)
g2 <- function(x) x
g3 <- function(x) exp(x)-1

###################
n=1000 
g=g1
gindex=1 
crindex=1 
censoringvalues<-c(0.1, 0.3)
cvalue=censoringvalues[crindex] 

taulist1<-c(9.01,2.5) 
taulist2<-c(14.2,4.1) 
taulist3<-c(11.5,3.4) 

taulist=taulist1 # generate censor variable for specified censoring rate

h <- calh <- 0.01
p <- 1/h+1 
t <- seq(0,1,calh)

lambda=c(1,1/2,1/4,1/8)

q<-4 
mygamma<-rep(0.2,q) 

#eigenfunction
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

LengthAlpha=q

############ functions defined ##############
generator <- function(g,crindex,taulist) 
{
  #Xi(s) sigma
  Sigma<-matrix(0,4,4)
  for (j in 1:4) Sigma[j,j]=sqrt(lambda[j])
  
  #Z sigma
  SigmaZ<-0.5^t(sapply(1:q, function(i, j) abs(i-j), 1:q))
  
  correlationmatrix<-matrix(0,nrow(Sigma),q)
  #allow some correlation between Zi and Xi(s)
  correlationmatrix[1,1:q]<-0.1
  correlationmatrix[1:q,1]<-0.1
  
  Bigsigma<-rbind(cbind(Sigma, correlationmatrix), 
                  cbind(t(correlationmatrix), SigmaZ))
  
  mu<-rep(0,nrow(Bigsigma)) #mean zero
  data<-mvrnorm(n,mu,Bigsigma)
  
  # generate Xi(s)
  B <- cbind(eigenf1(t), eigenf2(t), eigenf3(t), eigenf4(t))
  score <- matrix(0, nrow=n, ncol=length(lambda))
  
  for(j in 1:4)
    score[,j] <- data[,j]
  rawX <- t(matrix(rep(t,n),nrow=p,ncol=n)) + score %*% t(B)
  
  # generateZ
  ZS=data[,(nrow(Sigma)+1):(nrow(Sigma)+q)]
  
  # add noisy term
  error<-matrix(rnorm(n*p,0,sqrt(errorvariance)),n,p)
  observeX<-rawX+error #W=X+epsilon
  
  # smooth Xi(s)
  X<-matrix(NA,n,p)
  for (i in 1:n){
    dataframe <- data.frame(t)
    dataframe$observeX<-observeX[i,]
    lpfit <- locpol(observeX~t,dataframe, xeval=t)
    X[i,]<-lpfit$lpFit[,2]
  }
  
  # Estimate FPC
  scaleX<-scale(X,center=TRUE,scale=FALSE) #-μ(s)
  SVDresult=svd(scaleX)
  Eigvec=SVDresult$v #eigenvector
  Singval=SVDresult$d 
  est_Eigval=Singval^2*h/(n-1) #eigenvalue
  
  #cumulative variance
  for(i in 1:10){
    culVar=sum(est_Eigval[1:i])/sum(est_Eigval)
    if(culVar>=cumulative_var){VARselectindex=i;break}
  }
  
  tau=taulist[crindex]
  
  #generate time-to-event data
  parameter<-exp(g((rawX-t(matrix(rep(t,n),nrow=p,ncol=n)))%*%truebeta*calh)+
                   ZS%*%mygamma)
  
  failuretime<-rexp(n, rate = parameter)
  censoringtime<-runif(n,0,tau)
  
  event<-rep(0,n)
  event<-as.numeric(failuretime<censoringtime)
  
  time<-failuretime*event+censoringtime*(rep(1,n)-event)
  
  dataframe<-list(time=time,event=event,scaleX=scaleX,ZS=ZS,Eigvec=Eigvec,VARselectindex=VARselectindex)
}

parasimulation<-function(count,init_seed=1){
  #sink(paste('PWCB',count,'.txt',sep = ''))
  set.seed(thisseed+count)
  dataframe <- generator(g,crindex,taulist) #time,event,scaleX,ZS,Eigvec
  time=dataframe$time
  event=dataframe$event
  scaleX=dataframe$scaleX
  ZS=dataframe$ZS
  Eigvec=dataframe$Eigvec
  
  VARselectindex=dataframe$VARselectindex
  
  LengthBeta=VARselectindex
  
  est_Eigfun=Eigvec[,1:VARselectindex]/sqrt(h) # eigen function
  beta0 <- as.vector(truebeta%*%(est_Eigfun)*h) # determine the direction such that beta1>0
  if (beta0[1]<0) est_Eigfun <- -est_Eigfun
  
  eigenscore<-scaleX%*%est_Eigfun*h 
  designmatrix<-cbind(eigenscore,ZS)
  data<-data.frame(cbind(time, event, designmatrix))
  
  # following notation is named for utilizing the code of Sun et al.(2008) Polynomial spline estimation of partially linear single-index proportional hazards regression models
  mat <- t(data[ ,3:(2+VARselectindex)]) #FPC
  delta<-data[ ,2]
  v <- t(data[ , -1:-(2+VARselectindex)]) #scalar Z
  Z<-data[,1]   #observed time
  
  tm<-sort(unique(Z[delta==1]))
  m<-length(tm) # distinctive event time
  n<-length(Z) # sample size
  
  # initial
  set.seed(thisseed+count+init_seed)
  beta = rnorm(VARselectindex)
  beta = beta/sqrt(sum(beta^2)) # norm
  if(beta[1]<0) beta=-beta
  betanew=beta
  alpha = rnorm(q)
  gamma = rnorm(df)
  
  criterion=5
  loop=0
  
  while(criterion > crit && loop < 50)
  {
    ######## Beta update #######
    xbeta<-t(beta) %*% mat
    xbeta<-t(xbeta) 
    xbetas<-sort(xbeta)
    step_xbeta<- (range(xbetas)[2]-range(xbetas)[1]+corr)/(nknots-1)
    knotspos_xbeta<-rep(0,nknots)
    for(i in 1:nknots)  knotspos_xbeta[i]<- range(xbetas)[1]-corr/2+(i-1)*step_xbeta
    f1 <- function(x) {  bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0)[,1]  }
    f2 <- function(x) {  bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0)[,2]  }
    f3 <- function(x) {  bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0)[,3]  }
    f4 <- function(x) {  bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0)[,4]  }
    f5 <- function(x) {  bsplineS(x, breaks=knotspos_xbeta, norder=3, nderiv=0)[,5]  }
    
    Btilda<-matrix(0,n,df)
    ab<- max(range(xbeta)[1],0)
    
    f1integrate<-function(x) return(integrate(f1, ab, x)$value)
    f2integrate<-function(x) return(integrate(f2, ab, x)$value)
    f3integrate<-function(x) return(integrate(f3, ab, x)$value)
    f4integrate<-function(x) return(integrate(f4, ab, x)$value)
    f5integrate<-function(x) return(integrate(f5, ab, x)$value)
    
    Btilda[,1]=apply(xbeta,MARGIN = 1,f1integrate)
    Btilda[,2]=apply(xbeta,MARGIN = 1,f2integrate)
    Btilda[,3]=apply(xbeta,MARGIN = 1,f3integrate)
    Btilda[,4]=apply(xbeta,MARGIN = 1,f4integrate)
    Btilda[,5]=apply(xbeta,MARGIN = 1,f5integrate)
    
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
    HB <- H
    
    #############ScoreBeta###############
    s1p <- matvec(mat, (B %*% gamma)) 
    s1<- 0
    ex<-exp(Btilda %*% gamma+alphaV)* (B %*% gamma)  
    numer <- matvec(mat,ex)  
    s2<-0
    for(i in 1:m){
      indD <- as.numeric(Z==tm[i]) * delta # indicator for Di set
      Dsize <- sum(indD) 
      s1 <- s1 + matvec(s1p,indD) %*% rep(1,n)
      
      ind<-as.numeric(Z>=tm[i])  # 1*n indicator
      newnum <- matvec(numer, ind) %*% rep(1,n) 
      denom  <- sum(exp(Btilda %*% gamma+alphaV) * ind)
      s2<-s2 + Dsize*newnum/denom
    }
    s <- s1-s2
    SB <- s
    
    invHB = solve(HB)
    
    betanew = beta -invHB%*% SB 
    betanew = betanew/sqrt(sum(betanew**2))
    
    if(is.na(betanew[1])) break
    if(betanew[1] <= 0) betanew <- -betanew
    
    ######## gamma alpha update#######
    xbetanew<-t(betanew) %*% mat  # 1*n
    xbetanew<-t(xbetanew)  # n*1
    xbetanews<-sort(xbetanew)
    step_xbetanew<- (range(xbetanews)[2]-range(xbetanews)[1]+corr)/(nknots-1)
    knotspos_xbeta<-rep(0,nknots)
    for(i in 1:nknots)  knotspos_xbeta[i]<- range(xbetanews)[1]-corr/2+(i-1)*step_xbetanew
    
    Btildanew<-matrix(0,n,df)
    ab<- max(range(xbetanew)[1],0)
    
    Btildanew[,1]=apply(xbetanew,MARGIN = 1,f1integrate)
    Btildanew[,2]=apply(xbetanew,MARGIN = 1,f2integrate)
    Btildanew[,3]=apply(xbetanew,MARGIN = 1,f3integrate)
    Btildanew[,4]=apply(xbetanew,MARGIN = 1,f4integrate)
    Btildanew[,5]=apply(xbetanew,MARGIN = 1,f5integrate)
    
    ss <- try(coxph(Surv(Z,delta) ~ Btildanew + t(v),  
                    method="breslow", eps=1e-7, iter.max=20), TRUE)
    
    if(is.character(ss)) break
    coeff<-coef(ss)
    aa<-is.na(coeff)  
    if(any(aa)) coeff[aa] <-0
    
    gammanew <- coeff[1:df]
    alphanew <- coeff[(df+1): (LengthAlpha+df)]
    
    lastnorm2=sum((beta)**2)+sum((gamma)**2)+sum((alpha)**2)
    Thiserror2=sum((betanew-beta)**2)+sum((gammanew-gamma)**2)+sum((alphanew-alpha)**2)
    criterion <- sqrt(Thiserror2/lastnorm2)
    
    
    cat(count, "th simu: ","convergance criterion is: ", criterion,  
        "difference is:",criterion-crit," number of loops: ", loop,"\n")
    
    beta = betanew
    gamma = gammanew
    alpha = alphanew
    loop = loop+1
    if (is.na(criterion)) break     
  }  
  
  if(criterion<=crit) {
    # sink()
    B=Bnew <- bsplineS(xbetanew, breaks=knotspos_xbeta, norder=3, nderiv=0)
    Bder=Bdernew <- bsplineS(xbetanew, breaks=knotspos_xbeta, norder=3, nderiv=1)
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
    
    #############ScoreBeta_i###############
    scorebeta_i=matrix(NA,VARselectindex,m)
    
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
      newnum <- matvec(numer, ind) %*% rep(1,n)
      denom  <- sum(exp(Btilda %*% gamma+alphaV) * ind)
      s2_i= Dsize*newnum/denom
      
      scorebeta_i[,i]=s1_i-s2_i
    }
    
    ###############
    linkxbeta <- t(eigenscore %*% as.matrix(beta))
    linkres <- t(Btildanew %*% as.matrix(gamma)) 
    
    cindex<-concordance.index(predict(ss),surv.time = Z, 
                              surv.event = delta ,method = "noether")
    corindex<-c(cindex$c.index,cindex$se)
    cr=1-mean(delta)

    ##############  
    return(list(beta=beta,gamma=gamma,alpha=alpha,inv_hessianbt=invHB,
                scorebeta_i=scorebeta_i,eigfun=est_Eigfun,summary=summary(ss),
                num_distin=m,rn=VARselectindex,cr=cr,varbeta=varbeta,
                linkxbeta=linkxbeta,linkres=linkres,corindex=corindex))
    
  }
  else{
    return(paste0(count, "th simulation 's rn is", (VARselectindex),' and not converge'))
  }
}

#############################
n=200 # 500 or 1000
g=g1  # g2 or g3
gindex=1 # corresponding to g1
taulist=switch(gindex,taulist1,taulist2,taulist3)

crindex=1 # 1 is censoring rate 10%, 2 is censoring rate 30%
cvalue=censoringvalues[crindex] 

replicate_time=580
sfInit(parallel = TRUE,cpus = 80)
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
res=sfClusterApplyLB(1:replicate_time,parasimulation)
sfStop()
save(res, file=paste('newPointwise','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))
