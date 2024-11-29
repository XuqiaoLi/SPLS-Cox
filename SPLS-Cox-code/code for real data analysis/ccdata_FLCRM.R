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
scoredata=scoredata[-331,]
cleandata=read.csv('ccdata.csv',header = F)
CCinfo=read.csv('CCinfo.csv',header = T)

observeX=as.matrix(cleandata)
ZS=scoredata[,4:15]
info=apply(CCinfo[,2:6]/CCinfo[,7],2,scale)
colnames(info)=c()

ZS=cbind(ZS,info)
time=scoredata[,1]
event=scoredata[,2]

n=372;p=200;q=17;
t <- seq(0,2,length.out=200)
h=t[2]

###################################
X<-matrix(NA,n,p)
for (i in 1:n){
  dataframe <- data.frame(t)
  dataframe$observeX<-observeX[i,]
  lpfit <- locpol(observeX~t,dataframe, xeval=t)
  X[i,]<-lpfit$lpFit[,2]
}

scaleX<-scale(X,center=TRUE,scale=FALSE)
Eigvec=svd(scaleX)$v 
est_Eigval=svd(scaleX)$d^2*h/(n-1)

Rnselectindex<-16 

est_Eigfun=Eigvec[,1:Rnselectindex]/sqrt(h)

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

fullcoxcurve<-survfit(fullcoxresult)

finalresult<-summary(fullcoxresult)

summarycoefficient<-finalresult$coefficients

coefficient<-summarycoefficient[,1]
secoefficient<-summarycoefficient[,3]

#### Estimated Beta
AICestimatebeta<-est_Eigfun%*%as.matrix(coefficient[1:Rnselectindex])

#### Estimated Gamma
AICestimategamma<-coefficient[-(1:Rnselectindex)]

#### C-index
cindex<-concordance.index(predict(fullcoxresult),surv.time = time, 
                          surv.event = event ,method = "noether")
corindex<-c(cindex$c.index,cindex$se)
round(corindex,3)


coef_vec=round(summarycoefficient[-(1:Rnselectindex),1],3) # coef
se_vec=round(summarycoefficient[-(1:Rnselectindex),3],3) # se
pvalue_vec=round(summarycoefficient[-(1:Rnselectindex),5],3) # p-value
aa=bb=cc=''
for(i in 1:length(coef_vec)){
  aa=paste(aa,coef_vec[i],sep = ' & ')
  bb=paste(bb,se_vec[i],sep = ' & ')
  cc=paste(cc,pvalue_vec[i],sep = ' & ')
}


# ###### prediction using cross-validation analysis #####
# thisseed=2020
# n.test=122
# n.train=n-n.test
# FLCRM_estimate_predict<-function(init_seed){
#   set.seed(thisseed+init_seed)
#   index.train=sample(x = n,size = n.train,replace = F) 
#   newdataframe.train=newdataframe[index.train,]
#   newdataframe.test=newdataframe[-index.train,]
#   
#   fullcoxresult<-coxph(formula=totalformular, data=newdataframe.train)
#   
#   pre.new=predict(fullcoxresult,newdata = newdataframe.test[,-c(1,2)],type = 'lp')
#   cindex.pre<-concordance.index(pre.new,surv.time = newdataframe.test$time, 
#                                 surv.event = newdataframe.test$event ,method = "noether")
#   corindex.test=c(cindex.pre$c.index,cindex.pre$se)
#   
#   return(corindex.test=corindex.test)}
# 
# replication.num=500
# cindex.all=matrix(NA,nrow = replication.num,ncol = 2)
# for(i in 1:replication.num){
#   cindex.all[i,]=FLCRM_estimate_predict(i)
# }
# apply(cindex.all,2,mean)