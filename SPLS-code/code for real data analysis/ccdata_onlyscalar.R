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

observeX=(cleandata)
ZS=scoredata[,4:15]
info=apply(CCinfo[,2:6]/CCinfo[,7],2,scale)

ZS=cbind(ZS,info)
time=scoredata[,1]
event=scoredata[,2]

thisseed=1000
set.seed(thisseed)


#### only scalar covariate(clinical and CC size) ####
onlyscalar <- try(coxph(Surv(time,event) ~ ZS,  
                method="breslow", eps=1e-7, iter.max=20), TRUE)
finalresult=summary(onlyscalar)
finalresult$coefficients
round(finalresult$coefficients[,c(1,3,5)],3)
cindex=concordance.index(predict(onlyscalar),surv.time = time, 
                  surv.event = event, method = "noether")
round(c(cindex$c.index,cindex$se),3)



### only scalar covariate(clinical) ####
ZS=scoredata[,4:15]
onlyscalar <- try(coxph(Surv(time,event) ~ ZS,  
                        method="breslow", eps=1e-7, iter.max=20), TRUE)
finalresult=summary(onlyscalar)
finalresult$coefficients
round(finalresult$coefficients[,c(1,3,5)],3)
cindex=concordance.index(predict(onlyscalar),surv.time = time, 
                         surv.event = event, method = "noether")
round(c(cindex$c.index,cindex$se),3)