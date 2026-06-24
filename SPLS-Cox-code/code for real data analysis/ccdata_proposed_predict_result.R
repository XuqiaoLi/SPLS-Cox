library(parallel)
library(snowfall)
rm(list=ls())




##############  C index cross validated 202606 #######################

load('realdata_16CCfpca_Prediction26_test122.Rdata')
replication_num=length(res)
cindex.train.all=cindex.test.all=matrix(NA,nrow = replication_num,ncol = 2)
for(i in 1:replication_num){
  cindex.train.all[i,]=res[[i]]$corindex.train
  cindex.test.all[i,]=res[[i]]$corindex.test
}

cindex.train.all=cindex.train.all[1:1000,]
cindex.test.all=cindex.test.all[1:1000,]

apply(cindex.train.all,2,mean)
apply(cindex.test.all,2,mean)
sd(cindex.test.all[,1])


###### C index cross validated (separated)202606 ###############

load('realdata_separated_24CCfpca_Prediction2026_test122.Rdata')

replication_num=length(res)
valid=rep(NA,length(res));i=1
for(count in 1:length(res)){
  if(class(res[[count]])=='list'){valid[i]=count;i=i+1}else{next}
}
valid=na.omit(valid)

cindex.train.all=cindex.test.all=matrix(NA,nrow = length(valid),ncol = 2)

i=1
for (count in valid){
  res_count=res[[count]]
  cindex.train.all[i,]=res_count$corindex.train
  cindex.test.all[i,]=res_count$corindex.test
  i=i+1
}

### 1000 times cross validated results ###
cindex.train.all=cindex.train.all[1:1000,]
cindex.test.all=cindex.test.all[1:1000,]

apply(cindex.train.all,2,mean)
apply(cindex.test.all,2,mean)
sd(cindex.test.all[,1])