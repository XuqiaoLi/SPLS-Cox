rm(list=ls())
library(snowfall)
library(parallel)
library(MASS)
library(locpol)
simutimes=500

###### Copy from simulation_shape ###########
thisseed=2024
cleandata=read.csv('ccdata.csv',header = F)
observeX=as.matrix(cleandata)

n=372
p=200
q=6
h <- calh <-1  
t <- seq(1,ncol(cleandata),h) 
mygamma<-rep(0.2,q) 

gindex=1 # 2 or 3
crindex=1 
censoringvalues<-c(0.1) #censoring rate
cvalue=censoringvalues[crindex] 
load(paste('newShape','g',gindex,'cr',cvalue,'.Rdata',sep =''))


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


valid=rep(NA,length(res))
j=1
for(i in 1:length(res)){
  if(!is.character(res[[i]])){valid[j]=i;j=j+1}
}

valid=na.omit(valid)
if(length(valid)>simutimes){valid=valid[1:simutimes]}

######Multiplier process for pointwise confidence interval#########
get_PWCB<-function(count){
  set.seed(2024*gindex*count)
  res_count=res[[count]]
  M=10000
  m=res_count$num_distin
  appr_bias_M=appr_bias=matrix(0,M,length(t))
  
  for(i in 1:M){
    
    norm_vec=matrix(rep(rnorm(m),res_count$rn),ncol = m,byrow = T) 
    score_pertur=res_count$scorebeta_i*norm_vec 
    score_pertur=apply(score_pertur,1,sum)
    
    asym_linear_beta=res_count$varbeta%*%score_pertur
    appr_bias_M[i,]=res_count$eigfun%*%(asym_linear_beta)
  }
  
  get_quantile<-function(x,alpha){return(quantile(x,alpha))}
  
  q_alpha_std_left=apply(appr_bias_M,2,get_quantile,alpha=0.025)
  q_alpha_std_right=apply(appr_bias_M,2,get_quantile,alpha=0.975)
  
  upbound_std=res_count$eigfun%*%res_count$beta-q_alpha_std_left
  lowbound_std=res_count$eigfun%*%res_count$beta-q_alpha_std_right
  
  return(list(upbound_std=upbound_std,lowbound_std=lowbound_std))
}

sfInit(parallel = TRUE,cpus = 14)
sfLibrary(MASS)
sfExportAll()
res_pointwise=sfClusterApplyLB(valid,get_PWCB)
sfStop()
save(res_pointwise, file=paste('shape_pointwise_result','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))
#############

library(ggplot2)
get_value<-function(n,gindex,cvalue){
  load(paste('shape_pointwise_result','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))
  if_cover=matrix(NA,nrow=simutimes,ncol=p)
  upbound=lowbound=upbound_std=lowbound_std=matrix(0,nrow=simutimes,ncol=p)
  
  for(i in 1:simutimes){
    result=res_pointwise[[i]]
    
    upbound_std[i,]=result$upbound_std
    lowbound_std[i,]=result$lowbound_std
    
    if_cover[i,]=(upbound_std[i,]>truebeta)&(lowbound_std[i,]<truebeta)
  }
  
  rate=apply(if_cover,MARGIN = 2,FUN = mean)
  
  df_value=data.frame(rate)
  df_value$t=t
  df_value$linkfun=paste("psi","[",gindex,"]",sep = "")
  df_value$cr=paste("cr=",cvalue,sep="") 
  df_value$n=n
  
  return(df_value)
}
n372g1cr0.1=get_value(n,1,0.1)
n372g2cr0.1=get_value(n,2,0.1)
n372g3cr0.1=get_value(n,3,0.1)

# alldata=rbind(n372g1cr0.1,n372g2cr0.1,n372g3cr0.1)
alldata=rbind(n372g2cr0.1,n372g3cr0.1) 

filename='pointwise_CP_shape.pdf'
pdf(filename,width = 15,height = 5)
ggplot(alldata, aes(x=t,y=rate)) +
  geom_point(shape=1,size=1.5)+facet_grid(~linkfun,scales='free',labeller = labeller(linkfun = label_parsed))+
  geom_hline(aes(yintercept=0.95),colour='red',linetype="dashed",size=0.8)+
  xlab('s')+ylab('covarage probabilities') + theme_bw()
dev.off()


# overall coverage probabilities
round(mean(n372g1cr0.1$rate),3)
round(mean(n372g2cr0.1$rate),3)
round(mean(n372g3cr0.1$rate),3)