rm(list=ls())
library(snowfall)
library(parallel)
library(locpol)
library(MASS)
#############Pointwise coverage probability result#################
set.seed(2024)
simutimes=500
n=500
gindex=1 # 2 or 3
cvalue=0.1
load(paste('newPointwise','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))

h <- calh <- 0.01
p <- 1/h+1 
t <- seq(0,1,calh)

#eigenfunction, FPC
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
  appr_bias_M=matrix(0,M,length(t))
  
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
save(res_pointwise, file=paste('pointwise_result','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))


#####Plot the pointwise cp result#############
library(ggplot2)
filename='pointwise_CP.pdf'
get_value<-function(n,gindex,cvalue){
  load(paste('pointwise_result','g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = ''))
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

n500g1cr0.1=get_value(500,1,0.1)
n500g2cr0.1=get_value(500,2,0.1)
n500g3cr0.1=get_value(500,3,0.1)

# alldata=rbind(n500g1cr0.1,n500g2cr0.1,n500g3cr0.1)
alldata=rbind(n500g1cr0.1,n500g2cr0.1)

pdf(filename,width = 15,height = 5)
ggplot(alldata, aes(x=t,y=rate)) +  
  geom_point(shape=1,size=1.5)+facet_grid(~linkfun,scales='free',labeller = labeller(linkfun = label_parsed))+
  geom_hline(aes(yintercept=0.95),colour='red',linetype="dashed",size=0.8)+
  xlab('s')+ylab('covarage probabilities') + theme_bw() 
dev.off()


# overall coverage probabilities
round(mean(n500g1cr0.1$rate),3)
round(mean(n500g2cr0.1$rate),3)
round(mean(n500g3cr0.1$rate),3)
