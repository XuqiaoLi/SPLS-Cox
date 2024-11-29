library(parallel)
library(snowfall)
rm(list=ls())

################ Estimated ##########################
# Load the 6000 results and find the valid results, around 1600.
load(paste('realdata_',16,'CCfpca','.Rdata',sep = ''))  
valid=rep(NA,6000);i=1
for(count in 1:6000){
  if(class(res[[count]])=='list'){valid[i]=count;i=i+1}else{next}
}
valid=na.omit(valid)

cindex_fit=matrix(NA,length(valid),2)
loglik_fit=rep(NA,length(valid))
i=1
for (count in valid){
  res_count=res[[count]]
  cindex_fit[i,]=res_count$corindex
  loglik_fit[i]=res_count$loglik
  i=i+1
}
loglik_fit=na.omit(loglik_fit)
best_res=which.max(loglik_fit)

############## Generate the multiplier process ############
get_ci<-function(count){
  res_count=res[[count]]
  M=10000
  m=res_count$num_distin
  appr_bias_M=appr_bias=matrix(0,M,200)
  
  for(i in 1:M){
    
    norm_vec=matrix(rep(rnorm(m),res_count$rn),ncol = m,byrow = T) 
    score_pertur=res_count$scorebeta_i*norm_vec 
    score_pertur=apply(score_pertur,1,sum)
    
    asym_linear_beta=diag(res_count$varbeta)*score_pertur #multiply the diagonal matrix of varbeta

    appr_bias_M[i,]=res_count$eigfun%*%(asym_linear_beta)
  }

  get_quantile<-function(x,alpha){return(quantile(x,alpha))}

  q_alpha_std_left=apply(appr_bias_M,2,get_quantile,alpha=0.025)
  q_alpha_std_right=apply(appr_bias_M,2,get_quantile,alpha=0.975)

  upbound_std=res_count$eigfun%*%res_count$beta-q_alpha_std_left
  lowbound_std=res_count$eigfun%*%res_count$beta-q_alpha_std_right

  original=res_count$eigfun%*%res_count$beta
  return(list(upbound_std=upbound_std,lowbound_std=lowbound_std,
              original=original))
}

sfInit(parallel = TRUE,cpus = 14)
sfExportAll()
CI=sfClusterApplyLB(valid,get_ci)
sfStop()
save(CI,file=paste('realdata_smooth_CI','.Rdata',sep = ''))

### Use the top 250  of the valid results, for final CI result to transfer into CC, named by  "CCdata_PointwiseCI.Rdata" ######
load(paste('realdata_smooth_CI','.Rdata',sep = ''))
need=which(loglik_fit>=quantile(loglik_fit)[4]) # top 25% quantile of loglik
up_result=low_result=ori_result=matrix(NA,nrow=200,ncol=length(need))
j=1
for (i in need){
  up_result[,j]=CI[[i]]$upbound_std
  low_result[,j]=CI[[i]]$lowbound_std
  ori_result[,j]=CI[[i]]$original
  j=j+1
}

library(ggplot2)
filename=paste('CCdata_CI','.pdf',sep='')
# take the mean of the top 250 valid results
lower=apply(low_result,1,mean)
upper=apply(up_result,1,mean)
estimation=apply(ori_result,1,mean)
pdf(filename,width = 8,height = 6)
p1=ggplot() + 
  labs(x='s',y=expression(hat(beta)(s)),
       title="Estimated coefficient function",fill="")+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=seq(0,2,length.out=200), fill = "Pointwise Confidence Interval"), alpha = 0.4)+
  geom_line(aes(x=seq(0,2,length.out=200),y=estimation,color='Estimated Value')) + 
  #scale_fill_manual(values = c('Pointwise Confidence Interval'='grey'))+
  scale_color_manual(values = c("Estimated Value"="black"),labs(color = ''))+
  theme_bw() + geom_hline(aes(yintercept=0), colour="red", linetype="dashed") +
  theme(plot.title = element_text(hjust = 0.5,vjust=1,face = 'bold' ),
        axis.title.x=element_text(vjust=0.5, size=15),
        axis.title.y=element_text(vjust=0.5, size=15),
        axis.text.x=element_text(vjust=1,size=11),
        axis.text.y=element_text(vjust=1,size=11),
        # legend.position='bottom',
        legend.position='none',
        legend.text = element_text(size=15))
print(p1)
dev.off()

CCdata_PointwiseCI=list(upper=upper,lower=lower,estimation=estimation)
save(CCdata_PointwiseCI,file = 'CCdata_PointwiseCI.Rdata')

#### Estimated psi function, using the best result ####

linkxbeta=res[[valid[best_res]]]$linkxbeta
linkres=res[[valid[best_res]]]$linkres

library(ggplot2)
filename='CCfpca_linkfun.pdf'
pdf(filename,width = 8,height = 6)
p1 <- ggplot() +
  labs(x=expression(hat(xi)^T*hat(b)),
       y=expression(hat(psi)("·")),title="Estimated link function",fill="")+
  geom_line(aes(x=linkxbeta,y=linkres),alpha=1,size=0.8)+
  ylim(-2.5,2.5)+xlim(-0.065,0.065)+theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5,vjust=1,face = 'bold' ),
        axis.title.x=element_text(vjust=0.5, size=15),
        axis.title.y=element_text(vjust=0.5, size=15),
        axis.text.x=element_text(vjust=1,size=11),
        axis.text.y=element_text(vjust=1,size=11))
#jpeg(file = "CCfpca_linkfun.png",width=1400,height=800,res=100)
print(p1)
dev.off()



############### gamma and C-index ############
coefficient=round(res[[valid[best_res]]]$summary$coefficients,3)
coef_vec=pvalue_vec=se_vec=''
for(i in 6:22){
  coef_vec=paste(coef_vec,coefficient[i,1],sep = ' & ')
  pvalue_vec=paste(pvalue_vec,coefficient[i,5],sep = ' & ')
  se_vec=paste(se_vec,coefficient[i,3],sep = ' & ')
  }


round(res[[valid[best_res]]]$corindex,3) # or cindex_fit[best_res,]

