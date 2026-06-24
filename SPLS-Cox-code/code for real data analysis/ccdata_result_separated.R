library(parallel)
library(snowfall)
rm(list=ls())

################ Original(Estimated) ##########################
# Load the 6000 results and find the valid results

load(paste('realdata_varbeta_separated24CCfpca1-27000','.Rdata',sep = ''))

# length(res)
valid=rep(NA,27000);i=1
for(count in 1:27000){
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

# pointwise CI for 1:100 points
get_ci1<-function(count){
  res_count=res[[count]]
  M=10000
  m=res_count$num_distin
  appr_bias_M=matrix(0,M,100)
  
  
  eigfun1=res_count$eigfun
  eigfun1[,(res_count$rn1+1):ncol(res_count$eigfun)]=0
  
  for(i in 1:M){
    
    norm_vec=matrix(rep(rnorm(m),res_count$rn),ncol = m,byrow = T) 
    score_pertur=res_count$scorebeta_i*norm_vec  
    score_pertur=apply(score_pertur,1,sum)
    
    asym_linear_beta=diag(res_count$varbeta)*score_pertur 
    
    appr_bias_M[i,]=eigfun1%*%(asym_linear_beta)
  }

  get_quantile<-function(x,alpha){return(quantile(x,alpha))}

  q_alpha_std_left=apply(appr_bias_M,2,get_quantile,alpha=0.025)
  q_alpha_std_right=apply(appr_bias_M,2,get_quantile,alpha=0.975)

  upbound_std=res_count$eigfun[,1:res_count$rn1]%*%res_count$beta[1:res_count$rn1]-q_alpha_std_left
  lowbound_std=res_count$eigfun[,1:res_count$rn1]%*%res_count$beta[1:res_count$rn1]-q_alpha_std_right

  original=res_count$eigfun[,1:res_count$rn1]%*%res_count$beta[1:res_count$rn1]
  return(list(upbound_std=upbound_std,lowbound_std=lowbound_std,
              original=original))
}

sfInit(parallel = TRUE,cpus = 14)
sfExportAll()
CI=sfClusterApplyLB(valid,get_ci1)
sfStop()
save(CI,file=paste('realdata_smooth_separatedCI1','.Rdata',sep = ''))


# pointwise CI for 101:200 points
get_ci2<-function(count){
  res_count=res[[count]]
  M=10000
  m=res_count$num_distin
  appr_bias_M=matrix(0,M,100)
  
  
  eigfun2=res_count$eigfun
  eigfun2[,1:res_count$rn1]=0
  
  for(i in 1:M){
    
    norm_vec=matrix(rep(rnorm(m),res_count$rn),ncol = m,byrow = T) 
    score_pertur=res_count$scorebeta_i*norm_vec  
    score_pertur=apply(score_pertur,1,sum)
    
    asym_linear_beta=diag(res_count$varbeta)*score_pertur 
    
    appr_bias_M[i,]=eigfun2%*%(asym_linear_beta)
  }
  
  get_quantile<-function(x,alpha){return(quantile(x,alpha))}
  
  q_alpha_std_left=apply(appr_bias_M,2,get_quantile,alpha=0.025)
  q_alpha_std_right=apply(appr_bias_M,2,get_quantile,alpha=0.975)
  
  upbound_std=res_count$eigfun[,(res_count$rn1+1):ncol(eigfun2)]%*%res_count$beta[(res_count$rn1+1):ncol(eigfun2)]-q_alpha_std_left
  lowbound_std=res_count$eigfun[,(res_count$rn1+1):ncol(eigfun2)]%*%res_count$beta[(res_count$rn1+1):ncol(eigfun2)]-q_alpha_std_right
  
  original=res_count$eigfun[,(res_count$rn1+1):ncol(eigfun2)]%*%res_count$beta[(res_count$rn1+1):ncol(eigfun2)]
  return(list(upbound_std=upbound_std,lowbound_std=lowbound_std,
              original=original))
}

sfInit(parallel = TRUE,cpus = 14)
sfExportAll()
CI=sfClusterApplyLB(valid,get_ci2)
sfStop()
save(CI,file=paste('realdata_smooth_separatedCI2','.Rdata',sep = ''))


### Use the top 25%  of the valid results, for final CI result to transfer into CC ######

whichCI=2

load(paste('realdata_smooth_separatedCI',whichCI,'.Rdata',sep = ''))
need=which(loglik_fit>=quantile(loglik_fit)[4])
up_result=low_result=ori_result=matrix(NA,nrow=100,ncol=length(need))
j=1
for (i in need){
  up_result[,j]=CI[[i]]$upbound_std
  low_result[,j]=CI[[i]]$lowbound_std
  ori_result[,j]=CI[[i]]$original
  j=j+1
}

library(ggplot2)
filename=paste('CCdata_separatedCI',whichCI,'.pdf',sep='')
# take the mean of the top 250 valid results
lower=apply(low_result,1,mean)
upper=apply(up_result,1,mean)
estimation=apply(ori_result,1,mean)
pdf(filename,width = 8,height = 6)
p1=ggplot() + 
  labs(x='s',y=expression(hat(beta)(s)),
       title="Estimated coefficient function",fill="")+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=seq(0,1,length.out=100), fill = "Pointwise Confidence Interval"), alpha = 0.4)+
  geom_line(aes(x=seq(0,1,length.out=100),y=estimation,color='Estimated Value')) + 
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
save(CCdata_PointwiseCI,file = paste0('CCdata_separated_PointwiseCI',whichCI,'.Rdata'))


### subplots ###
load("CCdata_separated_PointwiseCI1.Rdata") 
res1 <- CCdata_PointwiseCI

load("CCdata_separated_PointwiseCI2.Rdata")
res2 <- CCdata_PointwiseCI

#
plot_data <- rbind(
  data.frame(
    s = seq(0, 1, length.out = length(res1$estimation)),
    lower = res1$lower,
    upper = res1$upper,
    estimation = res1$estimation,
    coef = "hat(beta)[1](s)"
  ),
  data.frame(
    s = seq(0, 1, length.out = length(res2$estimation)),
    lower = res2$lower,
    upper = res2$upper,
    estimation = res2$estimation,
    coef = "hat(beta)[2](s)"
  )
)

pdf("CCdata_separatedCI_both.pdf", width = 10, height = 5)
p <- ggplot(plot_data, aes(x = s)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#F8766D", alpha = 0.4) +
  geom_line(aes(y = estimation), color = "black") +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  facet_wrap(~ coef, nrow = 1, labeller = label_parsed) +
  labs(
    x = "s",
    y = NULL,
    title = "Estimated coefficient functions"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "none"
  )
print(p)
dev.off()


#### Estimated psi function, using the best result ####

linkxbeta=res[[valid[best_res]]]$linkxbeta
linkres=res[[valid[best_res]]]$linkres

library(ggplot2)
filename='CCfpca_linkfun.pdf'
pdf(filename,width = 8,height = 6)
p1 <- ggplot() +
  labs(x=expression(hat(tilde(xi))^T*hat(tilde(b))),y=expression(hat(psi)("·")),title="Estimated link function",fill="")+
  geom_line(aes(x = as.vector(linkxbeta),y = as.vector(linkres),group = 1),alpha = 1, linewidth = 0.8)+
  #geom_line(aes(x=linkxbeta,y=linkres),alpha=1,linewidth=0.8)+
  theme_bw() + ylim(-2.5,2.5)+ xlim(-0.04,0.04)+
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

