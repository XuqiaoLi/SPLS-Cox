
rm(list=ls())
library(snowfall)
library(parallel)
library(locpol)
library(MASS)
library(ggplot2)


load("clinical.dat") 
scoredata=scoredata[-331,]
cleandata=read.csv('ccdata.csv',header = F)
observeX=as.matrix(cleandata)
############### RMSE table   ###################
simutimes=500
cvalue=0.1

n.shape=372
n=1000 # 500 #372
p=200;q=6 
h <- calh <-1  # 1/100
t <- seq(1,ncol(cleandata),h) #seq(0.01,2,calh)->t
mygamma<-rep(0.2,q)
LengthAlpha=q


### CC data for FPCA
X<-matrix(NA,n.shape,p)
for (i in 1:n.shape){
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
shape_est_Eigval=shape_Singval^2*h/(n.shape-1) 
shape_est_Eigfun=shape_Eigvec/sqrt(h)
shape_est_Eigfun_number=length(shape_est_Eigval)

# Only take the first ** eigen functions that contribute to total 70%/85% variance, to generate data
shape_est_Eigfun_number=9
shape_est_Eigfun=shape_est_Eigfun[,1:shape_est_Eigfun_number]

#beta coefficient
# beta_coef=rep(1,shape_est_Eigfun_number)/seq(1,shape_est_Eigfun_number)
beta_coef=rep(1,shape_est_Eigfun_number)/(seq(1,shape_est_Eigfun_number))^2
# beta_coef[4:shape_est_Eigfun_number]=0
beta_coef=beta_coef/sqrt(sum(beta_coef**2)) 


truebeta <- c(beta_coef %*% t(shape_est_Eigfun))



get.name='newShape_smooth'




# n_list <- c(372, 500, 1000)
n_list <- c(372, 744)
g_list <- c(1, 2, 3)

plot_data_list <- list()
idx <- 1

for (n in n_list) {
  for (gindex in g_list) {
    
    file_name <- paste(get.name,'g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = '')
    

    e <- new.env()
    load(file_name, envir = e)
    res <- e$res
    

    valid <- rep(NA, length(res))
    j <- 1
    for (i in 1:length(res)) {
      if (!is.character(res[[i]])) {
        valid[j] <- i
        j <- j + 1
      }
    }
    valid <- na.omit(valid)
    
    if (length(valid) > simutimes) {
      valid <- valid[1:simutimes]
    }
    
    if (length(valid) < simutimes) {
      warning(paste("For", file_name, "only", length(valid), "valid results are available."))
    }
    

    beta_hat <- matrix(NA, nrow = length(valid), ncol = p)
    
    for (count in 1:length(valid)) {
      beta_hat[count, ] <- t(res[[valid[count]]]$eigfun %*% res[[valid[count]]]$beta)
    }
    

    mean_beta <- colMeans(beta_hat)
    

    df_true <- data.frame(
      s = t,
      beta = truebeta,
      type = "True beta(s)",
      n_label = factor(n, levels = n_list),
      g_label = factor(gindex, levels = g_list)
    )
    
    df_est <- data.frame(
      s = t,
      beta = mean_beta,
      type = "Mean estimated beta(s)",
      n_label = factor(n, levels = n_list),
      g_label = factor(gindex, levels = g_list)
    )
    
    plot_data_list[[idx]] <- rbind(df_true, df_est)
    idx <- idx + 1
  }
}

plot_data <- do.call(rbind, plot_data_list)


plot_data$n_label <- factor(plot_data$n_label,
                            levels = c(372, 744),
                            labels = c(expression(n==372),
                                       expression(n==744)))

plot_data$g_label <- factor(plot_data$g_label,
                            levels = c(1, 2, 3),
                            labels = c(expression(psi[1]),
                                       expression(psi[2]),
                                       expression(psi[3])))


plot_data$n_lab_char <- factor(as.character(plot_data$n_label),
                               levels = as.character(unique(plot_data$n_label)))

plot_data$g_lab_char <- factor(as.character(plot_data$g_label),
                               levels = as.character(unique(plot_data$g_label)))


p1 <- ggplot(plot_data, aes(x = s, y = beta, linetype = type, color = type)) +
  geom_line(linewidth = 0.8) +
  facet_grid(n_lab_char ~ g_lab_char,
             scales = "fixed",
             labeller = labeller(
               n_lab_char = label_parsed,
               g_lab_char = label_parsed
             )) +
  labs(
    x = "s",
    y = expression(hat(beta)(s)),
    linetype = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "none"
  )
print(p1)

ggsave(paste0("True_vs_mean_estimated_beta_facet_",get.name,".pdf"), p1, width = 12, height = 9)