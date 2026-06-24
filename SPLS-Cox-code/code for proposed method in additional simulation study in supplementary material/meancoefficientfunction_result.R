
rm(list=ls())
library(snowfall)
library(parallel)
library(locpol)
library(MASS)
library(ggplot2)


set.seed(2024)
simutimes=500
gindex=1
cvalue=0.1

# get.name='newPointwise'
# get.name='newPointwise.weibull'
get.name='newPointwise.piecewise'

h <- calh <- 0.01
p <- 1/h+1 #p=101
t <- seq(0,1,calh)#

#eigenfunction, FPC
eigenf1 <- function(x) sin(2*pi*x)*sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)*sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)*sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)*sqrt(2)

#beta0(s)
betafun <- function(x){
  index=c(1,1/4,1/9,1/16)
  true_index=index/sqrt(sum(index^2))
  beta=true_index[1]*eigenf1(x)+true_index[2]*eigenf2(x)+
    true_index[3]*eigenf3(x)+true_index[4]*eigenf4(x)
}
truebeta<-betafun(t)




n_list <- c(200, 500, 1000)
g_list <- c(1, 2, 3)

plot_data_list <- list()
idx <- 1

for (n in n_list) {
  for (gindex in g_list) {
    
    file_name <- paste(get.name,'g',gindex,'n',n,'cr',cvalue,'.Rdata',sep = '')
    

    e <- new.env()
    load(file_name, envir = e)
    res <- e$res
    
    # 找 valid
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
                            levels = c(200, 500, 1000),
                            labels = c(expression(n==200),
                                       expression(n==500),
                                       expression(n==1000)))

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