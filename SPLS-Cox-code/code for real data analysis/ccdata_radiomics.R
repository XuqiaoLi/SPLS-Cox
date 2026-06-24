rm(list=ls())
library(survival)
library(glmnet)

######### 202606 Cox model with variable selection ############

### clinical + CC volumetric data ###
load("clinical.dat") 
scoredata=scoredata[-331,]
CCinfo=read.csv('CCinfo.csv',header = T)

ZS=scoredata[,4:15]
info=apply(CCinfo[,2:6]/CCinfo[,7],2,scale)

ZS=cbind(ZS,info)
time=scoredata[,1]
event=scoredata[,2]


## radiomics features
cleandata=read.csv('radiomics_features.csv',header = T)
cleandata=cleandata[,-1]
radiomics_features=as.matrix(cleandata)

designmatrix <- cbind(radiomics_features, ZS)


x <- as.matrix(designmatrix)

y <- Surv(time = as.numeric(time), event = as.numeric(event))

set.seed(123)

cvfit <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1
)

plot(cvfit)

## lambda.min
coef_min <- coef(cvfit, s = "lambda.min")

selected_result <- data.frame(
  variable = rownames(coef_min),
  coefficient = as.numeric(coef_min)
)

selected_result <- subset(selected_result, coefficient != 0)
selected_result$HR <- exp(selected_result$coefficient)

selected_result


# ### coxph ###
# coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
# selected <- which(coef_min[, 1] != 0)
# 
# selected_data <- data.frame(
#   time = time,
#   event = event,
#   designmatrix[, selected]
# )
# 
# fit_refit <- coxph(Surv(time, event) ~ ., data = selected_data)
# 
# summary(fit_refit)



############  Cox model with glmnet, for cross-validated Cindex ############
library(glmnet)
library(survival)
library(survcomp)



cleandata=read.csv('radiomics_features.csv',header = T)
cleandata=cleandata[,-1]
radiomics_features=as.matrix(cleandata)

designmatrix <- cbind(radiomics_features, ZS)


x <- as.matrix(designmatrix)
y <- Surv(time, event)

set.seed(123)

B <- 1000

n.train=250
n.test=n-n.train

corindex.test <- matrix(
  NA_real_,
  nrow = B,
  ncol = 2,
  dimnames = list(NULL, c("C-index", "SE"))
)

for (b in 1:B) {

  ## 1.
  train.id <- sample(1:nrow(x), size = n.train, replace = FALSE)
  
  ## 
  test.id <- setdiff(1:nrow(x), train.id)
  
  x.train <- x[train.id, , drop = FALSE]
  x.test  <- x[test.id, , drop = FALSE]
  
  y.train <- Surv(
    time  = time[train.id],
    event = event[train.id]
  )
  
  ## 2. 
  cvfit <- cv.glmnet(
    x = x.train,
    y = y.train,
    family = "cox",
    alpha = 1,             # LASSO
    nfolds = 5
  )
  
  ## 3. 
  pre.new <- as.numeric(
    predict(
      cvfit,
      newx = x.test,
      s = "lambda.min",
      type = "link"
    )
  )
  
  ## 4. 
  cindex.pre <- concordance.index(
    x = pre.new,
    surv.time = time[test.id],
    surv.event = event[test.id],
    method = "noether"
  )
  
  corindex.test[b, ] <- c(
    cindex.pre$c.index,
    cindex.pre$se
  )
  
  if (b %% 100 == 0) {
    cat("Finished:", b, "\n")
  }
}

mean(corindex.test[,1])
sd(corindex.test[,1])