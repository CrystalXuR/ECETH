
##---------------------- HTE calibration paper data application using SPRINT data -------------------##
##----------------------           Crystal Xu               10/8/2021             -------------------##
library(survival); library(glmnet); library(grf)
# R-learner lasso approach 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Analyses Stanford Team/HTE_experiments/rlearner_surv")
files.sources = list.files()
sapply(files.sources, source)
base_surv <- function(fit, Y, D, x, lambda){
  data <- data.frame(t_event=Y, event=D, x)
  tab <- data.frame(table(data[data$event == 1, "t_event"])) 
  y <- as.numeric(as.character(sort(unique(tab[,1]))))
  d <- tab[,2]  # number of events at each unique time                               
  
  betaHat <- as.vector((fit$glmnet.fit$beta)[,fit$lambda==lambda])
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(x[data$t_event >= y[l], rownames(fit$glmnet.fit$beta)] %*% betaHat))
  }    
  
  S0 <- exp(-cumsum(h0))
  outcome <- data.frame(time=y,survival=S0)
  outcome
}
pred_surv <- function(fit, S0, x, times, lambda){
  link <- predict(fit$glmnet.fit,x,type = "link")[,fit$lambda==lambda]
  colnames(link) <- NULL
  
  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }
  
  surv <- S0_t^exp(link)
  surv
}
pred_surv_preval <- function(fit, S0, times, lambda){
  link <- fit$fit.preval[,!is.na(colSums(fit$fit.preval))][, fit$lambda[!is.na(colSums(fit$fit.preval))] == lambda] 
  colnames(link) <- NULL
  
  if(length(times)>1){
    S0_t <- rep(NA, length(times))
    for (i in 1:length(times)){
      S0_t[i] <- S0$survival[S0$time>=times[i]][1]
    }
  }else{
    S0_t <- S0$survival[S0$time>=times][1]
  }
  
  surv <- S0_t^exp(link)
  surv
}

# Load data 
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE/Data Application")
blvars <- read.csv("blvars.csv")
outcome <- read.csv("outcome.csv")
tmpdat <- merge(blvars,outcome, by="MASKID",all=T)
tmpdata <- tmpdat[complete.cases(tmpdat),]; dim(tmpdata) # 8969 43
row.names(tmpdata) <- 1:dim(tmpdata)[1]
t50 <- floor(365.25*3.26)

# Split data to training and testing 
set.seed(13)
index <- sample(1:dim(tmpdata)[1], floor(dim(tmpdata)[1]*0.5))
traindat <- tmpdata[-index, ]
testdat <- tmpdata[index, ]

# Estimate ATE using Cox PH 
cox_ate <- coxph(Surv(t_cvds, cvd) ~ INTENSIVE, data = traindat)
base_haz <- basehaz(cox_ate, centered = FALSE)[basehaz(cox_ate)$time==t50,1]
testdat1 <- data.frame(rep(1,dim(testdat)[1]));names(testdat1) <- "INTENSIVE" 
testdat0 <- data.frame(rep(0,dim(testdat)[1]));names(testdat0) <- "INTENSIVE" 

# survival probabilities
est_S1_cvd <- exp(-base_haz)^exp(testdat1$INTENSIVE * cox_ate$coefficients)
est_S1_UB <- exp(-base_haz)^exp(testdat1$INTENSIVE * (cox_ate$coefficients - 1.96*0.1088))
est_S1_LB <- exp(-base_haz)^exp(testdat1$INTENSIVE * (cox_ate$coefficients + 1.96*0.1088))
est_S0_cvd <- exp(-base_haz)

ate <- est_S1_cvd-est_S0_cvd  
ate_LB <- est_S1_LB-est_S0_cvd
ate_UB <- est_S1_UB-est_S0_cvd
print(round(c(ate[1], ate_LB[1], ate_UB[1]),4)) # 0.0164 (0.0063 0.0246)

# S-learner Cox PH w/ all interactions
scoxph_fit <- scoxph(x = as.matrix(traindat[, 3:39]),
                     w = traindat$INTENSIVE,
                     y = traindat$t_cvds, 
                     D = traindat$cvd, 
                     times = t50)
Scoxph_cate <- predict(scoxph_fit, newx = as.matrix(testdat[, 3:39]), times = t50)
summary(Scoxph_cate)

# S-learner lasso w/ all interactions
slasso_fit <- slasso_surv(x = as.matrix(traindat[, 3:39]),
                          w = traindat$INTENSIVE,
                          y = traindat$t_cvds, 
                          D = traindat$cvd, 
                          lambda_choice = "lambda.min",
                          times = t50)
Slasso_cate <- predict(slasso_fit, newx = as.matrix(testdat[, 3:39]), times = t50)
summary(Slasso_cate) # 0.01791 0.02486 0.03044 0.03695 0.04077 0.19569 

# R-learner lasso
# fit censoring weight 
set.seed(13)
x <- as.matrix(traindat[, 3:39])
penalty_factor_nuisance_c <- rep(1, (2*dim(x)[2] + 1))
x_int <- data.frame(cbind(as.numeric(traindat$INTENSIVE - 0.5) * cbind(1, x), x))
names(x_int)[1] <- "INTENSIVE"
x_int <- as.matrix(x_int)
foldid <- sample(rep(seq(10), length = length(traindat$INTENSIVE)))
c_fit <- cv.glmnet(x_int, 
                   Surv(traindat$t_cvds, 1-traindat$cvd),
                   family = "cox",
                   foldid = foldid,
                   keep = TRUE,
                   lambda = NULL,
                   alpha = 1,
                   penalty.factor = penalty_factor_nuisance_c)
c_lambda_min <- c_fit$lambda[which.min(c_fit$cvm[!is.na(colSums(c_fit$fit.preval))])]
S0 <- base_surv(c_fit, traindat$t_cvds, 1-traindat$cvd, x_int, lambda = c_lambda_min)
cent <- rep(t50, length(traindat$t_cvds))
cent[traindat$cvd==1] <- traindat$t_cvds[traindat$cvd==1]
c_hat <- pred_surv_preval(c_fit, S0, times = cent, lambda = c_lambda_min)
summary(c_hat)
#c_hat[c_hat==0] <- min(c_hat[c_hat!=0])
#c_hat[c_hat<0.1] <- 0.1

rlasso_fit <- rlasso(x = as.matrix(traindat[, 3:39]), 
                     w = traindat$INTENSIVE, 
                     y = traindat$t_cvds, 
                     D = traindat$cvd,
                     alpha = 1,
                     p_hat = rep(0.5, dim(traindat)[1]),
                     k_folds = 10,
                     c_hat = c_hat,
                     lambda_choice = "lambda.min",
                     times = t50)
rlasso_cate <- predict(object = rlasso_fit, newx = as.matrix(testdat[,3:39])) 
summary(rlasso_cate)


#-------------------------- Evaluate model calibration using proposed metric ----------------#
source("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/risk-vs-hte-R/survival_aipw_score.R")
rownames(testdat) <- 1:dim(testdat)[1]
nboot <- 1000
theta_hat_scoxph <- theta_hat_rlasso <- theta_hat_slasso <- est_ate <- Scoxph_ate <- rlasso_ate <- Slasso_ate <- rep(NA, nboot)
for (B in 1:nboot){
  index <- sample(1:dim(testdat)[1], dim(testdat)[1], replace = TRUE)
  tmpdat <- testdat[unique(index),]
  unique_times <- sort(intersect(unique(tmpdat[tmpdat$INTENSIVE==1, "t_cvds"]), unique(tmpdat[tmpdat$INTENSIVE==0, "t_cvds"])))
  unique_times <- sort(unique(tmpdat$t_cvds))
  
  # survival curve for treated 
  S1 <- survival_forest(X = tmpdat[tmpdat$INTENSIVE==1, 3:39],
                        Y = tmpdat[tmpdat$INTENSIVE==1, "t_cvds"], 
                        D = tmpdat[tmpdat$INTENSIVE==1, "cvd"],
                        failure.times = unique_times)
  surv1 <- predict(S1, tmpdat[, 3:39], failure.times = unique_times)$predictions
  
  S0 <- survival_forest(X = tmpdat[tmpdat$INTENSIVE==0, 3:39],
                        Y = tmpdat[tmpdat$INTENSIVE==0, "t_cvds"], 
                        D = tmpdat[tmpdat$INTENSIVE==0, "cvd"],
                        failure.times = unique_times)
  surv0 <- predict(S0, tmpdat[, 3:39], failure.times = unique_times)$predictions
  
  C1 <- survival_forest(X = tmpdat[tmpdat$INTENSIVE==1, 3:39],
                        Y = tmpdat[tmpdat$INTENSIVE==1, "t_cvds"], 
                        D = 1 - tmpdat[tmpdat$INTENSIVE==1, "cvd"],
                        failure.times = unique_times)
  cen1 <- predict(C1, tmpdat[, 3:39], failure.times = unique_times)$predictions
  
  C0 <- survival_forest(X = tmpdat[tmpdat$INTENSIVE==0, 3:39],
                        Y = tmpdat[tmpdat$INTENSIVE==0, "t_cvds"], 
                        D = 1- tmpdat[tmpdat$INTENSIVE==0, "cvd"],
                        failure.times = unique_times)
  cen0 <- predict(C0, tmpdat[, 3:39], failure.times = unique_times)$predictions
  
  Gamma_i <- compute_scores(treatment = tmpdat$INTENSIVE, 
                            study_time = tmpdat$t_cvds, 
                            status = tmpdat$cvd, 
                            propensity_score = rep(0.5, dim(tmpdat)[1]),
                            treated_survival_curve = surv1, 
                            treated_censoring_survival_curve = cen1,
                            control_survival_curve = surv0,
                            control_censoring_survival_curve = cen0,
                            unique_times = unique_times, 
                            end_time = t50, 
                            rmst = FALSE)#; mean(Gamma_i)
  
  est_ate[B] <- mean(Gamma_i)
  # sum(tmpdat$t_cvds[tmpdat$cvd==1]<=median(tmpdat$t_cvds[tmpdat$cvd==1]))  # check how close P(T>t0) to 0
  # sum(tmpdat$t_cvds[tmpdat$cvd==0]>median(tmpdat$t_cvds))                  # check how close P(C>t0) to 0 
  
  g <- 20
  N <- length(Gamma_i)
  
  # Scoxph calibration error 
  Scoxph_cate_sub <- Scoxph_cate[unique(index)]
  Scoxph_ate[B] <- mean(Scoxph_cate_sub)
  groups <- cut(Scoxph_cate_sub, breaks = unique(quantile(Scoxph_cate_sub, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_i <- rep(NA, N)
  for (k in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
    gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[k]   # leave-one-out gamma_hat estimates
  }
  theta_hat_scoxph[B] <- mean((Gamma_i-Scoxph_cate_sub)*(gamma_Delta_hat_i - Scoxph_cate_sub))
  
  # Slasso calibration error 
  Slasso_cate_sub <- Slasso_cate[unique(index)]
  Slasso_ate[B] <- mean(Slasso_cate_sub)
  groups <- cut(Slasso_cate_sub, breaks = unique(quantile(Slasso_cate_sub, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_i <- rep(NA, N)
  for (k in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
    gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[k]   # leave-one-out gamma_hat estimates
  }
  theta_hat_slasso[B] <- mean((Gamma_i-Slasso_cate_sub)*(gamma_Delta_hat_i - Slasso_cate_sub))
  
  
  # rlasso calibration error 
  rlasso_cate_sub <- rlasso_cate[unique(index)]
  rlasso_ate[B] <- mean(rlasso_cate_sub)
  groups <- cut(rlasso_cate_sub, breaks = unique(quantile(rlasso_cate_sub, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_i <- rep(NA, N)
  for (k in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
    gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[k]   # leave-one-out gamma_hat estimates
  }
  theta_hat_rlasso[B] <- mean((Gamma_i-rlasso_cate_sub)*(gamma_Delta_hat_i - rlasso_cate_sub))
  print(B)
}


# 95% bootstrapped CI 
Scoxph_CI <- quantile(theta_hat_scoxph, probs = c(0.025, 0.975))
Slasso_CI <- quantile(theta_hat_slasso, probs = c(0.025, 0.975))
Rlasso_CI <- quantile(theta_hat_rlasso, probs = c(0.025, 0.975))
aipw_ate_CI <- quantile(est_ate, probs = c(0.025, 0.975))
Scoxph_ate_CI <- quantile(Scoxph_ate, probs = c(0.025, 0.975))
Slasso_ate_CI <- quantile(Slasso_ate, probs = c(0.025, 0.975))
rlasso_ate_CI <- quantile(rlasso_ate, probs = c(0.025, 0.975))


Scoxph_est <- c(mean(theta_hat_scoxph), Scoxph_CI)
Slasso_est <- c(mean(theta_hat_slasso), Slasso_CI)
Rlasso_est <- c(mean(theta_hat_rlasso), Rlasso_CI)
aipw_ate_est <- c(mean(est_ate), aipw_ate_CI)
Scoxph_ate_est <- c(mean(Scoxph_ate), Scoxph_ate_CI)
Slasso_ate_est <- c(mean(Slasso_ate), Slasso_ate_CI)
rlasso_ate_est <- c(mean(rlasso_ate), rlasso_ate_CI)


out <- round(rbind(Scoxph_est, Slasso_est, Rlasso_est, aipw_ate_est, Scoxph_ate_est, Slasso_ate_est, rlasso_ate_est),4)
write.csv(out, "calibration_error.csv")


#                         2.5%  97.5%
# Scoxph_est     0.0022  0.0003 0.0044
# Slasso_est     0.0045  0.0019 0.0077
# Rlasso_est     0.0005 -0.0007 0.0016
# aipw_ate_est   0.0226  0.0083 0.0354
# Scoxph_ate_est 0.0213  0.0203 0.0225
# Slasso_ate_est 0.0646  0.0632 0.0662
# rlasso_ate_est 0.0166  0.0166 0.0167 

