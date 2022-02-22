# library(devtools)
# install_github("erikcs/grf@RATE", subdir = "r-package/grf")
#remove(list = ls())
library(readr)
library(dplyr)
library(grf)
#setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE/Data Application/Criteo")

data_loc = "criteo-uplift-v2.1.csv.gz"
seed = 20211013
ns = 160000
num.trees = 500
mtry = 12
honestly = F
min_node_size = 5000

data = read.csv(data_loc)
set.seed(seed)
subsamp = data %>% mutate(traintest = sample(0:1, n(), replace = TRUE)) %>% group_by(traintest, treatment) %>% sample_n(ns) %>% ungroup() %>% group_by(traintest) %>% group_split()
test_subsamp = subsamp[[2]]
subsamp = subsamp[[1]]
train_features = model.matrix(~.-conversion-treatment-exposure-visit-traintest, data=subsamp)

# save(test_subsamp, file="test_subsamp.Rdata")
load("test_subsamp.Rdata")

# print("Average visit:")
# print(mean(test_subsamp$visit))
# print(t.test(test_subsamp$visit)$"conf.int")
# print("Average conversion:")
# print(mean(test_subsamp$conversion))
# print(t.test(test_subsamp$conversion)$"conf.int")
# print("Per arm train")
# print(aggregate(subsamp[,c("visit", "conversion")], list(subsamp$treatment), mean))
# print("Per arm test")
# print(aggregate(test_subsamp[,c("visit", "conversion")], list(test_subsamp$treatment), mean))

test_features = model.matrix(~.-conversion-treatment-exposure-visit-traintest, data=test_subsamp)

# model = grf::regression_forest(X=train_features[subsamp$treatment==0, ], Y=subsamp$conversion[subsamp$treatment==0], 
#                                mtry=mtry, honesty=honestly, min.node.size=min_node_size)
# save(model, file="baseline_model.Rdata")
# load("baseline_model.Rdata")
# preds = predict(model, newdata=test_features)$predictions

### Causal forest model 
# Trained with 2000 versus 500 trees made not much difference on the CATE
# cate_model = grf::causal_forest(X=train_features, Y=subsamp$conversion, W=subsamp$treatment, W.hat=rep(mean(subsamp$treatment), ns*2),
#                                 num.trees = num.trees, mtry=mtry, honesty=honestly, min.node.size=min_node_size)
# save(cate_model, file="cate_model.Rdata")
load("cate_model.Rdata")
preds_cate = predict(cate_model, newdata=test_features)$predictions
# -0.000779816  0.000008152  0.000071896  0.001099639  0.000655435  0.011811772 

### S-learner of random forest 
# covs <- data.frame(train_features, subsamp$treatment); names(covs)[14] <- "treatment"
# s_rf_model <- grf::regression_forest(X = covs, Y = subsamp$conversion, num.trees = num.trees, 
#                                      mtry = mtry, honesty = honestly, min.node.size = min_node_size)
# save(s_rf_model, file="s_rf_model.Rdata")
load("s_rf_model.Rdata")
covs_test <- data.frame(test_features, test_subsamp$treatment); names(covs_test)[14] <- "treatment"
covs1 <- covs_test; covs1$treatment <- 1
covs0 <- covs_test; covs0$treatment <- 0
u1 <- predict(s_rf_model, covs1)$predictions 
u0 <- predict(s_rf_model, covs0)$predictions 
preds_srf_cate <- u1 - u0
# -0.0001333804 -0.0000005095  0.0000080070  0.0000145442  0.0000270746  0.0001925123 

# y_eval = test_subsamp$visit

# Evaluate the calibration error of causal forest estimates
full_test <- model.matrix(~.-exposure-visit-traintest, data=test_subsamp)
rownames(full_test) <- 1:dim(full_test)[1]
nboot <- 1000
theta_hat_cf <- theta_hat_srf <- Gamma_i_mean <- rep(NA, nboot)
for (B in 1:nboot){
  print(B)
  ## AIPW-based Gamma_i
  index <- sample(1:dim(full_test)[1], dim(full_test)[1], replace = TRUE)
  X_tmpdat <- full_test[unique(index), 1:(dim(full_test)[2]-2)]
  Y_tmpdat <- full_test[unique(index), "conversion"]
  W_tmpdat <- full_test[unique(index), "treatment"]
  
  # Out-of-bag mean outcome estimates
  Xs <- data.frame(X_tmpdat, W_tmpdat)
  Xs1 <- Xs;  Xs1$W_tmpdat <- 1
  Xs0 <- Xs;  Xs0$W_tmpdat <- 0
  N <- dim(Xs)[1]
  regfit <- grf::regression_forest(X = Xs, Y = Y_tmpdat, num.trees = num.trees, 
                                   mtry = mtry, honesty = honestly, min.node.size = min_node_size)
  u  <- predict(regfit)$predictions                    
  u1 <- predict(regfit, Xs1)$predictions 
  u0 <- predict(regfit, Xs0)$predictions 
  
  ps <- mean(W_tmpdat)
  Gamma_i <- unlist(u1 - u0 + (W_tmpdat - ps) / (ps * (1 - ps)) * (Y_tmpdat - u)) 
  Gamma_i_mean[B] <- mean(Gamma_i)
  
  ## Causal_forest calibration error
  g <- 20    # round(20*(N/500)^(2/5))
  preds_cate_sub <- preds_cate[unique(index)]
  groups <- cut(preds_cate_sub, breaks = unique(quantile(preds_cate_sub, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_cf_i <- (N/g)/(N/g-1)*gamma_Delta_hat[as.numeric(groups)] - 1/(N/g-1)*Gamma_i
  theta_hat_cf[B] <- mean((Gamma_i - preds_cate_sub)*(gamma_Delta_hat_cf_i - preds_cate_sub))
  
  pred_cate_cf <- aggregate(preds_cate_sub, by = list(groups), FUN = "mean")$x
  obs_cate_cf  <- aggregate(gamma_Delta_hat_cf_i, by = list(groups), FUN = "mean")$x
    
  ## S-learner RF calibration error
  preds_srfcate_sub <- preds_srf_cate[unique(index)]
  groups <- cut(preds_srfcate_sub, breaks = unique(quantile(preds_srfcate_sub, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_srf_i <- (N/g)/(N/g-1)*gamma_Delta_hat[as.numeric(groups)] - 1/(N/g-1)*Gamma_i
  theta_hat_srf[B] <- mean((Gamma_i - preds_srfcate_sub)*(gamma_Delta_hat_srf_i - preds_srfcate_sub))
  
  pred_cate_srf <- aggregate(preds_srfcate_sub, by = list(groups), FUN = "mean")$x
  obs_cate_srf  <- aggregate(gamma_Delta_hat_srf_i, by = list(groups), FUN = "mean")$x
}

cf_CI <- quantile(theta_hat_cf, probs = c(0.025, 0.975), na.rm=T)
Srf_CI <- quantile(theta_hat_srf, probs = c(0.025, 0.975), na.rm=T)

cf_est <- c(mean(theta_hat_cf, na.rm=T), cf_CI)
Srf_est <- c(mean(theta_hat_srf, na.rm=T), Srf_CI)

ate_aipw_CI <- quantile(Gamma_i_mean, probs = c(0.025, 0.975))
ate_est <- c(mean(Gamma_i_mean), ate_aipw_CI)

out <- rbind(cf_est, Srf_est, ate_est)
write.csv(out, "calibration_error_criteo_20bins.csv")


# # calibration plot for causal forest
# cf_data <- data.frame(pred_cate_cf, obs_cate_cf)
# cf <- ggplot(cf_data, aes(x = pred_cate_cf, y = obs_cate_cf)) +
#   theme_bw()+
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-0.001, 0.014)+
#   ylim(-0.001, 0.014)+
#   xlab("Predicted CATE using causal forest")+
#   ylab("Calibration function estimates")+
#   geom_abline(intercept = 0, slope = 1, col = 'red') +
#   geom_point(); cf
# 
# # calibration plot for S-learner RF
# srf_data <- data.frame(pred_cate_srf, obs_cate_srf)
# srf <- ggplot(srf_data, aes(x = pred_cate_srf, y = obs_cate_srf)) +
#   theme_bw() + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#   xlim(-0.001, 0.014)+
#   ylim(-0.001, 0.014)+
#   xlab("Predicted CATE using S-learner of random forest")+
#   ylab("Calibration function estimates")+
#   geom_abline(intercept = 0, slope = 1, col = 'red') +
#   geom_point(); srf
# 
# library(cowplot)
# png("criteo_calib_plot_10bins.png", width = 4, height = 6, units = 'in', res = 300)
# plot_grid(cf, srf, labels=c("A", "B"), ncol = 1, nrow = 2)
# dev.off()

