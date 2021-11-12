
## RCT setting + HD covariates + compare SE using an AIPW estimator  
library(grf)
setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")
parameter <- read.csv("./parameter_ndim.csv")

if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}

N <- parameter[caseid, 1]
alpha <- parameter[caseid, 2]
ndim <- parameter[caseid, 3]
B <- 1000
nboot <- 500
out <- matrix(NA, B, 6)
for (j in 1:B){
  
  # independent covariates (p=100)
  X_ind <- matrix(NA, N, ndim)
  for (k in 1:ndim){
    X_ind[,k] <- rnorm(N, 0, 1)
  }
  
  X1 <- rnorm(N, 0, 1)
  W <- rbinom(N, 1, 0.5)
  Delta <- runif(N, -1, 1)
  gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2
  Y0 <- X1 + rnorm(N,0,1)
  Y1 <- Y0 + gamma_Delta
  
  # Observed outcome 
  Y <- W*Y1 + (1-W)*Y0
  
  fulldata <- data.frame(X_ind, W, X1, Delta, Y)
  
  # Estimated calibration curve using IPW 
  g <- 20
  groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
  ps <- 0.5
  Xs <- data.frame(W, X1, X_ind)
  Xs1 <- Xs; Xs1$W <- 1
  Xs0 <- Xs; Xs0$W <- 0
  regfit <- regression_forest(X = Xs, Y = Y, num.trees = 500, num.threads = 4)
  u  <- predict(regfit)$predictions                   # out-of-bag estimates 
  u1 <- predict(regfit, Xs1)
  u0 <- predict(regfit, Xs0)
  Gamma_i <- unlist(u1 - u0 + (W-ps)/(ps*(1-ps))*(Y-u))   # AIPWE
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_i <- rep(NA, N)
  for (h in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[h]]
    gamma_Delta_hat_i[h] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[h]   # leave-one-out gamma_hat estimates
  }
  
  # True calibration error 
  integrand <- function(x) {0.5*x^2*(1-x)^2}
  integral <- integrate(integrand, -1, 1)$value
  theta <- alpha^2*integral
  
  # Estimated calibration error 
  theta_hat <- mean((gamma_Delta_hat_i - Delta)^2)
  
  # Robust estimates
  theta_hat_robust <- mean((Gamma_i-Delta)*(gamma_Delta_hat_i - Delta))
  
  # Bias 
  out[j,1] <- theta_hat - theta
  out[j,2] <- theta_hat_robust - theta 
  print(paste0("bootstrap #", j))
  
  
  #--------------------------------------
  # Bootstrapped calibration error
  bootout <- matrix(NA, nboot, 2)
  for (z in 1:nboot){
    bootindex <- sample(1:N, N, replace = T)
    bootdata <- fulldata[unique(bootindex),]
    bootdata$ps <- 0.5
    
    # Calculate Gamma_i using cross-fitting estimates
    Xs <- bootdata[,1:(ndim+2)]
    Xs1 <- Xs;  Xs1$W <- 1
    Xs0 <- Xs;  Xs0$W <- 0
    regfit <- regression_forest(X = Xs, Y = bootdata$Y, num.trees = 500, num.threads = 4)
    u  <- predict(regfit)$predictions                   # out-of-bag estimates 
    u1 <- predict(regfit, Xs1)
    u0 <- predict(regfit, Xs0)
    Gamma_i <- unlist(u1 - u0 + (bootdata$W-bootdata$ps)/(bootdata$ps*(1-bootdata$ps))*(bootdata$Y-u))   # AIPWE
    
    # Estimated calibration curve using IPW 
    g <- 20
    groups <- cut(bootdata$Delta, breaks = unique(quantile(bootdata$Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
    gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
    gamma_Delta_hat_i <- rep(NA, dim(bootdata)[1])
    for (i in 1:dim(bootdata)[1]){
      gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[i]]
      gamma_Delta_hat_i[i] <- (dim(bootdata)[1]/g)/(dim(bootdata)[1]/g-1)*gamma_hat_i - 1/(dim(bootdata)[1]/g-1)*unlist(Gamma_i)[i]   # leave-one-out gamma_hat estimates
    }
    
    # Estimated calibration error - L2
    theta_hat <- mean((gamma_Delta_hat_i - bootdata$Delta)^2)
    
    # Robust estimates
    theta_hat_robust <- mean((Gamma_i-bootdata$Delta)*(gamma_Delta_hat_i - bootdata$Delta))
    
    # Bias 
    bootout[z,1] <- theta_hat - theta
    bootout[z,2] <- theta_hat_robust - theta 
  }
  out[j,3] <- mean(bootout[,1])
  out[j,4] <- mean(bootout[,2])
  out[j,5] <- sqrt(var(bootout[,1]))
  out[j,6] <- sqrt(var(bootout[,2]))
}
out <- data.frame(caseid, out)
write.csv(out, paste0(caseid,"S1B_RCT_psC_BootSE.csv"))
