
## Setting 3: Replicate the simulations in Setting 2 using an AIPW estimator 
library(grf)
parameter <- read.csv("./parameter_ndim.csv")

if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}
# caseid <- 26
N <- parameter[caseid, 1]
alpha <- parameter[caseid, 2]
ndim <- parameter[caseid, 3]
B <- 10
nboot <- 500
out <- matrix(NA, B, 6)
for (j in 1:B){
  
  # independent covariates (p=100)
  X_ind <- matrix(NA, N, ndim)
  for (k in 1:ndim){
    X_ind[,k] <- rnorm(N, 0, 1)
  }
  
  X1 <- rnorm(N, 0, 1)
  X0 <- rnorm(N, 0, 1)
  logitp <- 0.3*X0
  p <- exp(logitp)/(1+exp(logitp))
  W <- rbinom(N, 1, p)
  Delta <- 0.5*X0     # Delta ~ N(0, 0.25)
  
  gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2
  Y0 <- X1 + rnorm(N, 0, 1)
  Y1 <- Y0 + gamma_Delta
  
  # Observed outcome 
  Y <- W*Y1 + (1-W)*Y0
  
  # Split data into 2 halves
  fulldata <- data.frame(X0, X_ind, W, X1, Delta, Y)
  index <- sample(1:N, N/2, replace = F)
  dat1 <- fulldata[index,]
  dat2 <- fulldata[-as.numeric(index),]
  
  # Out-of-bag PS estimates using split samples
  tmpdat1 <- data.frame(dat1[,1:(ndim+2)])
  tmpdat2 <- data.frame(dat2[,1:(ndim+2)])
  
  fitps1 <- glm(W ~ ., family = binomial, data=tmpdat1)
  ps1 <- predict(fitps1, newdata = tmpdat2[,1:(ndim+1)], type = "response")
  
  fitps2 <- glm(W ~ ., family = binomial, data = tmpdat2)
  ps2 <- predict(fitps2, newdata = tmpdat1[,1:(ndim+1)], type = "response")
  
  fulldata$ps <-rep(NA, dim(fulldata)[1])
  fulldata$ps[-as.numeric(index)] <- ps1
  fulldata$ps[index] <- ps2 
  
  # Calculate Gamma_i using cross-fitting estimates
  Xs <- data.frame(W, X0, X1, X_ind)
  Xs1 <- Xs;  Xs1$W <- 1
  Xs0 <- Xs;  Xs0$W <- 0
  regfit <- regression_forest(X = Xs, Y = Y, num.trees = 500, num.threads = 4)
  u  <- predict(regfit)$predictions                   # out-of-bag estimates 
  u1 <- predict(regfit, Xs1)
  u0 <- predict(regfit, Xs0)
  Gamma_i <- unlist(u1 - u0 + (W-fulldata$ps)/(fulldata$ps*(1-fulldata$ps))*(Y-u))   # AIPWE
  
  # Estimated calibration curve using IPW 
  g <- 20
  groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
  gamma_Delta_hat_i <- rep(NA, N)
  for (h in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[h]]
    gamma_Delta_hat_i[h] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*unlist(Gamma_i)[h]   # leave-one-out gamma_hat estimates
  }
  
  # True calibration error 
  integrand <- function(x) {dnorm(x, 0, 0.5)*x^2*(1-x)^2}    
  integral <- integrate(integrand, -Inf, Inf)$value
  theta <- alpha^2*integral
  
  # Estimated calibration error - L2
  theta_hat <- mean((gamma_Delta_hat_i - Delta)^2)
  
  # Robust estimates
  theta_hat_robust <- mean((Gamma_i-Delta)*(gamma_Delta_hat_i - Delta))
  
  # Bias 
  out[j,1] <- theta_hat - theta
  out[j,2] <- theta_hat_robust - theta 
  print(paste0("replicate #", j))
  
  #--------------------------------------
  # Bootstrapped calibration error
  bootout <- matrix(NA, nboot, 2)
  for (z in 1:nboot){
    bootindex <- sample(1:N, N, replace = T)
    bootdata <- fulldata[unique(bootindex),]
    rownames(bootdata) <- 1:dim(bootdata)[1]
    index <- sample(1:dim(bootdata)[1], floor(dim(bootdata)[1]/2), replace = F)
    dat1 <- bootdata[index,]
    dat2 <- bootdata[-as.numeric(index),]
    
    # Out-of-bag PS estimates using split samples
    tmpdat1 <- data.frame(dat1[,1:(ndim+2)])
    tmpdat2 <- data.frame(dat2[,1:(ndim+2)])
    
    fitps1 <- glm(W ~ ., family = binomial, data=tmpdat1)
    ps1 <- predict(fitps1, newdata = tmpdat2[,1:(ndim+1)], type = "response")
    
    fitps2 <- glm(W ~ ., family = binomial, data = tmpdat2)
    ps2 <- predict(fitps2, newdata = tmpdat1[,1:(ndim+1)], type = "response")
    
    bootdata$ps <-rep(NA, dim(bootdata)[1])
    bootdata$ps[-as.numeric(index)] <- ps1
    bootdata$ps[index] <- ps2 
    
    # Calculate Gamma_i using cross-fitting estimates
    Xs <- data.frame(bootdata$W, bootdata$X0, bootdata$X1, bootdata[,2:101]);names(Xs)[1:3]<-c("W","X0","X1")
    Xs1 <- Xs;  Xs1$W <- 1
    Xs0 <- Xs;  Xs0$W <- 0
    regfit <- regression_forest(X = Xs, Y = bootdata$Y, num.trees = 500, num.threads = 4);regfit
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
out <- data.frame(caseid, out); names(out) <- c("caseid","ce","robust_ce","boot_ce","boot_robust_ce","boot_se","boot_robust_se")
write.csv(out, paste0(caseid,"S3_psC_split_HD_BootSE.csv"))
