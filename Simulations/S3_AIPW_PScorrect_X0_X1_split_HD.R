
## Setting 3: Replicate the simulations in Setting 2 using an AIPW estimator 
library(grf)
parameter <- read.csv("./parameter_ndim.csv")

if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}

N <- parameter[caseid, 1]
alpha <- parameter[caseid, 2]
ndim <- parameter[caseid, 3]
B <- 1000
out <- rep(NA, B,2)
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
  
  # # Split data into 2 halves
  # fulldata <- data.frame(X0, X_ind, W, X1, Delta, Y)
  # index <- sample(1:N, N/2, replace = F)
  # dat1 <- fulldata[index,]
  # dat2 <- fulldata[-as.numeric(index),]
  # 
  # # Out-of-bag PS estimates using split samples
  # tmpdat1 <- data.frame(dat1[,1:(ndim+2)])
  # tmpdat2 <- data.frame(dat2[,1:(ndim+2)])
  # 
  # fitps1 <- glm(W ~ ., family = binomial, data=tmpdat1)
  # ps1 <- predict(fitps1, newdata = tmpdat2[,1:(ndim+1)], type = "response")
  # 
  # fitps2 <- glm(W ~ ., family = binomial, data = tmpdat2)
  # ps2 <- predict(fitps2, newdata = tmpdat1[,1:(ndim+1)], type = "response")
  # 
  # fulldata$ps <-rep(NA, dim(fulldata)[1])
  # fulldata$ps[-as.numeric(index)] <- ps1
  # fulldata$ps[index] <- ps2 
  
  psfit <- regression_forest(X = data.frame(X0, X1, X_ind), Y = W, num.trees = 500, num.threads = 4)
  ps  <- predict(psfit)$predictions
  
  # Calculate Gamma_i using cross-fitting estimates
  Xs <- data.frame(W, X0, X1, X_ind)
  Xs1 <- Xs;  Xs1$W <- 1
  Xs0 <- Xs;  Xs0$W <- 0
  regfit <- regression_forest(X = Xs, Y = Y, num.trees = 500, num.threads = 4)
  u  <- predict(regfit)$predictions                   # out-of-bag estimates 
  u1 <- predict(regfit, Xs1)
  u0 <- predict(regfit, Xs0)
  Gamma_i <- unlist(u1 - u0 + (W-ps)/(ps*(1-ps))*(Y-u))   # AIPWE
  
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
  print(paste0("bootstrap #", j))
}
out <- data.frame(caseid, out)
write.csv(out, paste0(caseid,"S3_psC_forest_HD.csv"))
