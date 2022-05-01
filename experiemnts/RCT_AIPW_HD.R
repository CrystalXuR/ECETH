
## RCT setting + HD covariates using an AIPW estimator 
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
out <- matrix(NA, B, 4)
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
  
  # Estimated calibration curve using IPW 
  ps <- mean(W)
  Xs <- data.frame(W, X1, X_ind)
  Xs1 <- Xs; Xs1$W <- 1
  Xs0 <- Xs; Xs0$W <- 0
  regfit <- regression_forest(X = Xs, Y = Y, num.trees = 200, num.threads = 4,
                              honesty=F, min.node.size=(N/65)^(4/5))
  u  <- predict(regfit)$predictions                   
  u1 <- predict(regfit, Xs1)
  u0 <- predict(regfit, Xs0)
  Gamma_i <- unlist(u1 - u0 + (W-ps)/(ps*(1-ps))*(Y-u))   # AIPWE
  
  g <- round(20*(N/500)^(2/5))
  groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
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
  
  # Estimates
  out[j,3] <- theta_hat 
  out[j,4] <- theta_hat_robust
  
  print(paste0("bootstrap #", j))
}
out <- data.frame(caseid, out)
write.csv(out, paste0(caseid,"RCT_AIPW_HD.csv"))
