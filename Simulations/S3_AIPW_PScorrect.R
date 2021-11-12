
## Setting 3: Replicate the simulations in Setting 2 using an AIPW estimator 
library(grf)
parameter <- read.csv("./parameter.csv")

if(length(args <- commandArgs(T))>0){
  caseid <- as.integer(args[[1]])
  message("running for parameter set ", caseid)
}

N <- parameter[caseid, 1]
alpha <- parameter[caseid, 2]
B <- 1000
out <- matrix(NA, B, 2)
for (j in 1:B){

  X0 <- rnorm(N, 0, 1)
  logitp <- 0.3*X0
  p <- exp(logitp)/(1+exp(logitp))
  W <- rbinom(N, 1, p)
  Delta <- 0.5*X0     # Delta ~ N(0, 0.25)
  
  gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2
  Y0 <- rnorm(N, 0, 1)
  Y1 <- Y0 + gamma_Delta
  
  # Observed outcome 
  Y <- W*Y1 + (1-W)*Y0
  
  # Estimated calibration curve using IPW 
  g <- 20
  groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
  
  # Estimated propensity score (correct model)
  fitps <- glm(W ~ X0, family = binomial)
  ps <- predict(fitps, data.frame(X0), type = "response")
  
  # Estimated PS with a wrong model 
  # fitps <- glm(W ~ 1, family = binomial)
  # ps <- predict(fitps, data.frame(X0), type = "response")
  
  # Use an AIPWE instead of IPWE
  Xs <- data.frame(W, X0)
  Xs1 <- Xs
  Xs1$W <- 1
  Xs0 <- Xs
  Xs0$W <- 0
  
  regfit <- regression_forest(X = Xs, Y = Y, num.trees = 500, num.threads = 4)
  u  <- predict(regfit)$predictions                   # out-of-bag estimates 
  u1 <- predict(regfit, Xs1)
  u0 <- predict(regfit, Xs0)
  
  Gamma_i <- u1 - u0 + (W-ps)/(ps*(1-ps))*(Y-u)   # AIPWE
  gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$predictions
  gamma_Delta_hat_i <- rep(NA, N)
  for (k in 1:N){
    gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
    gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*unlist(Gamma_i)[k]   # leave-one-out gamma_hat estimates
  }
  
  # True calibration error 
  integrand <- function(x) {dnorm(x, 0, 0.5)*x^2*(1-x)^2}    
  integral <- integrate(integrand, -Inf, Inf)$value
  theta <- alpha^2*integral
  
  # Estimated calibration error - L2
  theta_hat <- mean((gamma_Delta_hat_i - Delta)^2)
  
  # Robust estimates
  theta_hat_robust <- mean((as.numeric(unlist(Gamma_i))-Delta)*(gamma_Delta_hat_i - Delta))
  
  # Bias 
  out[j,1] <- theta_hat - theta
  out[j,2] <- theta_hat_robust - theta 
  print(paste0("bootstrap #", j))
}
out <- data.frame(caseid, out)
write.csv(out, paste0(caseid,"setting3_pscorrect_output.csv"))
