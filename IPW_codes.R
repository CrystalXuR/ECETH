#----------------------- HTE calibration paper  ----------------------#
#--------------------- Yizhe Xu & Steve Yadlowsky --------------------#
#---------------------------------------------------------------------#

setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")

### Setting 1A: a randomized trial with a continuous outcome 
B <- 1000
Ns <- c(500, 1000, 2000, 4000)
alphas <- c(0, 0.15, 0.3)
bias <- bias_robust <- theta_plugin <- theta_robust <- array(NA, dim=c(B, length(alphas),length(Ns)), 
                             dimnames = list(1:B, paste0("alpha = ",alphas), paste0("N = ",Ns)))
for (h in 1:length(Ns)){
  
  N <- Ns[h]
  
  for (i in 1:length(alphas)){
    
    alpha <- alphas[i]
    print(paste0("N=",N," alpha=", alpha))
    for (j in 1:B){
      
      X1 <- rnorm(N, 0, 1)
      W <- rbinom(N, 1, 0.5)
      Delta <- runif(N, -1, 1)  # predicted ARR
      gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2    # true ARR
      Y0 <- X1 + rnorm(N,0,1)
      Y1 <- Y0 + gamma_Delta
      
      # Observed outcome 
      Y <- W*Y1 + (1-W)*Y0
      
      # Estimated calibration curve using IPW 
      g <- round(20*(N/500)^(2/5))
      groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
      ps <- mean(W)
      Gamma_i <- W*Y/ps - (1-W)*Y/(1-ps)
      gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
      gamma_Delta_hat_i <- rep(NA, N)
      for (k in 1:N){
        gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
        gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[k]   # leave-one-out gamma_hat estimates
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
      bias[j,i,h] <- theta_hat - theta
      bias_robust[j,i,h] <- theta_hat_robust - theta
      
      # Estimates
      theta_plugin[j,i,h] <- theta_hat
      theta_robust[j,i,h] <- theta_hat_robust
    }
  }
}
out <- list(NA, length(Ns))
for (z in 1:length(Ns)){
  out[[z]] <- data.frame(rbind(colMeans(bias[,,z]), colMeans(bias_robust[,,z]),
                               apply(theta_plugin[,,z],2,sd), apply(theta_robust[,,z],2,sd)))
  colnames(out[[z]]) <- paste0("alpha = ",alphas)
  rownames(out[[z]]) <- c("L2 bias", "Robust bias", "L2 se", "Robust se")
}
names(out) <- paste0("N = ",Ns); out
outlong <- data.frame(c(rep(500, 3),rep(1000, 3),rep(2000, 3),rep(4000, 3)), 
                      rep(c(0, 0.15,0.3),4),
                      rbind(t(out[[1]]), t(out[[2]]), t(out[[3]]), t(out[[4]])))
names(outlong)[1:2] <- c("N", "alpha")
outlong$L2.std.bias <- outlong$L2.bias/outlong$L2.se
outlong$Robust.std.bias <- outlong$Robust.bias/outlong$Robust.se
outlong$L2.MSE <- outlong$L2.bias^2 + outlong$L2.se^2
outlong$Robust.MSE <- outlong$Robust.bias^2 + outlong$Robust.se^2
outlong <- round(outlong, 4)
write.csv(outlong, "RCT_IPW.csv")


### Setting 2C: an observational study with a continuous outcome 
### Added: 1. A risk factor X1 that only predicts the outcome Y
###        2. split the data into 2 halves, with one half to estimate the PS used in Gamma_i,
###        and the other half to estimate the PS used in the calibration function
# Scenario 1: ps is correctly modeled     
# Scenario 2: ps is NOT correctly modeled
B <- 1000
Ns <- c(500, 1000, 2000, 4000)
alphas <- c(0, 0.15, 0.3)
bias <- bias_robust <- theta_plugin <- theta_robust <- array(NA, dim=c(B, length(alphas),length(Ns)), 
                             dimnames = list(1:B, paste0("alpha = ",alphas), paste0("N = ",Ns)))
for (h in 1:length(Ns)){
  
  N <- Ns[h]
  
  for (i in 1:length(alphas)){
    
    alpha <- alphas[i]
    print(paste0("N=", N, " alpha=", alpha))
    for (j in 1:B){
      
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
      Y <- W*Y1 + (1-W)*Y0     # Y = Y0 + W*X0*beta
      
      # Split data into 2 halves
      fulldata <- data.frame(X0, X1, W, Delta, Y)
      index <- sample(1:N, N/2, replace = F)
      dat1 <- fulldata[index,]
      dat2 <- fulldata[-as.numeric(index),]
      
      # Out-of-bag PS estimates using split samples
      # fitps1 <- glm(W ~ X0, family = binomial, data=dat1)
      # ps1 <- predict(fitps1, newdata = data.frame(X0=dat2$X0), type = "response")
      # 
      # fitps2 <- glm(W ~ X0, family = binomial, data = dat2)
      # ps2 <- predict(fitps2, newdata = data.frame(X0=dat1$X0), type = "response")
      
      # Estimated PS with a wrong model
      fitps1 <- glm(W ~ 1, family = binomial, data=dat1)
      ps1 <- predict(fitps1, newdata = data.frame(X0=dat2$X0), type = "response")

      fitps2 <- glm(W ~ 1, family = binomial, data = dat2)
      ps2 <- predict(fitps2, newdata = data.frame(X0=dat1$X0), type = "response")
      
      fulldata$ps <-rep(NA, dim(fulldata)[1])
      fulldata$ps[-as.numeric(index)] <- ps1
      fulldata$ps[index] <- ps2 
      
      # Calculate Gamma_i using cross-fitting estimates
      Gamma_i <- W*Y/fulldata$ps - (1-W)*Y/(1-fulldata$ps)
      
      # Estimated calibration curve using IPW 
      g <- round(20*(N/500)^(2/5))
      groups <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
      gamma_Delta_hat  <- aggregate(Gamma_i, by = list(groups), FUN = "mean")$x
      gamma_Delta_hat_i <- rep(NA, N)
      for (k in 1:N){
        gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
        gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*Gamma_i[k]   # leave-one-out gamma_hat estimates
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
      bias[j,i,h] <- theta_hat - theta
      bias_robust[j,i,h] <- theta_hat_robust - theta 
      
      # Estimates 
      theta_plugin[j,i,h] <- theta_hat
      theta_robust[j,i,h] <- theta_hat_robust
    }
  }
}
out <- list(NA, length(Ns))
for (z in 1:length(Ns)){
  out[[z]] <- data.frame(rbind(colMeans(bias[,,z]), colMeans(bias_robust[,,z]), 
                               apply(theta_plugin[,,z],2,sd), apply(theta_robust[,,z],2,sd)))
  colnames(out[[z]]) <- paste0("alpha = ",alphas)
  rownames(out[[z]]) <- c("L2 bias", "Robust bias", "L2 se", "Robust se")
}
names(out) <- paste0("N = ",Ns); out
outlong <- data.frame(c(rep(500, 3),rep(1000, 3),rep(2000, 3),rep(4000, 3)), 
                      rep(c(0, 0.15,0.3),4),
                      rbind(t(out[[1]]), t(out[[2]]), t(out[[3]]), t(out[[4]])))
names(outlong)[1:2] <- c("N", "alpha")
outlong$L2.std.bias <- outlong$L2.bias/outlong$L2.se
outlong$Robust.std.bias <- outlong$Robust.bias/outlong$Robust.se
outlong$L2.MSE <- outlong$L2.bias^2 + outlong$L2.se^2
outlong$Robust.MSE <- outlong$Robust.bias^2 + outlong$Robust.se^2
outlong <- round(outlong, 4)
write.csv(outlong, "obs_psW_IPW.csv")
# plot(Delta, gamma_Delta, col='red')
# points(Delta, gamma_Delta_hat_i)


