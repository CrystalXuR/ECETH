#----------------------- HTE calibration paper  ----------------------#
#--------------------- Yizhe Xu & Steve Yadlowsky --------------------#
#---------------------------------------------------------------------#

setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")

### Setting 1: a randomized trial with a continuous outcome 
set.seed(123)
N <- 10000
W <- rbinom(N, 1, 0.5)
Delta <- runif(N, -1, 1)

# Scenario1: alpha = 0 -- HTE prediction model is well calibrated -- bias = 0.003
# Scenario2: alpha = 1 -- HTE prediction model is miscalibrated   -- bias = -0.49
alpha <- 0
gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2
Y0 <- rnorm(N, 0, 1)
Y1 <- Y0 + gamma_Delta

# Observed outcome 
Y <- W*Y1 + (1-W)*Y0

# Estimated calibration curve using IPW 
g <- 10
deciles <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
p <- 0.5
gamma_Delta_hat  <- aggregate(W*Y/p - (1-W)*Y/(1-p), by = list(deciles), FUN = "mean")$x

# Estimated calibration error 
agg_Delta  <- aggregate(Delta, by = list(deciles), FUN = "mean")$x
theta_hat <- mean((gamma_Delta_hat - agg_Delta)^2)

# True calibration error 
integrand <- function(x) {x^2*(1-x)^2}
integral <- integrate(integrand, -1, 1)$value
theta <- alpha^2*integral

# Bias 
bias <- theta_hat - theta; bias


### Setting 2: an observational study with a continuous outcome 
set.seed(123)
N <- 10000
X0 <- rnorm(N, 1, 2)
logitp <- 0.3*X0
p <- exp(logitp)/(1+exp(logitp))
W <- rbinom(N, 1, p)
Delta <- rnorm(N, -1, 1)

# Scenario 1: alpha = 0 & ps is correctly modeled     -- bias = 0.0039
# Scenario 2: alpha = 0 & ps is NOT correctly modeled -- bias = 0.051
# Scenario 3: alpha = 1 & ps is correctly modeled     -- bias = 18.72
# Scenario 4: alpha = 1 & ps is NOT correctly modeled -- bias = 22.33
alpha <- 1
gamma_Delta <- (1-alpha)*Delta + alpha*Delta^2
Y0 <- rnorm(N, 0, 1)
Y1 <- Y0 + gamma_Delta 

# Observed outcome 
Y <- W*Y1 + (1-W)*Y0

# Estimated calibration curve using IPW 
g <- 10
deciles <- cut(Delta, breaks = unique(quantile(Delta, probs=seq(0,1,1/g))), include.lowest=TRUE)

# Estimated propensity score (correct model)
fitps <- glm(W ~ X0, family = binomial)
p <- predict(fitps, data.frame(X0), type = "response")  

# Estimated PS with a wrong model 
# p <- 0.5   

gamma_Delta_hat  <- aggregate(W*Y/p - (1-W)*Y/(1-p), by = list(deciles), FUN = "mean")$x

# Estimated calibration error 
agg_Delta  <- aggregate(Delta, by = list(deciles), FUN = "mean")$x
theta_hat <- mean((gamma_Delta_hat - agg_Delta)^2)

# True calibration error 
integrand <- function(x) {x^2*(1-x)^2}
integral <- integrate(integrand, -1, 1)$value
theta <- alpha^2*integral

# Bias 
bias <- theta_hat - theta; bias

# Check if sharing nuisance parameter (PS) in the robust formula has an impact in terms of bias
# Calibration curve in robust formula? 





