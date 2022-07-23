# Estimated Calibration Error for Treatment Heterogeneity (ECETH)
This repository contains the implementation code for paper entitled "Calibration Error for Heterogeneous Treatment Effect". In this work, we define an analogue of the L2 expected calibration error for Heterogeneous Treatment Effects (HTEs) and propose a robust estimator. Our approach is motivated by doubly robust treatment effect estimators, making it unbiased, and resilient to confounding, overfitting, and high-dimensionality issues. Furthermore, our method is straightforward to adapt to many structures under which treatment effects can be identified, including randomized trials, observational studies, and survival analysis.

Some links for getting started

* The [R](https://github.com/CrystalXuR/ECETH/tree/master/R) folder contains all the functions needed to apply the ECETH metric
* The replication code of the our simulation study for showing the validity of the ECETH metric is available in [experiments](https://github.com/CrystalXuR/ECETH/tree/master/experiemnts)

### Usage Examples

The following script demonstrates the application of the ECETH metric for evaluating the estimated HTE on a binary and time-to-event outcome separately. The calibration curve is estimated using an augmented-inverse-probability-weighted (AIPW) score with a histogram estimator, and a generalized AIPW estimator is used under the survival setting.

```R
source("ECETH.R")

# Generate data with a binary outcome
p <- 5
X <- matrix(runif(n * p), n, p)
W <- rbinom(n, 1, 0.5)
event.prob <- 1 / (1 + exp(2*(pmax(2*X[, 1], 0) * W - X[, 2])))
Y <- rbinom(n, 1, event.prob)
train <- sample(1:n, n / 2)

# Fit a causal forest for estimating HTEs
cf.fit <- causal_forest(X[train, ], Y[train], W[train])
cf.cate <- predict(cf.fit, X[-train, ])$predictions

# Compute the AIPW score
source("scores.R")
cate.score <- aipw_score(X = X[-train, ], W = W[-train], Y = Y[-train])

# ECETH
caliberror <- eceth(cf.cate, cate.score)


# Generate data with a time-to-event outcome
horizon <- 1
failure.time <- pmin(rexp(n) * X[, 1] + W, horizon)
censor.time <- 2 * runif(n)
Y <- round(pmin(failure.time, censor.time), 2)
D <- as.integer(failure.time <= censor.time)
t0 <- median(Y)

# Fit a causal survival forest for estimating HTEs
source("survival_aipw_score.R")
cs.forest <- causal_survival_forest(X[train, ], Y[train], W[train], D[train], target = "survival.probability", horizon = t0)
csf.cate <- predict(cs.forest, X[-train, ])$predictions

# Compute the AIPW score for survival outcomes
cate.surv.score <- aipw_surv_score(X = X[-train, ], W = W[-train], Y = Y[-train], D = D[-train], t0 = t0)

# ECETH
caliberror <- eceth(csf.cate, cate.surv.score)
```

### Funding

This work was supported by R01 HL144555 from the National Heart, Lung, and Blood Institute (NHLBI).

### References
Yizhe Xu, Steve Yadlowsky (2022). **Calibration Error for Heterogeneous Treatment Effect** Proceedings of the 25th International Conference on Artificial Intelligence and Statistics (AISTATS) 2022, Valencia, Spain. PMLR: Volume 151. [<a href="https://proceedings.mlr.press/v151/xu22c.html">paper</a>]

Steve Yadlowsky, Scott Fleming, Nigam Shah, Emma Brunskill, and Stefan Wager. **Evaluating treatment prioritization rules via rank-weighted average treatment effects.** arXiv preprint arXiv:2111.07966, 2021. [<a href="https://arxiv.org/abs/2111.07966">arXiv</a>]

James M Robins, Andrea Rotnitzky, and Lue Ping Zhao. **Estimation of regression coefficients when some regressors are not always observed.** Journal of the American statistical Association, 89(427): 846–866, 1994. [<a href="https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476818">paper</a>]
