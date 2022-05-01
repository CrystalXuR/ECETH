#' @title Expected calibration error for treatment heterogeneity (ECETH)
#'
#' @description  Evaluate model calibration performance on estimating conditional average treatment effects
#'
#' @param cate.est A vector of estimated conditional average treatment effects
#' @param cate.score A vector of score that are surrogates of individual treatment effects
#' @param g The number of bins used to estimate calibration function. Default is 10
#' @examples
#' \donttest{
#' n <- 1500
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#'
#' # binary outcome
#' event.prob <- 1 / (1 + exp(2*(pmax(2*X[, 1], 0) * W - X[, 2])))
#' Y <- rbinom(n, 1, event.prob)
#' train <- sample(1:n, n / 2)
#' cf.fit <- causal_forest(X[train, ], Y[train], W[train])
#' cf.cate <- predict(cf.fit, X[-train, ])$predictions
#' # Compute the AIPW score
#' cate.score <- aipw_score(X = X[-train, ], W = W[-train], Y = Y[-train])
#' # ECETH
#' caliberror <- eceth(cf.cate, cate.score)
#'
#' # time-to-event outcome
#' horizon <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, horizon)
#' censor.time <- 2 * runif(n)
#' # Discretizing continuous events decreases runtime.
#' Y <- round(pmin(failure.time, censor.time), 2)
#' D <- as.integer(failure.time <= censor.time)
#' t0 <- median(Y)
#' # Compute the AIPW score
#' cs.forest <- causal_survival_forest(X[train, ], Y[train], W[train], D[train], target = "survival.probability", horizon = t0)
#' csf.cate <- predict(cs.forest, X[-train, ])$predictions
#' # Compute the AIPW score for survival outcomes
#' cate.surv.score <- aipw_surv_score(X = X[-train, ], W = W[-train], Y = Y[-train], D = D[-train], t0 = t0)
#' # ECETH
#' caliberror <- eceth(csf.cate, cate.surv.score)
#' }
#' @return
#' @export
eceth <- function(cate.est, cate.score, g = 10){
  if (length(unique(cate.est)) == 1) {
    warning("CATE predictions are identical, expected squared error of average treatment effect is computed instead")
    error <- mean((cate.score - cate.est)^2)
  } else {
    groups <- cut(cate.est, breaks = unique(quantile(cate.est, probs=seq(0,1,1/g))), include.lowest=TRUE)
    gamma_Delta_hat  <- aggregate(cate.score, by = list(groups), FUN = "mean")$x
    if (length(gamma_Delta_hat) < g){
      warning(paste0(length(gamma_Delta_hat), " groups are used instead of ", g, " due to not enough distinct values"))
    }
    N <- length(cate.score)

    if (N <= g) {
      stop("The number of groups must be smaller than number of predictions")
    }

    gamma_Delta_hat_i <- rep(NA, N)
    for (k in 1:N){
      gamma_hat_i <- gamma_Delta_hat[as.numeric(groups)[k]]
      gamma_Delta_hat_i[k] <- (N/g)/(N/g-1)*gamma_hat_i - 1/(N/g-1)*cate.score[k]
    }
    error <- mean((cate.score - cate.est) * (gamma_Delta_hat_i - cate.est))
  }
  return(error)
}
