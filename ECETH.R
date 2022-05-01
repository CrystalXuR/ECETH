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
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' event.prob <- 1 / (1 + exp(2*(pmax(2*X[, 1], 0) * W - X[, 2])))
#' Y <- rbinom(n, 1, event.prob)
#' train <- sample(1:n, n / 2)
#' cf.fit <- causal_forest(X[train, ], Y[train], W[train])
#' cf.cate <- -1 * predict(cf.fit, X[-train, ])$predictions
#'
#' # Compute the AIPW score
#' cate.score <- aipw.score(X = X[-train, ], W = W[-train], Y = Y[-train])
#'
#' # ECETH
#' caliberror <- eceth(cf.cate, cate.score)
#' }
#' @return
#' @export

# AIPW-based score for RCT, non-survival outcome
aipw.score <-  function(X, W, Y){
  fit <- grf::regression_forest(X = cbind(X, W), Y = Y)   # outcome estimates
  u  <- predict(fit)$predictions
  u1 <- predict(fit, cbind(1, X))$predictions
  u0 <- predict(fit, cbind(0, X))$predictions
  ps <- mean(W)
  score <- unlist(u1 - u0 + (W - ps) / (ps * (1 - ps)) * (Y - u))
  return(score)
}

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
