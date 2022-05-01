#' @title AIPW score
#'
#' @description  Compute AIPW scores as surrogate individual treatment effects
#'
#' @param X Covariate matrix
#' @param W A vector of randomized treatment indicator
#' @param Y Outcome variable (observed follow-up time for survival outcomes)
#' @param D Event indicator
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
#'
#' # Compute the AIPW score
#' cate.score <- aipw_score(X = X, W = W, Y = Y)
#'
#' # time-to-event outcome
#' horizon <- 1
#' failure.time <- pmin(rexp(n) * X[, 1] + W, horizon)
#' censor.time <- 2 * runif(n)
#' # Discretizing continuous events decreases runtime.
#' Y <- round(pmin(failure.time, censor.time), 2)
#' D <- as.integer(failure.time <= censor.time)
#' t0 <- median(Y)
#'
#' # Compute the AIPW score for survival outcomes
#' cate.surv.score <- aipw_surv_score(X = X, W = W, Y = Y, D = D, t0 = t0)
#' }
#' @return
#' @export
#'
#' # AIPW-based score for RCT, non-survival outcome
aipw_score <- function(X, W, Y){
  fit <- grf::regression_forest(X = cbind(X, W), Y = Y)   # outcome estimates
  u  <- predict(fit)$predictions
  u1 <- predict(fit, cbind(1, X))$predictions
  u0 <- predict(fit, cbind(0, X))$predictions
  ps <- mean(W)
  score <- unlist(u1 - u0 + (W - ps) / (ps * (1 - ps)) * (Y - u))
  return(score)
}

# AIPW-based score for RCT, survival outcome
aipw_surv_score <- function(X, W, Y, D, t0){
  # survival and censoring curves
  S1 <- survival_forest(X = X[W==1, ], Y = Y[W==1], D = D[W==1], prediction.type = "Nelson-Aalen", failure.times = sort(unique(Y)))
  surv1 <- predict(S1, X, failure.times = sort(unique(Y)))$predictions

  S0 <- survival_forest(X = X[W==0, ], Y = Y[W==0], D = D[W==0], prediction.type = "Nelson-Aalen", failure.times = sort(unique(Y)))
  surv0 <- predict(S0, X, failure.times = sort(unique(Y)))$predictions

  C1 <- survival_forest(X = X[W==1,], Y = Y[W==1], D = 1 - D[W==1], prediction.type = "Nelson-Aalen", failure.times = sort(unique(Y)))
  cen1 <- predict(C1, X, failure.times = sort(unique(Y)))$predictions

  C0 <- survival_forest(X = X[W==0,], Y = Y[W==0], D = 1 - D[W==0], prediction.type = "Nelson-Aalen", failure.times = sort(unique(Y)))
  cen0 <- predict(C0, X, failure.times = sort(unique(Y)))$predictions

  # AIPW score
  score <- compute_scores(treatment = W,
                          study_time = Y,
                          status = D,
                          propensity_score = rep(mean(W), nrow(X)),
                          treated_survival_curve = surv1,
                          treated_censoring_survival_curve = cen1,
                          control_survival_curve = surv0,
                          control_censoring_survival_curve = cen0,
                          unique_times = sort(unique(Y)),
                          end_time = t0,
                          rmst = FALSE)
  return(score)
}
