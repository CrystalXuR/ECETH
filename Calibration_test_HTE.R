### Hosmer-Lemeshow test for HTE calibration (may need to be revised) 
# yhat -- predicted event risk
# delta -- predicted risk difference 
# Y -- observed follow up time
# W -- treatment variable 
# time -- time of interest 
# g -- number of group

hosmerlem <- function(yhat, delta, Y, W, time, g=10) {
  
  cut.delta <- cut(delta, breaks = unique(quantile(delta, probs=seq(0,1,1/g))), include.lowest=TRUE)
  
  sfit.t0.delta <- summary(survfit(Surv(Y[W==0], D[W==0]==1) ~ cut.delta[W==0], se.fit = TRUE), times = time)
  sfit.t1.delta <- summary(survfit(Surv(Y[W==1], D[W==1]==1) ~ cut.delta[W==1], se.fit = TRUE), times = time)
  
  obs.risk.t0.tmp <- 1- sfit.t0.delta$surv
  obs.risk.t1.tmp <- 1- sfit.t1.delta$surv
  
  obs.sd.t0.tmp <- sfit.t0.delta$std.err
  obs.sd.t1.tmp <- sfit.t1.delta$std.err
  
  # observed ARD
  obs.risk.delta <- obs.risk.t1.tmp - obs.risk.t0.tmp
  obs.var.delta <- obs.sd.t0.tmp^2 + obs.sd.t1.tmp^2
  
  # expected ARD
  exp.delta   <- aggregate(delta, by = list(cut.delta), FUN = "mean")$x
  #ng.delta    <- as.numeric(unlist(table(cut.delta)))
  
  # GND for HTE statistic
  hl.delta <- sum((obs.risk.delta - exp.delta)^2 /obs.var.delta)
  pval.delta <- 1 - pchisq(hl.delta, g-1)
  
  return(list(chisq=hl.delta, p.value=pval.delta))
}