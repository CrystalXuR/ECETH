library(matrixStats)

# Return the weighted average and standard deviation.
# values, weights -- Numpy nparrays with the same shape.
weighted_avg_and_std <- function(values, weights){
  out <- list(NA,2)
  out$average  <- weighted.mean(values, weights)
  out$std <- sqrt(weighted.mean((values-out$average)**2, weights))
  return(out)
}

# Potentially useful
# Takes in a kaplan meier object and a set of time points to integrate against
# Computes the RMST of that object up to the end time
# Mean survival time up to some point (i.e., min(survival time, end time))
rmst_km <- function(km_obj, unique_times, end_time){
  S = cbind(summary(km_obj,times=unique_times)$time, 
            summary(km_obj,times=unique_times)$surv)
  return(S)
}

rmst <- function(S, unique_times, end_time){
  # :param S: A [num_observations x num_unique_times] array representing the survival at time T for each unique time
  # :param unique_times:  An iterable containing the list of unique times
  # :param end_time:
  #   :return:
  
  col_indx <- length(which(unique_times<=end_time))
  S <- S[,1:col_indx]
  unique_times <- unique_times[unique_times <= end_time]
  time_diffs <- rep(NA,(length(unique_times)-1))
  for (i in 1:(length(unique_times)-1)){
    time_diffs[i] <- unique_times[i+1]-unique_times[i]
  }
  time_diffs <- c(unique_times[1], time_diffs) 
  rmst <- S%*%time_diffs
  return(rmst)
}

# Get end time
get_good_end_time <-function(study_time, status, treatment, censoring_survival_curve, unique_times){
  W <- treatment 
  dat <- data.frame(study_time,status,W)
  dat1 <- dat[dat$W==1,]
  kmt.fit <- survfit(Surv(study_time,1-status) ~ 1, data=dat1)
  S_time_C_t <- cbind(summary(kmt.fit,times=unique_times-1e-7)$time, 
                      summary(kmt.fit,times=unique_times-1e-7)$surv)
  censoring_survival_curve[W==1,] <- t(matrix(rep(S_time_C_t[,2], dim(dat1)[1]),length(S_time_C_t[,2]), dim(dat1)[1]))            
  
  dat0 <- dat[dat$W==0,]
  kmc.fit <- survfit(Surv(study_time,1-status) ~ 1, data=dat0)
  S_time_C_c <- summary(kmc.fit, times= unique_times-1e-7)$surv
  if (length(S_time_C_c) < length(unique_times)){
    S_time_C_c <- c(S_time_C_c, rep(S_time_C_c[length(S_time_C_c)], (length(unique_times) - length(S_time_C_c))))
  }
  censoring_survival_curve[W==0,] <- t(matrix(rep(S_time_C_c, dim(dat0)[1]),length(S_time_C_c), dim(dat0)[1]))
  
  LTt <- data.frame(summary(kmt.fit)$time, summary(kmt.fit)$surv); names(LTt) <- c("time","surv")
  index_t <- which(abs(LTt$surv-0.3)==min(abs(LTt$surv-0.3)))+1
  
  LTc <- data.frame(summary(kmc.fit)$time, summary(kmc.fit)$surv); names(LTc) <- c("time","surv")
  index_c <- which(abs(LTc$surv-0.3)==min(abs(LTc$surv-0.3)))+1
  
  # The .percentile(p) function in Python finds the corresponding time for a given survival prob. 
  # There is no such a function in R, so an approx. is used here
  et <- min(5*365, LTt$time[index_t], LTc$time[index_c])
  if(any(summary(kmt.fit,times=et)$surv<0.2)==TRUE|any(summary(kmc.fit,times=et)$surv<0.2)==TRUE){
    et <- et-1
  }
  out <- list(NA, 2)
  out[[1]] <- censoring_survival_curve
  out[[2]] <- et
  return(out) 
}
  
cond_mean <- function(S, t, end_time){
  valid <- (t<=end_time)
  tempdat <- rbind(valid, as.matrix(S))
  S <- tempdat[-1,tempdat[1,]==1]
  tempdat <- data.frame(t,valid)
  t <- tempdat[tempdat$valid==1,]$t
  time_diffs <- rep(NA,(length(t)-1))
  for (i in 1:(length(t)-1)){
    time_diffs[i] <- t[i+1]-t[i]
  }
  time_diffs <- c(t[1], time_diffs)
  Sd <- t(time_diffs*t(S))
  
  Sdcum <- rowCumsums(Sd[,order(ncol(Sd):1)],)
  cond_mean <- Sdcum[,order(ncol(Sdcum):1)]   #approx. RMST
  cond_mean <- t(t + t(cond_mean))
  B <- rep(end_time, dim(S)[1])
  return(data.frame(cond_mean,B))
}
  
cond_surv <- function(S, t, end_time){
  valid <- (t<=end_time)
  tempdat <- rbind(valid, as.matrix(S))
  S <- tempdat[-1,tempdat[1,]==1]
  cond_surv <- matrix(1,dim(S)[1],dim(S)[2])
  
  # TODO: Verify that this change is what we want.
  S.std <- S[,dim(S)[2]]/S
  cond_surv[(S > 0)] <- S.std[S>0]
  return(cond_surv)
}

compute_scores_aug <- function (treatment, study_time, status, propensity_score, survival_curve, 
                                censoring_survival_curve, unique_times, end_time, rmst=TRUE, 
                                ipw=TRUE, time_lost=FALSE){
  n <- length(treatment)
  W <- treatment
  e <- propensity_score
  tempdat <- rbind(unique_times,survival_curve)
  S = survival_curve[,tempdat[1,]<=end_time]
  tempdat <- rbind(unique_times,survival_curve)
  S_c = censoring_survival_curve[,tempdat[1,]<=end_time]
  
  if(rmst==TRUE){
    Q <- cond_mean(survival_curve, unique_times, end_time)
    if(time_lost==TRUE){
      Q <- end_time - Q
    }
  }else{
    Q <- cond_surv(survival_curve, unique_times, end_time)
  }
  
  unique_times <- unique_times[unique_times <= end_time]
  nt = length(unique_times)
  
  S_c_z <- S_c
  dH_c <- cbind(1-S_c[,1],-(S_c[,-1]-S_c[,1:(dim(S_c)[2]-1)]))
  dH_c <- dH_c/S_c_z
  
  if(rmst==TRUE){
    Y <- pmin(study_time, end_time)
    if(time_lost==TRUE){
      Y <- end_time-Y
    }
  }else{
    Y <- as.numeric(study_time>=end_time)
  }
  mu <- Q[,1]
  
  at_risk <- matrix(NA,length(Y),length(unique_times))
  for (i in 1:length(study_time)){
    at_risk[i,] <- as.numeric(study_time[i]>unique_times)
  }
  
  if(dim(Q)[2] > dim(at_risk)[2]){
    Q = Q[, 1:(dim(at_risk)[2])]
  }
  compensator <- apply((at_risk * dH_c * Q / S_c_z),1,sum)
  
  study_end_time <- pmin(study_time, end_time)
  U_index <- findInterval(study_end_time, unique_times)
  Q_U_X <- rep(NA, n)
  for (i in 1:n){
    Q_U_X[i] <- Q[i, pmin(U_index, nt)[i]]
  }
 
  S_U_X <- rep(NA,n)
  for (i in (1:n)){
    S_U_X[i] <- S_c[i, pmin((U_index + 1), nt)[i]]
  }
  
  D <- ((status==TRUE | (study_time > end_time)) & (S_U_X > 0))
  S_U_X[S_U_X == 0] <- 1
  
  if (ipw==TRUE){
    dr_score <- (W / e) * (D * Y / S_U_X) / mean((W / e) * (D / S_U_X))
  }else{
    dr_score <- mu + (W / e) * ((D * Y + (1 - D) * Q_U_X) / S_U_X - compensator - mu) / mean((W / e) * (D / S_U_X))
  }
  return(dr_score)
}

compute_scores <- function(treatment, study_time, status,
                           propensity_score,
                           treated_survival_curve, treated_censoring_survival_curve,
                           control_survival_curve, control_censoring_survival_curve,
                           unique_times, end_time, censoring_filter = 0.2,
                           ipw = FALSE, rmst = TRUE, censoring_free=FALSE, time_lost=FALSE){
  # Args:
  #   treatment: np.array with 1 if treated and 0 otherwise
  #   study_time: minimum
  # 
  # Computes the AIPW score with doubly-robust censoring adjustment for outcome Y = 1{T >= end_time}
  # or if rmst = True, Y = min(T, end_time).
  # Assumes survival_curve and censoring_survival_curve are matrices specifying for each
  # study unit (row) and time (column) in the list of times in unique_times (which is assumed to be
  # in increasing order), the probability P(T >= t | X=x) or P(C >= t | X=x), where x are the
  # covariates. Treatment is assumed binary, with 1 if treated and 0 if control and propensity score
  # is P(treated | X=x). End time should be soon enough that there isn't a subgroup of patients where
  # P(C >= end_time | X=x) is close to 0. status is 1 if study_time is the time of the event, and 0 if
  # it's the time of censorship.
  
  if(censoring_free==TRUE){
    treated_censoring_filter_threshold <- quantile(treated_censoring_survival_curve[,dim(treated_censoring_survival_curve)[2]], censoring_filter)
    control_censoring_filter_threshold <- quantile(control_censoring_survival_curve[,dim(control_censoring_survival_curve)[2]], censoring_filter)
    
    above_threshold <- rep(NA, dim(treated_censoring_survival_curve)[1])
    for (i in 1:dim(treated_censoring_survival_curve)[1]){
      above_threshold[i] <- (all(treated_censoring_survival_curve[i,] >= treated_censoring_filter_threshold)==TRUE &
                             all(control_censoring_survival_curve[i,] >= control_censoring_filter_threshold)==TRUE)
    }
    
    if(mean(as.numeric(above_threshold))<0.3){
      print(paste0("Warning: dropping ",round((1-mean(as.numeric(above_threshold)))*100,1),"% of the data due to censoring issue"))
    }
    treatment <- treatment[above_threshold]
    study_time <- study_time[above_threshold]
    status <- status[above_threshold]
    propensity_score <- propensity_score[above_threshold]
    treated_survival_curve <- treated_survival_curve[above_threshold,]
    # treated_cond_mean <- treated_cond_mean[above_threshold,]
    treated_censoring_survival_curve <- treated_censoring_survival_curve[above_threshold, ]
    control_survival_curve <- control_survival_curve[above_threshold,]
    # control_cond_mean <- control_cond_mean[above_threshold,]
    control_censoring_survival_curve <- control_censoring_survival_curve[above_threshold,]
  }
  print("\t\tComputing mu1...")
  mu1_score <- compute_scores_aug(treatment, study_time, status, propensity_score, treated_survival_curve,
                                  treated_censoring_survival_curve, unique_times, end_time, rmst, ipw, time_lost)
  print("\t\tComputing mu0...")
  mu0_score <- compute_scores_aug(1-treatment, study_time, status, 1-propensity_score, control_survival_curve,
                                  control_censoring_survival_curve, unique_times, end_time, rmst, ipw, time_lost)
  score <- mu1_score - mu0_score
  
  return(score)
}

# test <- function(){
#   set.seed(1)
#   n = 10000
#   treatment  <- rbinom(n,1,0.5)    
#   study_time <- runif(n,0,1)
#   study_cens <- runif(n,0,1)
#   status     <- as.numeric(study_time <= study_cens)
#   fu_time <- pmin(study_time, study_cens)
#   propensity_score <- rep(0,n) + 0.5
#   unique_times     <- sort(fu_time)
#   survival_curve   <- t(replicate(n,1 - unique_times)) 
#   censoring_survival_curve <- t(replicate(n,1 - unique_times)) 
#   end_time <- 0.5
#   
#   print(mean(compute_scores(treatment, fu_time, status, propensity_score, survival_curve,
#                             censoring_survival_curve, survival_curve, censoring_survival_curve, unique_times, end_time)))
# }

# data = read.csv("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/risk-vs-hte-R/example_survival_aipw_score.csv")
# attach(data)
# mean(compute_scores(treatment, study_time, status, propensity_score, data[,5:104],
#                     data[,105:204], data[,5:104], data[,105:204], unique_times, 0.5,
#                     rmst=FALSE))
# python code result = 0.09727313182556326, checked!


  
  
  











