setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")

# low-dimensional settings 
parameter <- read.csv("parameter.csv")
parameter <- parameter[1:12,]
outtC <- NULL
for (i in 1:12){
   out <- data.frame(read.csv(paste0("Results/updated_results_012722/",i,"RCT_AIPW.csv")))[,-1]
   bias <- colMeans(out[, 1:2])
   sds <- apply(out[,3:4], MARGIN = 2, FUN = "sd")
   MSEs <- bias^2 + sds^2
   std.bias <- bias / sds
   outtC <- rbind(outtC, cbind(rep(unlist(parameter[i,1]),2),rep(unlist(parameter[i,2]),2),rep(unlist(parameter[i,3]),2),
                               bias, sds, std.bias, MSEs))
}
outtC <- round(rbind(outtC[rownames(outtC)=="bias_plugin",], outtC[rownames(outtC)=="bias_robust",]),4)
write.csv(outtC, "RCT_AIPW.csv", row.names = F)

# high-dimensional settings 
parameter <- read.csv("parameter_ndim.csv")
outtC <- NULL
for (i in 1:30){
  out <- data.frame(read.csv(paste0("Results/updated_results_012722/",i,"RCT_AIPW_HD.csv")))
  out <- out[, !colnames(out) %in% c("X", "caseid")]
  names(out) <- c("bias_plugin",	"bias_robust",	"theta_plugin",	"theta_robust")
  bias <- colMeans(out[, 1:2])
  sds <- apply(out[,3:4], MARGIN = 2, FUN = "sd")
  MSEs <- bias^2 + sds^2
  std.bias <- bias / sds
  outtC <- rbind(outtC, cbind(rep(unlist(parameter[i,1]),2),rep(unlist(parameter[i,2]),2),rep(unlist(parameter[i,3]),2),
                              bias, sds, std.bias, MSEs))
}
outtC <- round(cbind(outtC[rownames(outtC)=="bias_robust",], outtC[rownames(outtC)=="bias_plugin",4:7]),4)
outtC <- data.frame(outtC)
names(outtC)[1:3] <- c("N","alpha","P")
outtC2 <- outtC[order(outtC$alpha, outtC$N), ]
write.csv(outtC2, "RCT_AIPW_HD.csv", row.names = F)







outtCsub <- data.frame(outtC[rownames(outtC)=="bias_robust",])
names(outtCsub)[1:3] <- c("N","alpha","P")
outtCformat <- outtCsub[order(outtCsub$N),c(1,3,4:7)]
write.csv(outtCformat, "obs_psC_HD_main.csv", row.names = F)



outtW <- matrix(NA, 16, 6)
for (i in 1:16){
  out <- data.frame(read.csv(paste0(i,"setting3_pswrong_output.csv")))[,3:4]
  means <- round(colMeans(out),4)
  sds   <- round(apply(out,2,sd),4)
  outtW[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtW, "outt_aipw_W.csv", row.names = F)


# AIPW high-dim covariates (p=10)
outtC <- matrix(NA, 16, 6)
for (i in 1:16){
   out <- data.frame(read.csv(paste0(i,"S3_psC_X0_X1_split_10HD.csv")))[,3:4]
   means <- round(colMeans(out),4)
   sds   <- round(apply(out,2,sd),4)
   outtC[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtC, "../S3D_psC_X0_X1_split_10HD.csv", row.names = F)

# AIPW high-dim covariates with p = 10
outtC <- matrix(NA, 16, 7)
for (i in 1:16){
   out <- data.frame(read.csv(paste0(i,"S3_psC_X0_X1_split_10HD.csv")))[,3:4]
   means <- round(colMeans(out),4)
   sds   <- round(apply(out,2,sd),4)
   outtC[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtC, "../S3D_psC_X0_X1_split_10HD.csv", row.names = F)

# AIPW high-dim covariates (p=50-400)
outtC <- matrix(NA, 40, 7)
for (i in 1:40){
   out <- data.frame(read.csv(paste0(i,"S3_psC_X0_X1_split_HD.csv")))[,3:4]
   means <- round(colMeans(out),4)
   sds   <- round(apply(out,2,sd),4)
   outtC[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtC, "../S3D_psC_X0_X1_split_HD.csv", row.names = F)


# AIPW high-dim covariates with large p = 20%N
outtC <- matrix(NA, 16, 7)
for (i in 1:16){
   out <- data.frame(read.csv(paste0(i,"S3_psC_X0_X1_split_HD_largeP.csv")))[,3:4]
   means <- round(colMeans(out),4)
   sds   <- round(apply(out,2,sd),4)
   outtC[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtC, "../S3D_psC_X0_X1_split_HD_largeP.csv", row.names = F)

out <- data.frame(read.csv("26S3_psC_split_HD_BootSE.csv"))


# Single cases 
outtC <- NULL
for (i in c(53:56, 69:72)){
   out <- data.frame(read.csv(paste0("single_results/", i, "S3_psC_forest_HD.csv")))[,4]
   means <- round(mean(out),4)
   sds   <- round(sd(out),4)
   outtC <- rbind(outtC, c(unlist(parameter[i,]), means, sds))
}
write.csv(outtC, "psC_forest_HD.csv")







