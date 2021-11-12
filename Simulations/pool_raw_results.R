setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE/raw results")
parameter <- read.csv("../parameter_ndim.csv")
outtC <- matrix(NA, 16, 6)
for (i in 1:16){
   out <- data.frame(read.csv(paste0(i,"S3_psC_X0_X1_split.csv")))[,3:4]
   means <- round(colMeans(out),4)
   sds   <- round(apply(out,2,sd),4)
   outtC[i,] <- c(unlist(parameter[i,]), means, sds)
}
write.csv(outtC, "../S3C_psC_X0_X1_split.csv", row.names = F)

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







