setwd("/Users/yizhe/Desktop/Crystal Xu/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")
library(ggplot2)
library(gridExtra)

# Load RCT results
ipw_RCT <- read.csv("RCT_IPW.csv"); ipw_RCT <- ipw_RCT[, -1]
ipw_RCT <- data.frame(rbind(as.matrix(ipw_RCT[, c(1:3, 5,7,9)]), as.matrix(ipw_RCT[, c(1:2, 4,6,8, 10)])))
ipw_RCT$Metric <- c(rep("plug-in", 12), rep("robust", 12))
aipw_RCT <- read.csv("RCT_AIPW_cf_all.csv")
aipw_RCT$Metric <- c(rep("plug-in", 12), rep("robust", 12))
names(aipw_RCT) <- names(ipw_RCT) <- c("N", "alpha", names(aipw_RCT)[3:7])
RCT_dat <- rbind(ipw_RCT, aipw_RCT)
RCT_dat$Estimator <- c(rep("IPW", 24), rep("AIPW", 24))
RCT_dat <- RCT_dat[RCT_dat$alpha==0.15,c(1,5:8)]
fig1long <- data.frame(rbind(as.matrix(RCT_dat[,c(1,2,4,5)]), as.matrix(RCT_dat[,c(1,3,4,5)])))
fig1long$Ytype <- c(rep("Std.bias",dim(RCT_dat)[1]), rep("MSE",dim(RCT_dat)[1]))
# fig1long$Estimator <- ifelse(fig1long$Estimator=="AIPW", "AIPW_C_PS", "IPW_C_PS")
fig1long$group <- ifelse(fig1long$Estimator=="AIPW" & fig1long$Metric=="plug-in", 3, 
                         ifelse(fig1long$Estimator=="AIPW" & fig1long$Metric=="robust", 7, 
                                ifelse(fig1long$Estimator=="IPW" & fig1long$Metric=="plug-in", 1, 5)))
#fig1long$std.bias <- log(abs(as.numeric(fig1long$std.bias)))

# Load observational results
ipw_psC <- read.csv("obs_psC_IPW.csv"); ipw_psC <- ipw_psC[, -1]
ipw_psC <- data.frame(rbind(as.matrix(ipw_psC[, c(1:3, 5,7,9)]), as.matrix(ipw_psC[, c(1:2, 4,6,8, 10)])))
ipw_psC$Metric <- c(rep("plug-in", 12), rep("robust", 12))

ipw_psW <- read.csv("obs_psW_IPW.csv"); ipw_psW <- ipw_psW[, -1]
ipw_psW <- data.frame(rbind(as.matrix(ipw_psW[, c(1:3, 5,7,9)]), as.matrix(ipw_psW[, c(1:2, 4,6,8, 10)])))
ipw_psW$Metric <- c(rep("plug-in", 12), rep("robust", 12))

aipw_psW <- read.csv("obs_psW_cf_all.csv")
aipw_psW$Metric <- c(rep("plug-in", 12), rep("robust", 12))
aipw_psC <- read.csv("obs_psC_cf_all.csv")
aipw_psC$Metric <- c(rep("plug-in", 12), rep("robust", 12))

names(aipw_psW) <- names(aipw_psC) <- names(ipw_psW) <- names(ipw_psC) <- c("N", "alpha", names(aipw_psC)[3:7])
obs_dat <- rbind(ipw_psC, ipw_psW, aipw_psC, aipw_psW)
obs_dat$Estimator <- c(rep("IPW", 24), rep("IPW misspecified", 24), rep("AIPW", 24), rep("AIPW misspecified", 24))
obs_dat <- obs_dat[obs_dat$alpha==0.15,c(1,5:8)]
fig2long <- data.frame(rbind(as.matrix(obs_dat[,c(1,2,4,5)]), as.matrix(obs_dat[,c(1,3,4,5)])))
fig2long$Ytype <- c(rep("Std.bias",dim(obs_dat)[1]), rep("MSE",dim(obs_dat)[1]))
fig2long$group <- ifelse(fig2long$Estimator=="IPW" & fig2long$Metric=="plug-in", 1,
                         ifelse(fig2long$Estimator=="IPW misspecified" & fig2long$Metric=="plug-in", 2,
                                ifelse(fig2long$Estimator=="AIPW" & fig2long$Metric=="plug-in", 3,
                                       ifelse(fig2long$Estimator=="AIPW misspecified" & fig2long$Metric=="plug-in", 4,
                                              ifelse(fig2long$Estimator=="IPW" & fig2long$Metric=="robust", 5,
                                                     ifelse(fig2long$Estimator=="IPW misspecified" & fig2long$Metric=="robust", 6,
                                                            ifelse(fig2long$Estimator=="AIPW" & fig2long$Metric=="robust", 7, 8)))))))
#fig2long$std.bias <- log(abs(as.numeric(fig2long$std.bias)))

# Reorganize figures by metric type 
MSEdat <- cbind(rbind(fig1long[fig1long$Ytype=="MSE", ], fig2long[fig2long$Ytype=="MSE", ]),c(rep("RCT",16), rep("Obs",32)))
names(MSEdat)[c(2,dim(MSEdat)[2])] <- c("Y", "Study_Type")
MSEdat$Y <- as.numeric(MSEdat$Y)
p1 <- ggplot(MSEdat, aes(x=N, y=Y, col=Estimator, linetype=Metric, group = group)) + 
  geom_line(lwd=1.2) + 
  scale_color_manual(values = rep(c("red1", "goldenrod", "gray0", "gray50"),2)) +
  geom_point()+
  theme_classic()+
  xlab("Sample Size")+
  theme(axis.title.y = element_blank())+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  facet_wrap(~Study_Type); p1

stddat <- cbind(rbind(fig1long[fig1long$Ytype=="Std.bias", ], fig2long[fig2long$Ytype=="Std.bias", ]),c(rep("RCT",16), rep("Obs",32)))
names(stddat)[c(2,dim(stddat)[2])] <- c("Y", "Study_Type")
stddat$Y <- as.numeric(stddat$Y)
p2 <- ggplot(stddat, aes(x=N, y=Y, col=Estimator, linetype=Metric, group = group)) + 
  geom_line(lwd=1.2) + 
  scale_color_manual(values = rep(c("red1", "goldenrod", "gray0", "gray50"),2)) +
  geom_point()+
  theme_classic()+
  xlab("Sample Size")+
  theme(axis.title.y = element_blank())+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  facet_wrap(~Study_Type); p2


library(ggpubr)
png("ipw_aipw_comp.png", width = 5, height = 8, units = 'in', res = 300)
print(ggarrange(p1, p2, nrow=2, common.legend = TRUE, legend="bottom"))
dev.off()

