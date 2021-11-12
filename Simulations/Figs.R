setwd("C:/Users/cryst/Documents/Stanford Postdoc/NHLBI R01 Aim 2/Calibration_HTE")
library(ggplot2)
library(gridExtra)
ipw_aipw_fig <- read.csv("./ipw_aipw_fig.csv")
fig1 <- ipw_aipw_fig[ipw_aipw_fig$study.type=="RCT",]
fig1$PS.model <- c(rep("None", 4), rep("CY", 4),rep("None", 4), rep("CY", 4))
fig1$S.bias <- abs(fig1$S.bias)
fig1long <- data.frame(rbind(as.matrix(fig1[, c(1,3:6)]), as.matrix(fig1[, c(1,3:5,7)])))
fig1long$Ytype <- c(rep("Std.bias",16), rep("MSE",16))
fig1long$group <- ifelse(fig1long$estimator=="IPW" & fig1long$Metric=="plug-in", 1, 
                         ifelse(fig1long$estimator=="IPW" & fig1long$Metric=="robust", 2, 
                                ifelse(fig1long$estimator=="AIPW" & fig1long$Metric=="plug-in", 3, 4)))
names(fig1long)[2] <- "Estimator"

p1 <- ggplot(fig1long, aes(x=N, y=S.bias, col=Estimator, linetype=Metric, group = group)) + 
  geom_line(lwd=1.2) + 
  geom_point()+
  theme_classic()+
  xlab("Sample Size")+
  theme(axis.title.y=element_blank())+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  facet_wrap(~Ytype)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=8)); p1


fig2 <- ipw_aipw_fig[ipw_aipw_fig$study.type=="Obs",]
fig2$PS.model <- c(rep("C_PS",4), rep("W_PS",4), rep("C_PS",4), rep("W_PS",4),
                   rep("C_PS",4), rep("W_PS",4), rep("C_PS",4), rep("W_PS",4))
fig2$S.bias <- abs(fig2$S.bias)
fig2long <- data.frame(rbind(as.matrix(fig2[, c(1,3:6)]), as.matrix(fig2[, c(1,3:5,7)])))
fig2long$Ytype <- c(rep("Std.bias",32), rep("MSE",32))
fig2long$group <- ifelse(fig2long$estimator=="IPW" & fig2long$PS.model=="C_PS" & fig2long$Metric=="plug-in", 1, 
                     ifelse(fig2long$estimator=="IPW" & fig2long$PS.model=="C_PS" & fig2long$Metric=="robust", 2,
                            ifelse(fig2long$estimator=="IPW" & fig2long$PS.model=="W_PS" & fig2long$Metric=="plug-in", 3, 
                                   ifelse(fig2long$estimator=="IPW" & fig2long$PS.model=="W_PS" & fig2long$Metric=="robust", 4,
                                          ifelse(fig2long$estimator=="AIPW" & fig2long$PS.model=="C_PS" & fig2long$Metric=="plug-in", 5, 
                                                 ifelse(fig2long$estimator=="AIPW" & fig2long$PS.model=="C_PS" & fig2long$Metric=="robust", 6,
                                                        ifelse(fig2long$estimator=="AIPW" & fig2long$PS.model=="W_PS" & fig2long$Metric=="plug-in", 7, 
                                                               ifelse(fig2long$estimator=="AIPW" & fig2long$PS.model=="W_PS" & fig2long$Metric=="robust", 8,NA))))))))
names(fig2long)[2:3] <- c("Estimator", "Model")
p2 <- ggplot(fig2long, aes(x=N, y=S.bias, col=interaction(Estimator,Model), linetype=Metric, group = group)) + 
  geom_line(lwd=1.2) + 
  geom_point()+
  theme_classic()+
  xlab("Sample Size")+
  theme(axis.title.y=element_blank())+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
  guides(colour = guide_legend(nrow = 4))+
  facet_wrap(~Ytype)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(size=8));p2

png("./ipw_aipw_comp.png", width = 5, height = 10, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(p1, p2, nrow=2, heights=c(2,3)), nrow=2, heights=c(100,1)))
dev.off()

