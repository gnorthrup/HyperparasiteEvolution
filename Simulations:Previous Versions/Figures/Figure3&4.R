#setwd("~/HyperparasiteEvolution/Numerics")

library(ggplot2)

getESS <- function(crosses,parname,parval){ #Function to collect ESS from a grid of parameter values
  sub <- subset(crosses, crosses[[parname]]==parval) #for wild-type and mutant hyperparasites
  alphas <- unique(sub$alphaH)
  scores <- c()
  for(i in 1:length(alphas)){
    score <- sum(sub$alphaH == alphas[i] & sub$Invade == "No Invade") + sum(sub$alphaHM == alphas[i] & sub$Invade == "Invade")
    scores[i] <- score
  }
  return(alphas[which.max(scores)])
}

getESSvec <- function(crosses,parname){
  parvals <- unique(crosses[[parname]])
  ESSs <- c()
  for(i in 1:length(parvals)){
    ESSs[i] <- getESS(crosses,parname,parvals[i])
  }
  out <- data.frame(parvals,ESSs)
}

load("NumericalInvasionsFullTradeoffRangeAlphaI.Rdata")

alphaI_ESSs <- getESSvec(crosses_fulltrade_overalphaI,"alphaI")

load("NumericalInvasionsFullTradeoffRangeP.Rdata")

p_ESSs <- getESSvec(crosses_fulltrade_overp,"p")

load("NumericalInvasionsFullTradeoffRangeGammaI.Rdata")

gammaI_ESSs <- getESSvec(crosses_fulltrade_largegammaIsweep,"gammaI")

trimmed <- alphaI_ESSs[!(duplicated(alphaI_ESSs$ESSs)),]

p1 <- ggplot(trimmed, aes(parvals,ESSs)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Parasite Virulence", y = "Hyperparasite\nvirulence") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15))

print(p1)

p2 <- ggplot(p_ESSs, aes(parvals,ESSs)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Hitchhiking Probability", y = "") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15))

print(p2)

trimmed2 <- gammaI_ESSs[!(duplicated(gammaI_ESSs$ESSs)),]

p3 <- ggplot(trimmed2, aes(parvals,ESSs)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Host Recovery Rate", y = "") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15))

print(p3)

library(ggpubr)

ggarrange(p1,p2,p3,nrow=1,labels = c("A)","B)","C)"),widths=c(1.15,1,1))

###################
#####FIGURE 4######
###################

p4 <- ggplot(trimmed, aes(parvals,1-ESSs)) + #Invert virulence to get direction of selection on c_alpha
  geom_point(size=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Parasite Virulence", y = expression(paste("Parasite\nVirulence\nModifier",(c[alpha]),sep=""))) +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15),plot.margin = margin(10, 10, 10, 20))
  

print(p4)

p5 <- ggplot(p_ESSs, aes(parvals,1-ESSs)) +
  geom_point(size=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Hitchhiking Probability", y = "") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15),plot.margin = margin(10, 10, 10, 20))

print(p5)

trimmed2 <- gammaI_ESSs[!(duplicated(gammaI_ESSs$ESSs)),]

p6 <- ggplot(trimmed2, aes(parvals,1-ESSs)) +
  geom_point(size=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Host Recovery Rate", y="") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15))

print(p6)

library(ggpubr)

ggarrange(p4,p5,p6,nrow=1,labels = c("A)","B)","C)"),widths=c(1.4,1,1))
