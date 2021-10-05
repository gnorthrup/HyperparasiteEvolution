#setwd("~/HyperparasiteEvolution/Numerics")

library(ggplot2)

getESS <- function(crosses,parname,parval){ #Function to collect ESS from a grid of parameter values
  sub <- subset(crosses, crosses[[parname]]==parval) #for wild-type and mutant hyperparasites
  alphas <- unique(sub$a)
  scores <- c()
  for(i in 1:length(alphas)){
    score <- sum(sub$a == alphas[i] & sub$Invade == "No Invade") + sum(sub$aM == alphas[i] & sub$Invade == "Invade")
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

load("NumericalInvasionsFullTradeoffRangeD.Rdata")

d_ESSs <- getESSvec(crosses_fulltrade_overd2,"d")

trimmed <- d_ESSs[!(duplicated(d_ESSs$ESSs)),]

p1 <- ggplot(d_ESSs, aes(parvals,ESSs)) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Host Lifespan", y = "Hyperparasite\nvirulence") +
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=15)) +
  scale_x_reverse()

print(p1)
