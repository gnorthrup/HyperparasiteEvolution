#setwd("~/HyperparasiteEvolution/Numerics")

library(ggplot2)

load("NumericalInvasionsBasicTradeoff.Rdata") #Numerical soutions for invasion with basic virulence transmission tradeoff
                                          #Over a grid of parameter values for p and c_beta

sub <- subset(crosses_ABtrade_overpandc1, c1==1 & p %in% c(0.8,0.9,1))

formatp <- function(strings) {
  return(paste("p =",strings))
}

p.labs <- c("Low hitchhiking","Intermediate hitchiking","Complete hitchiking")
names(p.labs) <- c("0.8", "0.9", "1")

figure2 <- ggplot(sub, aes(alphaH,alphaHM)) +
  geom_raster(aes(fill = as.factor(Invade))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none', axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x = "Wild-type hyperparasite virulence", y = "Mutant\nhyperparasite\nvirulence") +
  #facet_wrap(~p,labeller=labeller(p=formatp))
  facet_wrap(~p,labeller=labeller(p=p.labs))+
  theme(axis.title.y = element_text(angle = 45,vjust = 0.5,size=20))
  
print(figure2)
