library(deSolve)
library(ggplot2)

rm(list = ls())
#setwd("~/HyperparasiteEvolution/Numerics")

SIHM<-function(t,y,p){ #Function specifying model, fed to deSolve
  S <- y["S"];
  I <- y["I"];
  H <- y["H"];
  M <- y["M"];
  with(as.list(p),{
    dS.dt <- b*S*(1-q*(S+I+H+M)) - d*S - betaI*S*I - c1*betaI*S*H - c1M*betaI*S*M + gammaI*I + (c3*gammaI + alphaH)*H + (c3M*gammaI + alphaHM)*M;
    dI.dt <- betaI*S*I + (1-p)*c1*betaI*S*H + (1-pM)*c1M*betaI*S*M + gammaH*H + gammaHM*M - (alphaI+gammaI+d)*I - betaH*I*H - betaHM*I*M;
    dH.dt <- betaH*I*H + (p)*c1*betaI*S*H - (c3*gammaI + alphaH + gammaH + c2*alphaI + d)*H;
    dM.dt <- betaHM*I*M + (pM)*c1M*betaI*S*M - (c3M*gammaI + alphaHM + gammaHM + c2M*alphaI + d)*M;
    return(list(c(dS.dt,dI.dt,dH.dt,dM.dt)));
  })
}
#Specify a lot of parameters
b <- 1
d <- 0.1
q <- 0.0005
alphaI <- 0.05
betaI <- 0.01
gammaI <- 1
alphaH <- 0.05
betaH <- 0.05
gammaH <- 0.01
c1 <- 1
c2 <- 1
c3 <- 1
p <- 0.2
alphaHM <- 0.05
betaHM <- 0.05
gammaHM <- 0.01
c1M <- 1
c2M <- 1
c3M <- 1
pM <- 0.2
hypt <- 200
mutt <- 50

a <- seq(0,1,0.05)
aM <- seq(0,1,0.05)
d <- c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)
crosses <- expand.grid(a,aM,d) #create grid of parameters to run model for

Pulse3 <- data.frame(var = c("M"),
                     time = c(mutt),
                     value = c(1),
                     method = c("add"))

tot_time <- 100
dt <- 0.1

t2 <- seq(from=0,to=tot_time,by=dt);

#Set ICs
S0 <- 900
I0 <- 100
H0 <- 1
M0 <- 0

N0 <- c(S=S0,I=I0,H=H0,M=M0)


for (i in 1:nrow(crosses)){
  if (i%%100 == 0 || i==1){
    print(paste("Run number",i))
  }
  ptemp <- list(b=b,d=d,q=q,betaI=betaI,betaH=betaH,
          gammaI=gammaI,gammaH=gammaH,alphaI=alphaI,alphaH=alphaH,hypt=hypt,
          c1=c1,c2=c2,c3=c3,p=p,
          c1M=c1M,c2M=c2M,c3M=c3M,pM=pM,
          betaHM=betaHM,gammaHM=gammaHM,alphaHM=alphaHM,mutt=mutt);
  
  ptemp[["alphaH"]] <- crosses[i,1] #set parameters for run according to tradeoffs
  ptemp[["betaH"]] <- (crosses[i,1])^(2/3)
  ptemp[["c1"]] <- 1 - crosses[i,1]
  ptemp[["c2"]] <- 1 - crosses[i,1]
  ptemp[["c3"]] <- crosses[i,1]
  ptemp[["alphaHM"]] <- crosses[i,2]
  ptemp[["betaHM"]] <- (crosses[i,2])^(2/3)
  ptemp[["c1M"]] <- 1 - crosses[i,2]
  ptemp[["c2M"]] <- 1 - crosses[i,2]
  ptemp[["c3M"]] <- crosses[i,2]
  ptemp[["d"]] <- crosses[i,3]

  

  odeSIHM <- data.frame(ode(y = N0,
                            times = t2,
                            func=SIHM, parms=ptemp,
                            method = "ode45", events = list(data = Pulse3)))

  colnames(odeSIHM) <- c("time", "S", "I","H","M")


  if (1 < odeSIHM[nrow(odeSIHM),"M"]){ #If M is resident strain, invasion occured
    crosses[i,4] <- "Invade"
  } else { #otherwise, no invasion occured
    crosses[i,4] <- "No Invade"
  }
  
}
colnames(crosses) <- c("a","aM","d","Invade")

crosses_fulltrade_overd2 <- crosses
save(crosses_fulltrade_overd2,
     file="NumericalInvasionsFullTradeoffRangeD.Rdata")


