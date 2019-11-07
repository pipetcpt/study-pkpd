setwd("D:/Rt/Gab/")
require(wnl)
dPK32 = read.csv("PK32.csv", skip=1)
colnames(dPK32) = c("TIME", "DV") ; dPK32

require(deSolve)

InfHist = cbind(TIME=c(0,           30.1, 125.2,                154.3, 260.1,   290.1), 
                RATE=c(1131.8/30.1, 0,    1884.4/(154.3-125.2), 0,     6300/30, 0)) ; InfHist

nBolus = 3
EventDat = data.frame(var = rep("y1", nBolus),
                      time = c(0, 125, 260),
                      value = c(1669/5.77348, 1701/5.77348, 1733/5.77348),
                      method = rep("add", nBolus)) ; EventDat

PKde1 = function(t, y, p)
{
  RateIn = InfHist[findInterval(t, InfHist[,"TIME"]),"RATE"]
  
  dy1dt = (RateIn + p["Kin"] - p["CL"]*y[1])/p["V"]
  return(list(c(dy1dt)))
}

PKde2 = function(t, y, p)
{
  RateIn = InfHist[findInterval(t, InfHist[,"TIME"]),"RATE"]
  
  dy1dt = (RateIn + p["Kin"] - p["Vmax"]*y[1]/(p["Km"] + y[1]))/p["V"]
  return(list(c(dy1dt)))
}

## calculate steady-state concentration for nonlinear elimination model
Vmax = 361.502
Km = 507.873
Kin = 14.9684
Cbase = Kin/Vmax*Km
for (i in 1:1000) {
  Cbase = Kin/Vmax*(Km + Cbase)
} ; Cbase # 21.93744
###

TIME = sort(c(0, 125, 260, dPK32[,"TIME"]))

y1 = lsoda(y=c(y1=8.586/0.438), times=TIME, func=PKde1, events=list(data=EventDat), parms=c(V=5.77348, CL=0.438, Kin=8.586))
y2 = lsoda(y=c(y1=21.93744), times=TIME, func=PKde2, events=list(data=EventDat), parms=c(V=5.94952, Km=507.873, Vmax=361.502, Kin=14.9684))

setdiff(y1[,"time"], dPK32[,"TIME"]) # 0, 125, 260
which(y1[,"time"] == 0)   # 1
which(y1[,"time"] == 125) # 28
which(y1[,"time"] == 260) # 52

cbind(y1[c(-1,-28,-52),], y2[c(-1,-28,-52),], dPK32)

# Figure 32.1, p 630
plot(dPK32[,"TIME"], dPK32[,"DV"], xlim=c(0, 450), ylim=c(0, 600), xlab="Time (min)", ylab="Concentration (ug/L)", pch=19, col="red")
lines(y1[,"time"], y1[,"y1"], lty=2)
lines(y2[,"time"], y2[,"y1"])
###

fPK32a = function(THETA)
{
  V   = THETA[1]
  CL  = THETA[2]
  Kin = THETA[3]

  EventDat = data.frame(var = rep("y1", nBolus),
                        time = c(0, 125, 260),
                        value = c(1669/V, 1701/V, 1733/V),
                        method = rep("add", nBolus)) ; EventDat

  y = lsoda(y=c(y1=Kin/V), times=TIME, func=PKde1, events=list(data=EventDat), parms=c(V=V, CL=CL, Kin=Kin))

  return(y[c(-1,-28,-52), "y1"])
}
fPK32a(c(5.77348, 0.438, 8.586))
length(fPK32a(c(5.77348, 0.438, 8.586)))
length(dPK32[,"DV"])


nlr(fPK32a, dPK32, pNames=c("V", "CL", "Kin"), IE=c(10, 1, 10), Error="P")


fPK32b = function(THETA)
{
  V    = THETA[1]
  Km   = THETA[2]
  Vmax = THETA[3]
  Kin  = THETA[4]

  EventDat = data.frame(var = rep("y1", nBolus),
                        time = c(0, 125, 260),
                        value = c(1669/V, 1701/V, 1733/V),
                        method = rep("add", nBolus)) ; EventDat

  y = lsoda(y=c(y1=Kin/V), times=TIME, func=PKde2, events=list(data=EventDat), parms=c(V=V, Km=Km, Vmax=Vmax, Kin=Kin))

  return(y[c(-1,-28,-52), "y1"])
}
fPK32b(c(5.94952, 507.873, 361.502, 14.9684))
length(fPK32a(c(5.94952, 507.873, 361.502, 14.9684)))
length(dPK32[,"DV"])


nlr(fPK32b, dPK32, pNames=c("V", "Km", "Vmax", "Kin"), IE=c(5, 450, 300, 15), Error="P")

