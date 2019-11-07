setwd("D:/Rt/Gab/")
require(wnl)
dPK17 = read.csv("PK17.csv", skip=1)
colnames(dPK17) = c("TIME", "DV") ; dPK17

Tau1 = 0.5
Tau2 = 0.5 + 39.68
Rate1 = 1800/0.5 ; Rate1
Rate2 = 5484.8/39.68 ; Rate2 # Erratum p553: 39.63 min -> 39.68 min

require(deSolve)
PKde1 = function(t, y, p)
{
  if (t < Tau1) {
    RateIn = Rate1
  } else if (t < Tau2) {
    RateIn = Rate2
  } else {
    RateIn = 0
  }

  dy1dt = (RateIn - p["Cl"]*y[1])/p["V"] # central, plasma
  return(list(dy1dt))
}

PKde2 = function(t, y, p)
{
  if (t < Tau1) {
    RateIn = Rate1
  } else if (t < Tau2) {
    RateIn = Rate2
  } else {
    RateIn = 0
  }

  dy1dt = (RateIn - p["Vmax"]/(p["Km"] + y[1])*y[1])/p["V"] # central, plasma
  return(list(dy1dt))
}

## Figure 17.1
Times = dPK17[, "TIME"]
y = lsoda(y=c(0), times=Times, func=PKde1, parms=c(V=1380, Cl=43.3)) ; y
plot(dPK17[,"TIME"], dPK17[,"DV"])
lines(y[,1], y[,2])

y2 = lsoda(y=c(0), times=Times, func=PKde2, parms=c(V=1450, Vmax=107, Km=0.56)) ; y2
lines(y2[,1], y2[,2], lty=2)
##

fPK17a = function(THETA)
{
  V  = THETA[1]
  Cl = THETA[2]
  y = lsoda(y=c(0), times=Times, func=PKde1, parms=c(V=V, Cl=Cl))
  return(y[,"1"])
}

nlr(fPK17a, dPK17, pNames=c("V", "Cl"), IE=c(1500, 50), Error="P")

fPK17b = function(THETA)
{
  V    = THETA[1]
  Vmax = THETA[2]
  Km   = THETA[3]
  y = lsoda(y=c(0), times=Times, func=PKde2, parms=c(V=V, Vmax=Vmax, Km=Km))
  return(y[,"1"])
}

nlr(fPK17b, dPK17, pNames=c("V", "Vmax", "Km"), IE=c(1500, 120, 1), Error="P")

