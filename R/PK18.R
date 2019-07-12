setwd("D:/Rt/Gab/")
require(wnl)
dPK18 = read.csv("PK18.csv", skip=1)
colnames(dPK18) = c("TIME", "DV") ; dPK18

Dinf = 28
Tau = 30
Rate = Dinf/Tau

require(deSolve)
PKde = function(t, y, p)
{
  if (t < Tau) {
    RateIn = Rate
  } else {
    RateIn = 0
  }
  Cl = p["Vmax"]/(p["Km"] + y[1])
  dy1dt = (RateIn - Cl*y[1] - p["Cld"]*y[1] + p["Cld"]*y[2])/p["Vc"] # central, plasma
  dy2dt = (p["Cld"]*y[1] - p["Cld"]*y[2])/p["Vt"] # peripheral, tissue
  return(list(c(dy1dt, dy2dt)))
}

Times = c(0, dPK18[, "TIME"])
iTime = 2:length(Times)

y = lsoda(y=c(0, 0), times=Times, func=PKde, parms=c(Vc=8.93, Vt=31.1174, Km=0.0125445, Vmax=0.0812189, Cld=1.29)) ; y
plot(dPK18[,"TIME"], dPK18[,"DV"], ylim=c(0.01, 10), log="y")
lines(y[,1], y[,2]) # Figure 18.2 p558

fPK18 = function(THETA)
{
  Vc   = THETA[1]
  Vt   = THETA[2]
  Km   = THETA[3]
  Vmax = THETA[4]
  Cld  = THETA[5]
  y = lsoda(y=c(0, 0), times=Times, func=PKde, parms=c(Vc=Vc, Vt=Vt, Km=Km, Vmax=Vmax, Cld=Cld))
  return(y[iTime,"1"])
}

nlr(fPK18, dPK18, pNames=c("Vc", "Vt", "Km", "Vmax", "Cld"), IE=c(20, 30, 0.03, 0.15, 1.5))
nlr(fPK18, dPK18, pNames=c("Vc", "Vt", "Km", "Vmax", "Cld"), IE=c(20, 30, 0.03, 0.15, 1.5), Error="POIS")
nlr(fPK18, dPK18, pNames=c("Vc", "Vt", "Km", "Vmax", "Cld"), IE=c(20, 30, 0.03, 0.15, 1.5), Error="P")
nlr(fPK18, dPK18, pNames=c("Vc", "Vt", "Km", "Vmax", "Cld"), IE=c(20, 30, 0.03, 0.15, 1.5), Error="C")
