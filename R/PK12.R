setwd("D:/Rt/Gab/")
require(wnl)
dPK12 = read.csv("PK12.csv", skip=1)
colnames(dPK12) = c("TIME", "DV") ; dPK12

Dpo = 2500

require(deSolve)
PKde = function(t, y, p)
{
  if (t >= 60 & t < 75 ) {
    RateIn = 33.33333 # 500/15
  } else {
    RateIn = 0
  }

  if (t < p["Tlag"]) {
    Ka = 0
  } else {
    Ka = p["Ka"]
  }

  dy1dt = -Ka*y[1]                           # absorption, gut comp
  dy2dt = (RateIn + Ka*y[1] - p["Cl"]*y[2] - p["Cld"]*y[2] + p["Cld"]*y[3])/p["Vc"] # central, plasma
  dy3dt = (p["Cld"]*y[2] - p["Cld"]*y[3])/p["Vt"] # peripheral, tissue
  return(list(c(dy1dt, dy2dt, dy3dt)))
}

Times = c(0, dPK12[,"TIME"]) # simultaneous dosing <- not for this case
iTime = 2:length(Times)
y = lsoda(y=c(0.0464*Dpo, 0, 0), times=Times, func=PKde, parms=c(Ka=0.104, Vc=0.12, Vt=0.2758, Cl=0.014566, Cld=0.02, Tlag=4.8694))
plot(dPK12[,"TIME"], dPK12[,"DV"])
lines(y[,1], y[,3])

fPK12 = function(THETA)
{
  Fa   = THETA[1]
  Ka   = THETA[2]
  Vc   = THETA[3]
  Vt   = THETA[4]
  Cl   = THETA[5]
  Cld  = THETA[6]
  Tlag = THETA[7]

  Fs = lsoda(y=c(Fa*Dpo, 0, 0), times=Times, func=PKde, parms=c(Ka=Ka, Vc=Vc, Vt=Vt, Cl=Cl, Cld=Cld, Tlag=Tlag))
  return(Fs[iTime,"2"])
}

nlr(fPK12, dPK12, pNames=c("F", "Ka", "Vc", "Vt", "Cl", "Cld", "Tlag"), 
    IE=c(0.05, 0.1, 0.2, 0.6, 0.02, 0.01, 2),
    LB=c(0.01, 0.01, 0.02, 0.06, 0.002, 0.001, 0),
    UB=c(1, 1, 100, 1000, 10, 100, 5)) 

wnl5(fPK12, dPK12, pNames=c("F", "Ka", "Vc", "Vt", "Cl", "Cld", "Tlag"), IE=c(0.05, 0.1, 0.2, 0.6, 0.02, 0.01, 2)) 

