setwd("D:/Rt/Gab/")
require(wnl)
dPK31 = read.csv("PK31.csv", skip=1)
colnames(dPK31) = c("TIME", "DV", "R0", "EXO") ; dPK31

require(deSolve)
PKde = function(t, y, p)
{
  if (t < 1/60) {
    RateIn = 36630*60
  } else {
    RateIn = 0
  }
  dy1dt = (RateIn + p["Syn"] - p["CL"]*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"]
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"]
  
  return(list(c(dy1dt, dy2dt)))
}

TIME = c(0, dPK31[,"TIME"])
iTime = 2:length(TIME)

y = lsoda(y=c(1531.87/76.6, 1531.87/76.6), times=TIME, func=PKde, parms=c(Syn=1531.87, Vc=8.8455, Vt=58.8, CL=76.6, CLd=56.8775))
y

fPK31 = function(THETA)
{
  Syn = THETA[1]
  Vc  = THETA[2]
  Vt  = THETA[3]
  CL  = THETA[4]
  CLd = THETA[5]
  
  y = lsoda(y=c(Syn/CL, Syn/CL), times=TIME, func=PKde, parms=c(Syn=Syn, Vc=Vc, Vt=Vt, CL=CL, CLd=CLd))
  return(y[iTime,"1"])
}
fPK31(c(1531.87, 8.8455, 58.8, 76.6, 56.8775))

nlr(fPK31, dPK31, pNames=c("Syn", "Vc", "Vt", "CL", "CLd"), IE=c(2000, 10, 100, 100, 100), Error="P")


