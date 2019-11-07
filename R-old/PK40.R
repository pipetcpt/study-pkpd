setwd("D:/Rt/Gab/")
require(wnl)
dPK40 = read.csv("PK40.csv", skip=1)
colnames(dPK40) = c("TIME", "DV") ; dPK40

Bolus = 5617.3

require(deSolve)
PKde = function(t, y, p)
{
  if (t >= 10 & t < (10 + p["Tau"])) {
    Rate = 1/p["Tau"]
  } else {
    Rate = 0
  }

	dy1dt = Rate*y[2] - p["Ka"]*y[1]          # eq 40:5, gut
	dy2dt = p["K1g"]*y[3]*p["Vc"] - Rate*y[2] # eq 40:4, bile
	dy3dt = (p["Ka"]*y[1] - p["CL"]*y[3] - p["CLd"]*y[3] + p["CLd"]*y[4])/p["Vc"] - p["K1g"]*y[3] # eq 40:1
	dy4dt = (p["CLd"]*y[3] - p["CLd"]*y[4])/p["Vt"]  # eq 40:2, peripheral
	return(list(c(dy1dt, dy2dt, dy3dt, dy4dt)))
}

cTIME = c(0, dPK40[,"TIME"])
iTime = 2:length(cTIME)
y = lsoda(y=c(0, 0, Bolus/12.82, 0), times=cTIME, func=PKde, parms=c(Vc=12.82, Vt=29, CL=0.8421, CLd=11.37, Ka=3, K1g=0.6, Tau=2.8))
y

# Figure 40.1, p 658
plot(dPK40[,"TIME"], dPK40[,"DV"], xlim=c(0, 40), ylim=c(1, 1000), log="y", xlab="Time (h)", ylab="Concentration (uM)", pch=19, col="red")
lines(cTIME[iTime], y[iTime,"3"])
###

fPK40 = function(THETA)
{
  Vc  = THETA[1]
  Vt  = THETA[2]
  CL  = THETA[3]
  CLd = THETA[4]
  Ka  = THETA[5]
  K1g = THETA[6]
  Tau = THETA[7]

  y = lsoda(y=c(0, 0, Bolus/Vc, 0), times=cTIME, func=PKde, parms=c(Vc=Vc, Vt=Vt, CL=CL, CLd=CLd, Ka=Ka, K1g=K1g, Tau=Tau))
  return(y[iTime,"3"])
}
fPK40(c(12.82, 29, 0.8421, 11.37, 3, 0.6, 2.8))

nlr(fPK40, dPK40, pNames=c("Vc", "Vt", "CL", "CLd", "Ka", "K1g", "Tau"), IE=c(13, 35, 1.2, 10, 2.3, 0.7, 2.8), Error="P")



