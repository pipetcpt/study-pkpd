setwd("D:/Rt/Gab/")
require(wnl)
dPK39 = read.csv("PK39.csv", skip=1)
colnames(dPK39) = c("TIME", "DV", "DV2", "DV3") ; dPK39

require(deSolve)

infHist = cbind(TIME = c(0, 0.25, 8, 24),
                RATE = c(26.9/0.25, 139/7.75, 138.95/(24 - 8), 0)) ; infHist

PKde = function(t, y, p)
{
  RateIn = infHist[findInterval(t, infHist[,"TIME"]), "RATE"]
  Vc  = p["Vc"]
  Vt  = p["Vt"]
  CL  = p["CL"]
  CLd = p["CLd"]
  dy1dt = (RateIn - CL*y[1] - CLd*y[1] + CLd*y[2])/Vc
  dy2dt = (CLd*y[1] - CLd*y[2])/Vt
  return(list(c(dy1dt, dy2dt)))
}

cTIME = c(0, dPK39[,"TIME"])
iTime = 2:length(cTIME)

y = lsoda(y=c(0, 0), times=cTIME, func=PKde, parms=c(Vc=0.32, Vt=2.12265, CL=0.417793, CLd=0.903188))
y

# Figure 39.1, p 654
plot(dPK39[,"TIME"], dPK39[,"DV"], xlim=c(0,60), ylim=c(0, 60), xlab="Time (h)", ylab="Concentration (ug/L)", pch=19, col="red")
lines(cTIME[iTime], y[iTime,"1"])
###


fPK39 = function(THETA)
{
  Vc  = THETA[1]
  Vt  = THETA[2]
  CL  = THETA[3]
  CLd = THETA[4]
  
  y = lsoda(y=c(0, 0), times=cTIME, func=PKde, parms=c(Vc=Vc, Vt=Vt, CL=CL, CLd=CLd))
  return(y[iTime,"1"])
}
fPK39(c(0.32, 2.12265, 0.417793, 0.903188))

nlr(fPK39, dPK39, pNames=c("Vc", "Vt", "CL", "CLd"), IE=c(0.5, 2, 0.4, 1), Error="POIS") # Cov step failure with "P"


