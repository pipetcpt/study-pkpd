setwd("D:/Rt/Gab/")
require(wnl)
dPK49 = read.csv("PK49.csv", as.is=TRUE)
colnames(dPK49) = c("TIME", "DV") # dPK49
dPK49 = dPK49[dPK49[,"DV"] != "Missing",] # dPK49
dPK49[,"DV"] = as.numeric(dPK49[,"DV"]) ; dPK49

AMT = 400
TAU = 19 / 60 # 19 min
infRate = 400/19*60 # 1263.158, 1261.83 in WinNonlin

fPK49 = function(THETA)
{
  V   = THETA[1]
  CL  = THETA[2]
  Syn = THETA[3]

  T1 = dPK49[,"TIME"][dPK49[,"TIME"] <= TAU] # only one point at TIME=0
  T2 = dPK49[,"TIME"][dPK49[,"TIME"] > TAU]
  y1 = Syn/CL + infRate/CL*(1 - exp(-CL/V*T1))
  y2 = Syn/CL + infRate/CL*(1 - exp(-CL/V*TAU))*exp(-CL/V*(T2 - TAU))
  
  return(c(y1, y2))
}
fPK49(c(3.59259, 0.14204, 16.2272))

nlr(fPK49, dPK49, pNames=c("V", "CL", "Syn"), IE=c(3.5, 0.12, 12))

# Figure 49.1, p 699
plot(dPK49[,"TIME"], dPK49[,"DV"], xlim=c(0, 160), ylim=c(80, 240), xlab="Time (h)", ylab="Concentration", pch=19)
lines(dPK49[,"TIME"], fPK49(e$PE[1:e$nTheta]))
###
