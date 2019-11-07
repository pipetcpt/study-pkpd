setwd("D:/Rt/Gab/")
require(wnl)
dPK38 = read.csv("PK38.csv", skip=1)
colnames(dPK38) = c("TIME", "DV", "DOSE", "ID") ; dPK38

require(deSolve)
PKde = function(t, y, p)
{
  dy1dt = -p["Vm1"]*y[1]/(p["Km1"] + y[1]) - p["Vm2"]*y[1]/(p["Km2"] + y[1])
  return(list(c(dy1dt)))
}

DOSEs = c(50, 30, 10, 3, 1)
## Figure 38.1, p 651
plot(0, 0.001, type="n", xlim=c(0, 35), ylim=c(0.001, 100), log="y", xlab="Time (min)", ylab="Concentration")
IDs = unique(dPK38[,"ID"])
nID = length(IDs) ; nID

y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  cAMT = DOSEs[i]
  cTIME = c(0, dPK38[dPK38$ID == cID,"TIME"])
  iTime = 2:length(cTIME)
  cy = lsoda(y=c(cAMT), times=cTIME, func=PKde, parms=c(Vm1=0.96, Km1=0.0896813, Vm2=1.01877, Km2=8.68))
  points(dPK38[dPK38$ID == cID,"TIME"], dPK38[dPK38$ID == cID,"DV"], pch=14+i)
  lines(cTIME, cy[,"1"]) 
  y = c(y, cy[iTime,"1"])
} ; y
###

fPK38 = function(THETA)
{
  Vm1 = THETA[1]
  Km1 = THETA[2]
  Vm2 = THETA[3]
  Km2 = THETA[4]
  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = DOSEs[i]
    cTIME = c(0, dPK38[dPK38$ID == cID,"TIME"])
    iTime = 2:length(cTIME)
    cy = lsoda(y=c(cAMT), times=cTIME, func=PKde, parms=c(Vm1=Vm1, Km1=Km1, Vm2=Vm2, Km2=Km2))
    y = c(y, cy[iTime,"1"])
  }
  return(y)
}
fPK38(c(0.96, 0.0896813, 1.01877, 8.68))

nlr(fPK38, dPK38, pNames=c("Vmax1", "Km1", "Vmax2", "Km2"), IE=c(1, 0.1, 1, 10), Error="P") # Cov step failure
e$PE # 0.92267345 0.07970444 1.10895754 9.52993449 0.01344223


