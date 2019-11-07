setwd("D:/Rt/Gab/")
require(wnl)
dPK26 = read.csv("PK26.csv", skip=1)
colnames(dPK26) = c("TIME", "DV", "ID") ; dPK26
dPK26 = dPK26[!is.na(dPK26[,"DV"]),] ; dPK26

require(deSolve)
PKde = function(t, y, p)
{
  CLmm = p["Vmax"]/(p["Km"] + y[1]) # eq 26:2
  dy1dt = (-CLmm*y[1] - p["CLr"]*y[1] - p["CLd"]*y[1] + p["CLd"]*y[2])/p["Vc"]
  dy2dt = (p["CLd"]*y[1] - p["CLd"]*y[2])/p["Vt"] # eq 26:1

  return(list(c(dy1dt, dy2dt)))
}

# Figure 26.1, p600
plot(0, 0.01, type="n", xlim=c(0, 45), ylim=c(0.01, 1000), xlab="Time (days)", ylab="Concentration (mg/L)", log="y")
AMT = c(0.05, 0.3, 1, 3, 10)
IDs = unique(dPK26[,"ID"]) ; IDs
nID = length(IDs) ; nID
y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  cAMT = AMT[cID]
  cTime = c(0, dPK26[dPK26$ID == cID, "TIME"])
  iTime = 2:length(cTime)
  cy = lsoda(y=c(cAMT/0.0729, 0), times=cTime, func=PKde, parms=c(Vc=0.0729, Vmax=0.0338, Km=0.076, Vt=0.029321, CLd=0.007, CLr=0.007))
  y = c(y, cy[iTime,"1"])
  points(dPK26[dPK26$ID == cID, "TIME"], dPK26[dPK26$ID == cID, "DV"], pch=19)
  lines(dPK26[dPK26$ID == cID, "TIME"], cy[iTime,"1"])
} ; y
####

fPK26 = function(THETA)
{
  Vc   = THETA[1]
  Vmax = THETA[2]
  Km   = THETA[3]
  Vt   = THETA[4]
  CLd  = THETA[5]
  CLr  = THETA[6]

  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMT[cID]
    cTime = c(0, dPK26[dPK26$ID == cID, "TIME"])
    iTime = 2:length(cTime)
    cy = lsoda(y=c(cAMT/Vc, 0), times=cTime, func=PKde, parms=c(Vc=Vc, Vmax=Vmax, Km=Km, Vt=Vt, CLd=CLd, CLr=CLr))
    y = c(y, cy[iTime,"1"])
  } 
  return(y)
}
fPK26(c(0.0729, 0.0338, 0.076, 0.029321, 0.007, 0.007))

nlr(fPK26, dPK26, pNames=c("Vc", "Vmax", "Km", "Vt", "CLd", "CLr"), IE=c(0.1, 0.1, 0.1, 0.1, 0.01, 0.01), Error="P")
# Table 26.1 seems wrong <- compare with WinNonlin and R results

