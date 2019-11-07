setwd("D:/Rt/Gab/")
require(wnl)
dPK42 = read.csv("PK42.csv", skip=1, as.is=TRUE)
colnames(dPK42) = c("TIME", "DV", "ID", "DOSE") ; dPK42
dPK42 = dPK42[!is.na(dPK42$DV),] ; dPK42

# oral 1 comp, MM absorption
AMTs = c(10000, 30000, 90000)

require(deSolve)
PKde = function(t, y, p)
{
	dy1dt = - p["Vm"]*y[1]/(p["Km"] + y[1])
	dy2dt = (p["Vm"]* y[1]/(p["Km"] + y[1]) - p["CL"]*y[2] - p["CLd"]*y[2] + p["CLd"]*y[3])/p["Vc"]
	dy3dt = (p["CLd"]*y[2] - p["CLd"]*y[3])/p["Vt"]
  return(list(c(dy1dt, dy2dt, dy3dt)))	
}

# Figure 42.1, p 666
plot(0, 0.01, xlim=c(0, 400), ylim=c(0.1, 1000), log="y", xlab="Time (min)", ylab="Concentration (ug/L)")
Bullets = c(19, 15, 17)
 
IDs = unique(dPK42$ID)
nID = length(IDs)

y = vector()
for (i in 1:nID) {
  cID = IDs[i]
  cAMT = AMTs[i]
  cTIME = c(0, dPK42[dPK42$ID == cID, "TIME"])
  iTime = 2:length(cTIME)
  cy = lsoda(y=c(cAMT, 0, 0), times=cTIME, func=PKde, parms=c(Vm=982.453, Km=9570.63, Vc=4.66257, Vt=34.9584, CL=2, CLd=0.986))
  points(dPK42[dPK42$ID == cID, "TIME"], dPK42[dPK42$ID == cID, "DV"], pch=Bullets[i])
  lines(cTIME[iTime], cy[iTime, "2"])
  y = c(y, cy[iTime,"2"]) 
} ; y
###

fPK42 = function(THETA)
{
  Vm  = THETA[1]
  Km  = THETA[2]
  Vc  = THETA[3]
  Vt  = THETA[4]
  CL  = THETA[5]
  CLd = THETA[6]

  y = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    cAMT = AMTs[i]
    cTIME = c(0, dPK42[dPK42$ID == cID, "TIME"])
    iTime = 2:length(cTIME)
    cy = lsoda(y=c(cAMT, 0, 0), times=cTIME, func=PKde, parms=c(Vm=Vm, Km=Km, Vc=Vc, Vt=Vt, CL=CL, CLd=CLd))  
    y = c(y, cy[iTime,"2"]) 
  }
  
  return(y)
}
fPK42(c(982.453, 9570.63, 4.66257, 34.9584, 2, 0.986))

nlr(fPK42, dPK42, pNames=c("Vm", "Km", "Vc", "Vt", "CL", "CLd"), IE=c(900, 11000, 5.2, 40, 2.2, 1.1), Error="P")

